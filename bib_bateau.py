import numpy as np
import multiprocessing as mp

import bib_dijkstra as dij

def produit(*Matrices):#produit matriciel
    P = np.identity(4)
    for M in Matrices:
        P = np.dot(P,M)
    return P

class Embarcation:
    def __init__(self, echelle, longitude_min, latitude_min):
        self._longitude_min = longitude_min
        self._latitude_min = latitude_min
        self._echelle = echelle
    def _convertir(self, longitude, latitude):#traduit (longitude, latitude) en coordonnées dans la matrice de risque
        xm = round((longitude - self._longitude_min)*self._echelle)
        ym = round((latitude - self._latitude_min)*self._echelle)
        return xm, ym


class Bateau(Embarcation):
    class Kalman:
        def __init__(self, Q, R, vecteur_kalman):#vecteur_kalman = [latitude, longitude, vitesse_latitudinale, vitesse_longitudinale]
            #constantes
            self._Q = Q
            self._R = R
            #variables
            self._X = vecteur_kalman
            self._P = np.identity(4)#matrice de covariance de l'état estimé, arbitrairement grande au départ
        def _predire(self, delta_t):#estime la position du bateau dans delta_t_pred secondes
            X = self._X
            P = self._P
            Q = self._Q
            F = np.array([[1,0,delta_t,0],[0,1,0,delta_t],[0,0,1,0],[0,0,0,1]]) #matrice représentant le modèle physique

            X_prime = produit(F, X)
            P_prime = produit(F, P, F.T) + Q
            return X_prime, P_prime
        def filtre(self, vecteurs):
            #"prédiction"
            t1 = vecteurs[0][0]
            t2 = vecteurs[1][0]
            delta_t = t2-t1
            X_prime, P_prime = self._predire(delta_t)

            #"mise à jour"
            R = self._R
            I = np.identity(4)
            Y = vecteurs[1][1:]#mesures reçues

            K = produit(P_prime, np.linalg.inv(P_prime + R))#gain de Kalman optimal

            self._X = X_prime + produit(K, Y - X_prime)
            self._P = produit(I-K, P_prime, I-K.T) + produit(K, R, K.T)
        def recuperer_etat(self):
            return self._X, self._P

    def __init__(self, mmsi, vecteur, echelle, longitude_min, latitude_min):#vecteur = [instant, latitude, longitude, vitesse_latitudinale, vitesse_longitudinale]
        super().__init__(echelle, longitude_min, latitude_min)
        Q = 1e3*np.array([[1e-6,0,0,0],[0,1e-6,0,0],[0,0,1e-11,0],[0,0,0,1e-11]])#Q
        R = np.array([[1e-6,0,0,0],[0,1e-6,0,0],[0,0,1e-11,0],[0,0,0,1e-11]])#R

        self._mmsi = mmsi
        self._r_risque = 50#rayon sur lequel on devra calculer le risque autour d'une position donnée

        self._kalman = Bateau.Kalman(Q, R, vecteur[1:])
        self._vecteurs =[vecteur]#cette liste contiendra les deux derniers vecteurs ajoutés, ce sera utile pour le filtre.
        self._Dernier_risque = []#[mini, maxi, MAT]
    def _calculer_delta_t(self, mon_bateau):
        vitesse = 1e-6
        vecteur = self._vecteurs[-1]
        position1 = mon_bateau.recuperer_position()
        lon0, lat0 = vecteur[1], vecteur[2]
        lon1, lat1 = position1[0], position1[1]
        distance = np.sqrt((lon1-lon0)**2 + (lat1-lat0)**2)#distance euclidienne
        delta_t = round(distance/vitesse)
        return delta_t
    def _maj_risque(self, MAT_risque, kalman, mon_bateau):
        def changer_repere(res, r, theta):
            #création du repere
            x = y = np.linspace(-res, res, num=res)
            X, Y = np.meshgrid(x, y)
            #rotation
            U, V = X*np.cos(theta) - Y*np.sin(theta), X*np.sin(theta) + Y*np.cos(theta)
            #translation
            U = U + r/2
            return U, V
        def Risque(X, Y, sigx, sigy, sigx_prime, sigy_prime):#crée une matrice avec une répartition normale
            T = np.linspace(0,1,len(X))
            Sx = (1-T)*sigx + T*sigx_prime
            Sy = (1-T)*sigy + T*sigy_prime
            return 100*np.exp(-X**2/(2*Sx**2) - Y**2/(2*Sy**2))

        X, P = kalman.recuperer_etat()

        delta_t_pred = self._calculer_delta_t(mon_bateau)
        X_prime, P_prime = kalman._predire(delta_t_pred)

        posx, posy = self._convertir(X[0], X[1])
        kalx, kaly = self._convertir(X_prime[0], X_prime[1])
        xmin, xmax, ymin, ymax = min(posx, kalx), max(posx, kalx), min(posy, kaly), max(posy, kaly)
        mini, maxi = min(xmin, ymin), max(xmax, ymax)

        x, y = kalx-posx, kaly-posy
        r = np.sqrt(x**2 + y**2)
        if not x+r == 0:
             theta = 2*np.arctan(y/(x+r))
        else:
             theta = 0
        sig = r/4

        varx, vary, varx_prime, vary_prime = P[0,0], P[1,1], P_prime[0,0], P_prime[1,1]
        sigx, sigy, sigx_prime, sigy_prime = np.sqrt(varx), np.sqrt(vary), np.sqrt(varx_prime), np.sqrt(vary_prime)
        sigu = np.cos(theta)*sigx + np.sin(theta)*sigy
        sigv = -np.sin(theta)*sigx + np.cos(theta)*sigy
        sigu_prime = np.cos(theta)*sigx_prime + np.sin(theta)*sigy_prime
        sigv_prime = -np.sin(theta)*sigx_prime + np.cos(theta)*sigy_prime

        R = self._r_risque
        U, V = changer_repere(maxi-mini+2*R, r, theta)

        G = Risque(U, V, sigu+sig, sigv, sigu_prime, sigv_prime)

        def tronquer(G, min, max, taille, taille_tot):
            taille = taille[0]
            if min < 0:
                G = G[-taille, -taille]
            elif max >= taille_tot:
                G = G[taille, taille]
            return G

        G = tronquer(G, mini-R, maxi+R, MAT_risque[mini-R:maxi+R, mini-R:maxi+R].shape, MAT_risque.shape[0])
        D = self._Dernier_risque
        if D:#si il y a le précedent risque de ce bateau à effacer, alors
            MAT_risque[D[0]-R:D[1]+R, D[0]-R:D[1]+R] = MAT_risque[D[0]-R:D[1]+R, D[0]-R:D[1]+R] - D[2]

        MAT_risque[mini-R:maxi+R, mini-R:maxi+R] = MAT_risque[mini-R:maxi+R, mini-R:maxi+R] + G
        self._Dernier_risque=[mini, maxi, G]


    def append(self, vecteur, MAT_risque, mon_bateau):#méthode appelée pour ajouter une mesure au bateau
        if len(self._vecteurs) == 2:#kalman utilise les instants t et t+1
            self._vecteurs[0], self._vecteurs[1] = self._vecteurs[1], vecteur
        else:#executé au premier appel de append()
            self._vecteurs.append(vecteur)

        kalman = self._kalman
        kalman.filtre(self._vecteurs)
        self._maj_risque(MAT_risque, kalman, mon_bateau)

class MonBateau(Embarcation):
    def __init__(self, n, echelle, longitude_min, latitude_min):
        super().__init__(echelle, longitude_min, latitude_min)

        self._Graphe = dij.graphe(n) #Le graphe pour Dijkstra
        self._depart, self._arrivee = (-68,43), (-66,45)

        manager = mp.Manager()
        self._return_dict = manager.dict()
        self._return_dict["LON"] = self._return_dict["LAT"] = self._return_dict["chemin"] = []
        self._dij_proc = mp.Process(target=None)
    def _dijkstra(self, M):#trouve le plus court chemin dans la matrice M
        depart, arrivee = self._convertir(self._depart[0], self._depart[1]), self._convertir(self._arrivee[0], self._arrivee[1])
        chemin = dij.trajet(self._Graphe, M, depart, arrivee)[1]
        LON, LAT = dij.preparer_trajet(chemin, self._echelle, self._longitude_min, self._latitude_min)
        self._return_dict["LON"], self._return_dict["LAT"], self._return_dict["chemin"] = LON, LAT, chemin
    def recuperer_position(self):
        return self._depart
    def calculer_plus_court_chemin(self, M):
        if not self._dij_proc.is_alive():#On relance le calcul seulement si le précédent est fini
            self._dij_proc = mp.Process(target=self._dijkstra, args=[M])
            self._dij_proc.start()
            chemin = self._return_dict["chemin"]
            #if len(chemin):
            #    dep = (chemin[vitesse_bateau_suivi][0], chemin[vitesse_bateau_suivi][1])#on fait avancer le bateau suivi
            dij.affiche_trajet(self._return_dict["LON"], self._return_dict["LAT"])
