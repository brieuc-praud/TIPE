import numpy as np


RT = 6311000 #rayon de la Terre (en m)
latitude_inter = 44#latitude intermédiaire
RTL = RT*np.cos(latitude_inter*np.pi/180)

def produit(ARRAY):#produit matriciel
    P = np.identity(4)
    for M in ARRAY:
        P = np.dot(P,M)
    return P

class Boat:
    delta_t_pred = 60*600
    class Kalman:
        def __init__(self, Q, R):
            #constantes
            self._Q = Q
            self._R = R
            #variables
            self._X = []#vecteur_kalman = [x, y, vx, vy]
            self._P = np.identity(4)#matrice de covariance de l'état estimé, arbitrairement grande au départ

        def predire(self, delta_t):

            X = self._X
            P = self._P
            Q = self._Q
            F = np.array([[1,0,delta_t,0],[0,1,0,delta_t],[0,0,1,0],[0,0,0,1]]) #matrice représentant le modèle physique

            X_prime = produit((F, X))
            P_prime = produit((F, P, F.T)) + Q
            return X_prime, P_prime

        def filtre(self, vecteurs):
            if self._X==[]:#pour le premier appel de la fonction
                self._X = vecteurs[-1][1:]

            #"prédiction"
            t1 = vecteurs[0][0]
            t2 = vecteurs[1][0]
            delta_t = t2-t1
            X_prime, P_prime = self.predire(delta_t)

            #"mise à jour"
            R = self._R
            I = np.identity(4)
            Y = vecteurs[1][1:]#mesures reçues

            K = produit((P_prime, np.linalg.inv(P_prime + R)))#gain de Kalman optimal

            self._X = X_prime + produit((K, Y - X_prime))
            self._P = produit((I-K, P_prime, I-K.T)) + produit((K, R, K.T))
        def recuperer_etat(self):
            return self._X, self._P

    def __init__(self, mmsi, echelle, longitude_min, longitude_max, latitude_min, latitude_max):#vecteur = [instant, latitude, longitude, vitesse_latitudinale, vitesse_longitudinale]
        Q = 5e10*np.array([[1e-10,0,0,0],[0,1e-10,0,0],[0,0,1e-8,0],[0,0,0,1e-8]])#Q
        R = 5e10*np.array([[1e-10,0,0,0],[0,1e-10,0,0],[0,0,1e-15,0],[0,0,0,1e-15]])#R

        self._mmsi = mmsi
        self._longitude_min, self._longitude_max = longitude_min, longitude_max
        self._latitude_min, self._latitude_max = latitude_min, latitude_max

        self._echelle = echelle
        self._r_risque = 0#rayon sur lequel on devra calculer le risque autour d'une position donnée

        self._kalman = Boat.Kalman(Q, R)
        self._vecteurs =[]#cette liste contiendra les deux derniers vecteurs ajoutés, ce sera utile pour le filtre.
        self._Dernier_risque = []#[mini, maxi, MAT]
    def _maj_risque(self, MAT_risque, kalman):
        def convertir(x, y):#traduit (x, y) en coordonnées dans la matrice de risque
            xm = round(x*self._echelle)
            ym = round(y*self._echelle)
            return xm, ym
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
        X_prime, P_prime = kalman.predire(Boat.delta_t_pred)#estime la position du bateau dans Boat.delta_t_pred secondes

        posx, posy = convertir(X[0], X[1])

        kalx, kaly = convertir(X_prime[0], X_prime[1])
        xmin, xmax, ymin, ymax = min(posx, kalx), max(posx, kalx), min(posy, kaly), max(posy, kaly)
        mini, maxi = min(xmin, ymin), max(xmax, ymax)

        x, y = kalx-posx, kaly-posy
        r = np.sqrt(x**2 + y**2)
        if not x+r == 0:
             theta = 2*np.arctan(y/(x+r))
        else:
             theta = 0
        sig = r/4
        #print("var="+str(var))

        varx, vary, varx_prime, vary_prime = P[1,1], P[0,0], P_prime[1,1], P_prime[0,0]
        sigx, sigy, sigx_prime, sigy_prime = np.sqrt(varx), np.sqrt(vary), np.sqrt(varx_prime), np.sqrt(vary_prime)
        sigu = np.cos(theta)*sigx + np.sin(theta)*sigy
        sigv = -np.sin(theta)*sigx + np.cos(theta)*sigy
        sigu_prime = np.cos(theta)*sigx_prime + np.sin(theta)*sigy_prime
        sigv_prime = -np.sin(theta)*sigx_prime + np.cos(theta)*sigy_prime

        R = self._r_risque
        U, V = changer_repere(maxi-mini+2*R, r, theta)
        # #print("varx="+str(varx_prime))
        G = Risque(U, V, sigu+sig, sigv, sigu_prime, sigv_prime)

        D = self._Dernier_risque
        if D:#si il y a le précedent risque de ce bateau à effacer, alors on l'efface
            MAT_risque[D[0]-R:D[1]+R, D[0]-R:D[1]+R] = MAT_risque[D[0]-R:D[1]+R, D[0]-R:D[1]+R] - D[2]

        MAT_risque[mini-R:maxi+R, mini-R:maxi+R] = MAT_risque[mini-R:maxi+R, mini-R:maxi+R] + G
        print(G)
        self._Dernier_risque=[mini, maxi, G]
    def append(self, infos, MAT_risque):#méthode appelée pour ajouter une mesure au bateau

        instant = infos[0]
        longitude = float(infos[2])
        latitude = float(infos[1])*np.pi/180
        SOG = float(infos[3])*0.514 #1 noeud ~= 0.514 m/s
        COG = float(infos[4])*np.pi/180 #1 rad = pi/180

        x = (longitude - self._longitude_min)*np.pi/180*RTL
        y = latitude*np.pi/180*RT

        vitx, vity = SOG*np.sin(COG), SOG*np.cos(COG)

        vecteur = [instant, x, y, vitx, vity]

        if len(self._vecteurs) == 0:
            self._vecteurs.append(vecteur)
        else:
            if len(self._vecteurs) == 2:#kalman utilise les instants t et t+1
                self._vecteurs[0], self._vecteurs[1] = self._vecteurs[1], vecteur
            else:
                self._vecteurs.append(vecteur)

            kalman = self._kalman
            kalman.filtre(self._vecteurs)
            self._maj_risque(MAT_risque, kalman)



