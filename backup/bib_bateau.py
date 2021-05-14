import numpy as np

def produit(ARRAY):#produit matriciel
    P = np.identity(4)
    for M in ARRAY:
        P = np.dot(P,M)
    return P

class Boat:
    delta_t_pred = 60*60
    class Kalman:
        def __init__(self, Q, R, vecteur_kalman):#vecteur_kalman = [latitude, longitude, vitesse_latitudinale, vitesse_longitudinale]
            #constantes
            self._Q = Q
            self._R = R
            #variables
            self._X = vecteur_kalman
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

    def __init__(self, mmsi, vecteur, echelle, latitude_min, longitude_min):#vecteur = [instant, latitude, longitude, vitesse_latitudinale, vitesse_longitudinale]
        Q = 5e10*np.array([[1e-10,0,0,0],[0,1e-10,0,0],[0,0,1e-8,0],[0,0,0,1e-8]])#Q
        R = 5e10*np.array([[1e-10,0,0,0],[0,1e-10,0,0],[0,0,1e-15,0],[0,0,0,1e-15]])#R

        self._mmsi = mmsi
        self._longitude_min = longitude_min
        self._latitude_min = latitude_min
        self._echelle = echelle
        self._r_risque = 20#rayon sur lequel on devra calculer le risque autour d'une position donnée

        self._kalman = Boat.Kalman(Q, R, vecteur[1:])
        self._vecteurs =[vecteur]#cette liste contiendra les deux derniers vecteurs ajoutés, ce sera utile pour le filtre.
        self._Dernier_risque = []#[mini, maxi, MAT]
    def _maj_risque(self, MAT_risque, kalman):
        def convertir(longitude, latitude):#traduit (longitude, latitude) en coordonnées dans la matrice de risque
            x = round((longitude - self._longitude_min)*self._echelle)
            y = round((latitude - self._latitude_min)*self._echelle)
            return x, y
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

        posx, posy = convertir(X[1], X[0])
        kalx, kaly = convertir(X_prime[1], X_prime[0])
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
        #print("varx="+str(varx_prime))
        G = Risque(U, V, sigu+sig, sigv, sigu_prime, sigv_prime)

        D = self._Dernier_risque
        if D:#si il y a le précedent risque de ce bateau à effacer, alors
            MAT_risque[D[0]-R:D[1]+R, D[0]-R:D[1]+R] = MAT_risque[D[0]-R:D[1]+R, D[0]-R:D[1]+R] - D[2]

        MAT_risque[mini-R:maxi+R, mini-R:maxi+R] = MAT_risque[mini-R:maxi+R, mini-R:maxi+R] + G
        self._Dernier_risque=[mini, maxi, G]
    def append(self, vecteur, MAT_risque):#méthode appelée pour ajouter une mesure au bateau
        if len(self._vecteurs) == 2:#kalman utilise les instants t et t+1
            self._vecteurs[0], self._vecteurs[1] = self._vecteurs[1], vecteur
        else:#executé au premier appel de append()
            self._vecteurs.append(vecteur)

        kalman = self._kalman
        kalman.filtre(self._vecteurs)
        self._maj_risque(MAT_risque, kalman)

#
# class Boat_backup:#bateau de base
#     kal_Q = 1e-6*np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) #matrice de covariance du bruit du modèle physique
#     kal_R = 1e-7*np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) #matricee de covariance liée aux bruits des capteurs (donné par le constructeur du capteur)
#     delta_t_prediction=300#temps pour lequel doit s'effecteur la prédiction
#     def __init__(self, mmsi, vecteur):#vecteur = [instant, latitude, longitude, vitesse_latitudinale, vitesse_longitudinale]
#         self._mmsi = mmsi
#         self._liste_vecteurs = []#avec cette liste, on garde l'historique des positions de chaque bateau
#         self._liste_vecteurs_kal = []#liste des positions successives enregistrées par le filtre
#         self._liste_vecteurs_kal_pred = []#liste des prédictions successives effectuées par le filtre
#         self._liste_kal_cov = []#liste des covariances successives du filtre
#         self._liste_kal_cov_pred = []#liste des covariances successives des predictions
#
#         #pour le filtre de Kalman
#         self._kal_vecteur=vecteur[1:] #[latitude, longitude, vitesse_latitudinale, vitesse_longitudinale]
#         self._kal_P = np.identity(4) #matrice de covariance de l'état estimé, arbitrairement grande au départ
#     def _prediction(self, delta_t):#la phase prédiction du filtre de Kalman
#         F = np.array([[1,0,delta_t,0],[0,1,0,delta_t],[0,0,1,0],[0,0,0,1]]) #matrice représentant le modèle physique
#         kal_vecteur_prime = produit((F, self._kal_vecteur))
#         kal_P_prime = produit((F, self._kal_P, F.T)) + Boat.kal_Q
#         return kal_vecteur_prime, kal_P_prime
#     def _kalman(self):#le filtre
#         #"__prediction"
#         t1 = self._liste_vecteurs[-2][0]
#         t2 = self._liste_vecteurs[-1][0]
#         delta_t = t2-t1
#
#         kal_vecteur_prime, kal_P_prime = self._prediction(delta_t)
#
#         #"mise a jour"
#         K=np.dot(kal_P_prime, np.linalg.inv(kal_P_prime + Boat.kal_R))#gain de Kalman optimal
#
#         self._kal_vecteur = kal_vecteur_prime + produit((K, self._liste_vecteurs[-1][1:]-kal_vecteur_prime))
#         self._kal_P = produit((np.identity(4)-K, kal_P_prime, np.identity(4)-K.T)) + produit((K, Boat.kal_R, np.transpose(K)))
#
#         self._liste_vecteurs_kal.append(self._kal_vecteur)
#         self._liste_kal_cov.append(self._kal_P)
#     def _predire(self, t):#prédit et affiche la position du bateau au bout de t secondes
#         kal_vecteur_prime, kal_P_prime = self._prediction(t)
#         self._liste_vecteurs_kal_pred.append(kal_vecteur_prime)
#         self._liste_kal_cov_pred.append(kal_P_prime)
#         return kal_vecteur_prime, kal_P_prime
#     def append(self, vecteur):#met à jour le bateau avec le vecteur passé en paramètre
#         self._liste_vecteurs.append(vecteur)
#         if len(self._liste_vecteurs) >= 2:#kalman utilise l'instant t+1
#             self._kalman()#calcule et affiche la position estimée du bateau à l'instant même
#             self._predire(Boat.delta_t_prediction)#prédit la position du bateau
#
#
# #redondances de code avec la classe Boat_graph, propreté à améliorer
# class Boat_risque(Boat):#bateau qui met a jour la carte des risques.
#     R=10
#     delta_t_prediction_L = [x*10 for x in range(0,50)]
#     def __init__(self, mmsi, vecteur, axes):#vecteur = [instant, latitude, longitude, vitesse_latitudinale, vitesse_longitudinale]
#
#         Boat.__init__(self, mmsi, vecteur)
#
#         self._dot, = axes.plot([],[], marker='o',color='red',markersize=0) #le point représentant la position fournie pas l'AIS sur la carte
#         self._kal_dot, = axes.plot([],[], marker='x',color='green',markersize=0) #le point représentant l'estimation du filtre de Kalman
#
#
#     def _prediction(self, delta_t):#la phase prédiction du filtre de Kalman
#         F = np.array([[1,0,delta_t,0],[0,1,0,delta_t],[0,0,1,0],[0,0,0,1]]) #matrice représentant le modèle physique
#         kal_vecteur_prime = produit((F, self._kal_vecteur))
#         kal_P_prime = produit((F, self._kal_P, F.T)) + Boat.kal_Q
#         return kal_vecteur_prime, kal_P_prime
#     def _kalman(self):#le filtre
#         Boat._kalman(self)
#         self._kal_dot.set_data(self._kal_vecteur[1], self._kal_vecteur[0])
#     def append(self, vecteur, M, scale, lat_min, lon_min):#met à jour le bateau avec le vecteur passé en paramètre
#         Boat.append(self, vecteur)#self.__liste_vecteurs.append(vecteur)
#         self._dot.set_data(vecteur[2], vecteur[1])#modifie la position AIS affichée
#         if len(self._liste_vecteurs) >= 2:#kalman utilise l'instant t+1
#             self._kalman()#calcule et affiche la position estimée du bateau à l'instant même
#             predictions=[]
#             for delta_t in Boat_risque.delta_t_prediction_L:
#                 predictions.append(self._predire(delta_t))#prédit la position du bateau
#             self._risque_update(M, scale, predictions, lat_min, lon_min)
#     def _risque_update(self, M, scale, predictions, lat_min, lon_min):
#         R=Boat_risque.R
#
#         for pred in predictions:
#
#             vect = pred[0]
#             COV = pred[1]
#
#             varx, vary = COV[0,0], COV[1,1]
#
#             x=np.linspace(-R,R, num=2*R)
#             y=x
#             x, y = np.meshgrid(x, y)
#
#             posx, posy = convertir(vect[1], vect[0], lon_min, lat_min, scale)
#             MAT = gaussienne(x, y, 1000, 5*10e0*varx, 5*10e0*vary)
#             M[posy-R:posy+R, posx-R:posx+R] = M[posy-R:posy+R, posx-R:posx+R] + MAT
#
#
