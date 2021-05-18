##bibliothèques
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs#à installer avec "pip install Cartopy-0.18.0-cp38-cp38-win_amd64.whl" (problèmes de dépendances si installée depuis les dépots de Python)


##localisation de l'étude
n, echelle = 8, 150#la largeur de la carte est d'environ n*echelle km
latitude_min = 42.5-3
latitude_max = latitude_min + n
longitude_min = -69.5-3
longitude_max = longitude_min + n
n = n*echelle

##couleurs
colormap='jet'#couleurs utilisées pour représenter le risque
norme = LogNorm()#None


fig = plt.figure()
axes = plt.axes()
axes.set_aspect('equal') #évite la déformation de la carte
plt.xlim(longitude_min, longitude_max)
plt.ylim(latitude_min, latitude_max)

def changer_repere(res, r, theta):
    #création du repere
    x = y = np.linspace(-res, res, num=res)
    X, Y = np.meshgrid(x, y)
    #rotation
    U, V = X*np.cos(theta) - Y*np.sin(theta), X*np.sin(theta) + Y*np.cos(theta)
    #translation
    U = U - r/2
    return U, V
def Risque(X, Y, sigx, sigy, sigx_prime, sigy_prime):#crée une matrice avec une répartition normale
    #T = np.linspace(0,1,len(X))
    xmax = X.max()
    T = X/(2*xmax)+0.5
    Sx = (1-T)*sigx + T*sigx_prime
    Sy = (1-T)*sigy + T*sigy_prime
    return 100*np.exp(-X**2/(2*Sx**2) - Y**2/(2*Sy**2))

M = N = O = P = Q = np.ones((n,n))
im = axes.imshow(M, extent=[longitude_min, longitude_max, latitude_min, latitude_max], cmap='seismic')#dessine la carte de risque

U,V = changer_repere(n, 1000/4, 0)
M += Risque(U, V, 1e1+1000/4, 1e1, 1e2, 1e2)
U,V = changer_repere(n, 720/4, -2*np.pi/3)
U, V = U, V+800
N += Risque(U, V, 1e1+360/4, 1e1, 1e2, 1e2)
U,V = changer_repere(n, 0, 0)
U, V = U-200, V+490
O += Risque(U, V, 3e1, 3e1, 3e1, 3e1)
U,V = changer_repere(n, 0, 0)
U, V = U-150, V+500
P += Risque(U, V, 3e1, 3e1, 3e1, 3e1)
U,V = changer_repere(n, 0, 0)
U, V = U, V+400
Q += Risque(U, V, 3e1, 3e1, 3e1, 3e1)


im = axes.imshow(M, extent=[longitude_min, longitude_max, latitude_min, latitude_max], cmap=colormap, norm=norme)
im = axes.imshow(N, extent=[longitude_min, longitude_max, latitude_min, latitude_max], cmap=colormap, norm=norme)
im = axes.imshow(O, extent=[longitude_min, longitude_max, latitude_min, latitude_max], cmap=colormap, norm=norme)
im = axes.imshow(P, extent=[longitude_min, longitude_max, latitude_min, latitude_max], cmap=colormap, norm=norme)

plt.show()

