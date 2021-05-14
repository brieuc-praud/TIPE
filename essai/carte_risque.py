##bibliothèques
import time
from math import ceil

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs#à installer avec "pip install Cartopy-0.18.0-cp38-cp38-win_amd64.whl" (problèmes de dépendances si installée depuis les dépots de Python)
import multiprocessing as mp
##bibliothèques personelles
import bib_ais as ais
import bib_bateau as bat
import bib_dijkstra as dij

##base de données à utiliser
base_de_données='bdd.csv'

##localisation de l'étude
n, echelle = 5, 3e-5#la largeur de la carte est alors d'environ 300 km
latitude_min = 42.5-1
latitude_max = latitude_min + n
longitude_min = -69.5-1
longitude_max = longitude_min + n
#n = n*echelle

RT = 6311000 #rayon de la Terre (en m)
latitude_inter = 44#latitude intermédiaire
RTL = RT*np.cos(latitude_inter*np.pi/180)
n *= ceil(RT*echelle)


##autres paramètres
colormap='jet'#couleurs utilisées pour représenter le risque
norme = None#LogNorm()
intervalle=1#temps d'attente entre deux tours de boucle (en millisecondes)

vitesse_bateau_suivi=20#"vitesse" du bateau suivi
dep, arr = (0,100), (n-1,200)#depart/arrivee du bateau suivi


##initialisation
data=[] #contient les données du fichier d'entrée
boats={} #se remplit des bateaux détectés au fur et à mesure
#counter=0#compteur pour numéroter les images sauvegardées automatiquement

if __name__ == '__main__':
    manager = mp.Manager()
    return_dict = manager.dict()
    return_dict["LON"],return_dict["LAT"],return_dict["chemin"]=[],[],[]

Graphe = dij.graphe(n) #Le graphe pour Dijkstra

M=np.ones((n,n))
x=np.linspace(0,n, num=n)
y=np.linspace(0,n, num=n)
x, y = np.meshgrid(x, y)

#charger données
data = ais.charger_ais(base_de_données)
len_data=len(data)


fig = plt.figure()
axes = plt.axes(projection=ccrs.PlateCarree(), autoscale_on=False)
axes.coastlines() #dessine la carte
axes.set_aspect('equal') #évite la déformation de la carte
plt.xlim(longitude_min, longitude_max)
plt.ylim(latitude_min, latitude_max)

im = axes.imshow(M, cmap='seismic')#dessine la carte de risque

##processus qui calcule le chemin le plus court
def Dijkstra(M , ret_dict, dep, arr):#trouve et affiche le plus court chemin dans la matrice M
    N=M.copy()
    chemin = dij.trajet(Graphe,N,dep,arr)[1]
    LON, LAT = dij.preparer_trajet(chemin, echelle, latitude_min, longitude_min)
    ret_dict["LON"], ret_dict["LAT"], ret_dict["chemin"] = LON, LAT, chemin

dij_proc = mp.Process(target=None)

##pour quitter le programme lorsque l'on quitte la fenêtre
trigger = False
def handle_close(evt):
    global trigger
    trigger = True
fig.canvas.mpl_connect('close_event', handle_close)

##fonction principale
if __name__ == '__main__':# and not trig:#on vérifie que l'on est pas dans un thread
    plt.ion()#mode interactif de pyplot : permet d'animer
    for temps in range(len_data):#on parcourt toutes les lignes de la base de données

        if trigger:#vérifie si la fenêtre à été fermée, auquel cas on termine la boucle
            break

        plt.title((str(temps//3600)+" h "+str((temps//60)%60)+" min "+str(temps%60)) + " ("+str(temps)+" secondes)")

        if not data[temps] == []: #si il y a des données "reçues" durant la seconde représentée par la frame
            for infos in data[temps]:
                mmsi=infos[0] #le bateau est identifié par son mmsi dans le programme
                if True:#mmsi=="367362680":
                    #print(infos[1])
                    if mmsi in boats: #Si le bateau est déjà enregistré,
                            boats[mmsi].append(infos[1], M) #on met a jour sa position, sa vitesse et son angle,
                    else:#sinon,
                        boats[mmsi]=bat.Boat(mmsi, echelle, longitude_min, longitude_max, latitude_min, latitude_max) #on en crée un nouveau
                        boats[mmsi].append(infos[1], M)#puis on met a jour sa position, sa vitesse et son angle.

        #affiche la carte des risques
        im.remove()
        im = axes.imshow(M, extent=[longitude_min, longitude_max, latitude_min, latitude_max], cmap=colormap, norm=norme, origin='lower')

        # if not dij_proc.is_alive():#Si on n'est plus en train de calculer un nouvel itinéraire, alors on recommence.
        #     dij_proc = mp.Process(target=Dijkstra, args=(M, return_dict, dep, arr))
        #     dij_proc.start()
        #     chemin = return_dict["chemin"]
        #     if len(chemin):
        #         dep = (chemin[vitesse_bateau_suivi][0], chemin[vitesse_bateau_suivi][1])#on fait avancer le bateau suivi
        #     dij.affiche_trajet(return_dict["LON"], return_dict["LAT"])
            #counter+=1
            #plt.savefig('image'+str(counter)+'.pdf')

        #M=M*(1-1e-2) + 0.01

        plt.pause(intervalle/1000)
