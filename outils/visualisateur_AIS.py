input_file='../bdd.csv'

RT = 6311000 #rayon de la Terre (en m)
pi = 3.141592

time_scale=1000 #time_scale=x => le temps s'écoule x fois plus vite qu'en vrai

#changer la localisation de l'étude :
# latitude_min=55
# latitude_max=65
# longitude_min=-150
# longitude_max=-145
latitude_min = -90
latitude_max = 90
longitude_min = -180
longitude_max = 180
#DEBUG
frame_offset=0 #temps en secondes


interval=int(1000/time_scale) #ne pas toucher, modifier time_scale au-dessus

import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.widgets as widget
import cartopy.crs as ccrs #à installer avec "pip install Cartopy-0.18.0-cp38-cp38-win_amd64.whl" (problèmes de dépendances si installé depuis les dépots de Python)
import time

data=[] #contient les données du fichier d'entrée
boats={} #se remplit des bateaux détectés au fur et à mesure

#charger données

data=[]
with open(input_file, 'r', newline='') as input:
    csvreader = csv.reader(input, delimiter=',')
    for row in csvreader:

        mmsi=row[0]

        hours=int(row[1][11:13])
        minuts=int(row[1][14:16])
        seconds=int(row[1][17:19])
        temps=int(hours*3600+minuts*60+seconds)


        longitude=float(row[3])
        latitude=float(row[2])

        #speed over ground (SOG) : vitesse par rapport au sol (pas à la mer)
        SOG=float(row[4])*0.514 #1 noeud ~= 0.514 m/s
        #course over ground (COG) : direction de déplacement du bateau (pas son orientation), par rapport au Nord
        COG=float(row[5])*pi/180 #l'angle est donné en degrés

        vitesse_lat=SOG*np.cos(COG)/RT
        vitesse_lon=(SOG*np.sin(COG))/(RT*np.cos(latitude*pi/180)) #on considèrera que cos(lat) ne varie pas entre 2 mesures consécutives

        while temps > len(data)-1:
            data.append([]) #secondes pendant lesquelles rien est reçu
        data[temps].append([mmsi,[temps, longitude, latitude, vitesse_lon, vitesse_lat]])


len_data=len(data)

class Boat:
    def __init__(self, mmsi, vecteur): #vecteur = [instant, longitude, latitude, vitesse_longitudinale, vitesse_latitudinale]
        self.__mmsi = mmsi
        self.__liste_vecteurs = []
        self.__dot, = ax.plot([],[], marker='o',color='blue',markersize=5) #le point représentant le bateau sur la carte
        self.append(vecteur)
    def append(self, vecteur):
        self.__liste_vecteurs.append(vecteur)
        self.tracer()
    def tracer(self):
        self.__dot.set_data(self.__liste_vecteurs[-1][1], self.__liste_vecteurs[-1][2])


def update(frame):
    global data
    global boats
    temps=(frame+frame_offset) % len_data #permet de recommencer l'animation au début plutôt que de planter #frame_offset : DEBUG
    plt.title((str(temps//3600)+" h "+str((temps//60)%60)+" min "+str(temps%60)) + " ("+str(temps)+" secondes)")
    if not data[temps] == []: #si il y a des données "reçues" durant la seconde représentée par la frame
        for infos in data[temps]:
            mmsi=infos[0] #le bateau est identifié par son mmsi dans le programme
            if mmsi in boats: #Si le bateau est déjà enregistré,
                boats[mmsi].append(infos[1]) #on met a jour sa position sa vitesse et son angle,
            else: #sinon,
                boats[mmsi]=Boat(mmsi,infos[1]) #on en crée un nouveau.


fig = plt.figure() #il faut créer une figure pour l'animation (j'ai pas tout compris au fonctionnement de mathplot.pyplot. J'ai l'impression qu'il y a beaucoup d'objets qui font sensiblement la même chose
ax = plt.axes(projection=ccrs.PlateCarree(), autoscale_on=False, xlim=(longitude_min, longitude_max), ylim=(latitude_min, latitude_max))
ax.coastlines() #dessine la carte
ax.set_aspect('equal') #évite la déformation de la carte


ani = animation.FuncAnimation(fig, update, interval=interval) #cette fonction permet d'appeler la fonction "update" tous les interval ms
plt.show() #affiche la fenetre