import csv
import numpy as np
import bib_bateau as bateau

RT = 6311000 #rayon de la Terre (en m)
pi = 3.141592

def charger_ais(input_file):#pour préparer la simulation en temps réel
    data=[]
    with open(input_file, 'r', newline='') as input:
        csvreader = csv.reader(input, delimiter=',')
        for row in csvreader:

            mmsi=row[0]

            hours=int(row[1][11:13])
            minuts=int(row[1][14:16])
            seconds=int(row[1][17:19])
            temps=int(hours*3600+minuts*60+seconds)


            latitude=float(row[2])
            longitude=float(row[3])

            #speed over ground (SOG) : vitesse par rapport au sol (pas à la mer)
            SOG=float(row[4])*0.514 #1 noeud ~= 0.514 m/s
            #course over ground (COG) : direction de déplacement du bateau (pas son orientation), par rapport au Nord
            COG=float(row[5])*pi/180 #l'angle est donné en degrés

            vitesse_lat=SOG*np.cos(COG)/RT
            vitesse_lon=(SOG*np.sin(COG))/(RT*np.cos(latitude*pi/180)) #on considèrera que cos(lat) ne varie pas entre 2 mesures consécutives

            while temps > len(data)-1:
                data.append([]) #secondes pendant lesquelles rien est reçu
            data[temps].append([mmsi,[temps, longitude, latitude, vitesse_lat, vitesse_lon]])

    return data