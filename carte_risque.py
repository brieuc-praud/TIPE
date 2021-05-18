if __name__ == '__main__':#nécessaire pour le multiprocessing

    ##bibliothèques
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import cartopy.crs as ccrs#à installer avec "pip install Cartopy-0.18.0-cp38-cp38-win_amd64.whl" (problèmes de dépendances si installée depuis les dépots de Python)
    ##bibliothèques personelles
    import bib_ais as ais
    import bib_bateau as bat
    import bib_dijkstra as dij

    ##base de données à utiliser
    base_de_données='bdd.csv'

    ##localisation de l'étude
    n, echelle = 4, 200#la largeur de la carte est d'environ 500 km
    latitude_min = 42.5-1
    latitude_max = latitude_min + n
    longitude_min = -69.5-1
    longitude_max = longitude_min + n
    n = n*echelle

    ##autres paramètres
    colormap='jet'#couleurs utilisées pour représenter le risque
    norme = LogNorm()#None
    intervalle=1#temps d'attente entre deux tours de boucle (en millisecondes)


    ##initialisation
    data=[] #contient les données du fichier d'entrée
    boats={} #se remplit des bateaux détectés au fur et à mesure
    #counter=0#compteur pour numéroter les images sauvegardées automatiquement


    M=np.ones((n,n))
    mon_bateau = bat.MonBateau(n, echelle, longitude_min, latitude_min)

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
    compteur = 0
    ##pour quitter le programme lorsque l'on quitte la fenêtre
    trigger = False
    def handle_close(evt):
        global trigger
        trigger = True
    fig.canvas.mpl_connect('close_event', handle_close)

    ##fonction principale
    plt.ion()#mode interactif de pyplot : permet d'animer
    for temps in range(len_data):#on parcourt toutes les lignes de la base de données

        if trigger:#vérifie si la fenêtre à été fermée, auquel cas on termine la boucle
            break

        plt.title((str(temps//3600)+" h "+str((temps//60)%60)+" min "+str(temps%60)) + " ("+str(temps)+" secondes)")

        if not data[temps] == []: #si il y a des données "reçues" durant la seconde représentée par la frame
            for infos in data[temps]:
                mmsi=infos[0] #le bateau est identifié par son mmsi dans le programme
                if not (mmsi in boats): #Si le bateau n'est pas déjà enregistré on l'ajoute
                    boats[mmsi]=bat.Bateau(mmsi, infos[1], echelle, longitude_min, latitude_min)
                boats[mmsi].append(infos[1], M, mon_bateau)#puis on met a jour sa position, sa vitesse et son angle.

        # if temps%100 == 99:
        #      mon_bateau.calculer_plus_court_chemin(M)
        #      plt.savefig("DIJ"+str(compteur)+".svg")
        #      compteur += 1
        mon_bateau.calculer_plus_court_chemin(M)

        #affiche la carte des risques
        im.remove()
        im = axes.imshow(M, extent=[longitude_min, longitude_max, latitude_min, latitude_max], cmap=colormap, norm=norme, origin='lower')


        plt.pause(intervalle/1000)
