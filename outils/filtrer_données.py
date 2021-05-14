#crée un fichier avec seulement les lignes répondant aux critères suivants

input_file='./zone19.csv'
output_file='../bdd.csv'

dates = ("2017-01-01")
heures=('00','01','02','03','04') #heures à sélectionner
#carré dans lequel les données vont être conservées
latitude_min=42.5
latitude_max=latitude_min+3
longitude_min=-69.5
longitude_max=longitude_min+3


import csv

with open(input_file, 'r', newline='') as input:
    with open(output_file, 'w', newline='') as output:
        input.readline() #la premiere ligne ne contient pas les données
        csvreader = csv.reader(input, delimiter=',')
        csvwriter = csv.writer(output, delimiter=',')
        for row in csvreader:
            mmsi=row[0]
            date=row[1][0:10]
            hour=row[1][11:13]
            lat=float(row[2])
            lon=float(row[3])
            if (date in dates) and (hour in heures) and lat>=latitude_min and lat<=latitude_max and lon>=longitude_min and lon<=longitude_max: #vérification de tous les critères
                csvwriter.writerow(row)

