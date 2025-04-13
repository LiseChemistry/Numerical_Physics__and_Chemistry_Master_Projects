import numpy as np
import matplotlib.pyplot as plt

# Charger les données
fichier = "energie_data_10_L40.txt"
energie = []
cycles = []

# Chargement des données à partir du fichier
with open(fichier, 'r') as file:
    for line in file:
        parts = line.split(':')
        energie.append(float(parts[1].strip()))  # L'énergie
        cycles.append(float(parts[0].strip()))   # Le nombre de cycles

# Convertir les listes en arrays numpy
energie = np.array(energie)
cycles = np.array(cycles)

# Vérifier que le tableau contient au moins 1000 valeurs
if len(energie) >= 1000:
    # Extraire les 1000 dernières valeurs
    dernieres_valeurs = energie[-1000:]
    # Calculer la variance
    variance = np.var(dernieres_valeurs)
    print(f"Variance des 1000 dernières valeurs: {variance}")
else:
    print("Le tableau contient moins de 1000 valeurs.")

