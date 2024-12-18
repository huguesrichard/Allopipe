import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Charger le fichier TSV dans un DataFrame
data = pd.read_csv('../data/distrib_AMS.tsv', sep='\t')

# Filtrer les données pour les "related"
data_related = data[data['Relatedness'] == 'related']

# Filtrer les données pour les "unrelated"
data_non_related = data[data['Relatedness'] == 'unrelated']

# Extraire la colonne d'intérêt (AMS) pour chaque groupe
AMS_related = data_related['AMS']  # Extraire la colonne 'AMS' du groupe "related"
AMS_unrelated = data_non_related['AMS']  # Extraire la colonne 'AMS' du groupe "unrelated"

# Tracer les distributions des deux groupes
plt.figure(figsize=(10, 6))

# Tracer la distribution pour le groupe "related"
sns.histplot(AMS_related, kde=True, color='blue', label='Related pairs', stat='density', bins=30)

# Tracer la distribution pour le groupe "unrelated"
sns.histplot(AMS_unrelated, kde=True, color='red', label='Unrelated pairs', stat='density', bins=30)

# Définir un répertoire pour la recherche des fichiers
repertoire = './AMS_run'

# Liste tous les fichiers dans le répertoire, en excluant les sous-répertoires
fichiers = [f for f in os.listdir(repertoire) if not os.path.isdir(os.path.join(repertoire, f))]

# Trouver un fichier .tsv ou .csv
fichier_tsv = next((f for f in fichiers if f.endswith('.tsv')), None)
fichier_csv = next((f for f in fichiers if f.endswith('.csv')), None)

# Vérifier la priorité du fichier .tsv si disponible, sinon utiliser .csv
if fichier_tsv:
    # Charger le fichier .tsv
    df_ams = pd.read_csv(os.path.join(repertoire, fichier_tsv), sep='\t')
    print(f"Le fichier {fichier_tsv} a été chargé.")

    # Vérifier que la colonne 'ams_giab' existe dans le fichier .tsv
    if 'ams_giab' in df_ams.columns:
        valeurs_speciales = df_ams['ams_giab'].values  # Extraire les valeurs de la colonne 'ams_giab'

        # Ajouter une ligne verticale pour chaque valeur spéciale et positionner le label
        for ams in valeurs_speciales:
            # Tracer la ligne verticale
            plt.axvline(ams, color='red', linestyle='--')
    else:
        print("La colonne 'ams_giab' n'a pas été trouvée dans le fichier .tsv.")

elif fichier_csv:
    # Charger le fichier .csv
    df_ams = pd.read_csv(os.path.join(repertoire, fichier_csv))
    print(f"Le fichier {fichier_csv} a été chargé.")

    # Vérifier que la colonne 'ams' existe dans le fichier .csv
    if 'ams' in df_ams.columns:
        valeurs_speciales = df_ams['ams'].values  # Extraire les valeurs de la colonne 'ams'

        # Ajouter une ligne verticale pour chaque valeur spéciale et positionner le label
        for ams in valeurs_speciales:
            # Tracer la ligne verticale
            plt.axvline(ams, color='red', linestyle='--', label=f'AMS of the pair : {ams}')
    else:
        print("La colonne 'ams' n'a pas été trouvée dans le fichier .csv.")

else:
    print("Aucun fichier .tsv ou .csv trouvé dans le répertoire.")

# Ajouter des labels et une légende
plt.title('AMS Distribution Related vs Unrelated + computed AMS')
plt.xlabel('AMS values')
plt.ylabel('Density')
plt.legend()

# Placer la légende en haut du graphique, centrée
plt.legend(loc='upper left')


# Enregistrer le graphique sous format PNG
plt.savefig("distrib-fambiomnis-shufflecryostem.png", format='png')

# Afficher le graphique
plt.show()
