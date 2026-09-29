"""
Extraction de données Google Earth Engine pour un fichier de points (lat/lon)
================================================================================

Ce script extrait, pour chaque point d'un fichier CSV :
  - Hansen Global Forest Change : treecover2000, loss, lossyear, gain, datamask
  - JRC Tropical Moist Forest (TMF) - Annual Changes : Dec1999, Dec2000,
    Dec2007, Dec2008, Dec2020, Dec2021

Les points sont traités par lots (chunks) pour éviter les timeouts / limites
de payload liées à getInfo() sur de grandes collections.

Pré-requis :
    pip install earthengine-api pandas

Auteur : généré avec Claude, relu et corrigé par Ghislain Vieilledent
"""

import time

import ee
import pandas as pd

# ---------------------------------------------------------------------------
# 1. CONFIGURATION - à adapter
# ---------------------------------------------------------------------------

PROJECT_ID = "deforisk"               # votre project ID GCP lié à Earth Engine
INPUT_CSV = "../data/df_surfor.csv"   # fichier d'entrée, colonnes 'lat' et 'lon'
OUTPUT_CSV = "resultats_gfc_tmf.csv"  # fichier de sortie
LAT_COL = "lat"
LON_COL = "lon"
CHUNK_SIZE = 500                          # points par lot (baisser si erreurs de timeout)
SCALE = 30                                # résolution en mètres (native Hansen/TMF)
MAX_RETRIES = 3                           # tentatives en cas d'erreur transitoire

# Bandes Hansen à extraire
HANSEN_ASSET = "UMD/hansen/global_forest_change_2025_v1_13"
HANSEN_BANDS = ["treecover2000", "loss", "lossyear", "gain", "datamask"]

# Bandes TMF Annual Changes à extraire
TMF_ASSET = "projects/JRC/TMF/v1_2025/AnnualChanges"
TMF_BANDS = ["Dec1999", "Dec2000", "Dec2007", "Dec2008", "Dec2020", "Dec2021"]


# ---------------------------------------------------------------------------
# 2. INITIALISATION EARTH ENGINE
# ---------------------------------------------------------------------------

def initialiser_ee(project_id):
    """Authentifie (si nécessaire) et initialise Earth Engine."""
    try:
        ee.Initialize(project=project_id)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=project_id)

        
# ---------------------------------------------------------------------------
# 3. CONSTRUCTION DE L'IMAGE COMBINÉE (Hansen + TMF)
# ---------------------------------------------------------------------------

def construire_image_combinee():
    hansen = ee.Image(HANSEN_ASSET).select(HANSEN_BANDS)
 
    # TMF Annual Changes est une ImageCollection tuilée (couverture mondiale
    # découpée en dalles) : il faut la mosaïquer pour obtenir une image unique
    # avant de l'assembler avec Hansen.
    tmf = ee.ImageCollection(TMF_ASSET).select(TMF_BANDS).mosaic()
 
    return hansen.addBands(tmf)


# ---------------------------------------------------------------------------
# 4. EXTRACTION PAR LOT
# ---------------------------------------------------------------------------

def extraire_chunk(df_chunk, image, lat_col, lon_col, scale):
    """Extrait les valeurs de `image` pour un sous-ensemble de points."""
    features = [
        ee.Feature(
            ee.Geometry.Point([row[lon_col], row[lat_col]]),
            # conserve l'index original du DataFrame pour la jointure
            {"id": idx, lat_col: row[lat_col], lon_col: row[lon_col]}
        )
        for idx, row in df_chunk.iterrows()
    ]
    fc = ee.FeatureCollection(features)

    resultats = image.reduceRegions(
        collection=fc,
        reducer=ee.Reducer.first(),
        scale=scale
    )

    info = resultats.getInfo()
    return [feat["properties"] for feat in info["features"]]


def extraire_tout(df, image, lat_col, lon_col, chunk_size, scale, max_retries):
    """Boucle sur tous les lots avec gestion des erreurs et retry."""
    tous_resultats = []
    n_chunks = (len(df) // chunk_size) + (1 if len(df) % chunk_size else 0)

    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i + chunk_size]
        num_lot = i // chunk_size + 1
        tentative = 0

        while tentative < max_retries:
            try:
                rows = extraire_chunk(chunk, image, lat_col, lon_col, scale)
                tous_resultats.extend(rows)
                print(f"Lot {num_lot}/{n_chunks} OK ({len(rows)} points)")
                break
            except Exception as e:
                tentative += 1
                print(f"Erreur lot {num_lot}, tentative {tentative}/{max_retries} : {e}")
                time.sleep(5)
        else:
            print(f"⚠️  Lot {num_lot} a échoué après {max_retries} tentatives — points ignorés")

    return tous_resultats


# ---------------------------------------------------------------------------
# 5. SCRIPT PRINCIPAL
# ---------------------------------------------------------------------------

def main():
    print("Initialisation Earth Engine...")
    initialiser_ee(PROJECT_ID)

    print(f"Lecture du fichier d'entrée : {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    df = df.reset_index(drop=True)
    print(f"{len(df)} points chargés.")

    print("Construction de l'image combinée (Hansen + TMF)...")
    image_combinee = construire_image_combinee()

    print("Début de l'extraction par lots...")
    resultats = extraire_tout(
        df,
        image_combinee,
        LAT_COL,
        LON_COL,
        CHUNK_SIZE,
        SCALE,
        MAX_RETRIES
    )

    result_df = pd.DataFrame(resultats)
    print(f"\nExtraction terminée : {result_df.shape[0]} lignes récupérées "
          f"sur {len(df)} points d'origine.")

    # Jointure avec le DataFrame original (lat/lon + attributs éventuels)
    if "id" in result_df.columns:
        result_df = result_df.set_index("id")
        final_df = df.join(result_df, how="left")
    else:
        final_df = result_df

    final_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Résultats sauvegardés dans : {OUTPUT_CSV}")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------
# ALTERNATIVE : export via Google Drive (plus robuste pour de très gros
# volumes, ne dépend pas de la stabilité de la connexion locale).
# ---------------------------------------------------------------------

initialiser_ee(PROJECT_ID)
df = pd.read_csv(INPUT_CSV).reset_index(drop=True)
image_combinee = construire_image_combinee()

features = [
    ee.Feature(
        ee.Geometry.Point([row[LON_COL], row[LAT_COL]]),
        {"id": idx, LAT_COL: row[LAT_COL], LON_COL: row[LON_COL]}
    )
    for idx, row in df.iterrows()
]
fc = ee.FeatureCollection(features)

resultats = image_combinee.reduceRegions(fc, ee.Reducer.first(), scale=SCALE)

task = ee.batch.Export.table.toDrive(
    collection=resultats,
    description="extracting_gfc_tmf",
    fileFormat="CSV"
)
task.start()

while task.active():
    print("En cours...", task.status()["state"])
    time.sleep(30)
print("Terminé :", task.status())

# Le fichier CSV apparaît ensuite dans votre Google Drive.
