"""Comparison with Birnbaum forest cover map in 2015."""

import os
from pathlib import Path

import zipfile
import requests
from osgeo import gdal, ogr
import geefcc
import ee

# Working directory
os.chdir("..")

# Defining region (grid with extent and resolution)
xmin = 349860
ymin = 7470900
xmax = 772200
ymax = 7841100
res = 30

# GDAL exceptions
ogr.UseExceptions()

# Reproject borders
data_raw = Path("data_raw")
gpkg_file = data_raw / "gadm-newcal" / "grande-terre.gpkg"
gpkg_file_proj = data_raw / "gadm-newcal" / "grande-terre_utm58s.gpkg"
gdal.VectorTranslate(
    gpkg_file_proj,   # fichier de sortie
    gpkg_file,          # fichier d'entrée
    dstSRS="EPSG:32758",    # UTM 58S
    format="GPKG"
)

# Rasterize borders of North and South Province
raster_file = data_raw / "gadm-newcal" / "grande-terre.tif"
land = gdal.Rasterize(
    raster_file,
    gpkg_file_proj,
    burnValues=[1],
    noData=0,
    outputBounds=[xmin, ymin, xmax, ymax],
    xRes=res,
    yRes=res,
    allTouched=False,
    targetAlignedPixels=True,
    outputType=gdal.GDT_Byte,
    creationOptions=["COMPRESS=DEFLATE"]
)
land = None

# ========================================================
# Get Birnbaum map (~2015)
# ========================================================

# Data download
url = ("https://zenodo.org/records/12731044/files/"
       "amap_carto_3k_20240715.zip?download=1")
destfile = data_raw / "amap_carto_3k_20240715.zip"

response = requests.get(url, stream=True)
with open(destfile, "wb") as f:
    for chunk in response.iter_content(chunk_size=8192):
        f.write(chunk)

# Uncompress
out_dir = os.path.join(data_raw, "amap_carto_3k_20240715")
os.makedirs(out_dir, exist_ok=True)
with zipfile.ZipFile(destfile, "r") as zip_ref:
    zip_ref.extractall(out_dir)

# Get shapefile
shp_file = None
for root, dirs, files in os.walk(out_dir):
    for file in files:
        if file.endswith(".shp"):
            shp_file = os.path.join(root, file)
            break

# Get CRS
ds = ogr.Open(shp_file)
layer = ds.GetLayer()
srs = layer.GetSpatialRef()
print(srs.ExportToWkt())
ds = None

# Shapefile -> GeoPackage
gpkg_file = os.path.join(out_dir, "amap_carto_3k_20240715.gpkg")
gdal.VectorTranslate(
    gpkg_file,
    shp_file,
    format="GPKG",
    layerName="forest_nc",
    dstSRS="EPSG:32758",
    reproject=True
)

# Rasterize
raster_file = os.path.join(out_dir, "forest_birnbaum.tif")
for2015 = gdal.Rasterize(
    raster_file,
    gpkg_file,
    layers=["forest_nc"],
    burnValues=[1],
    noData=0,
    outputBounds=[xmin, ymin, xmax, ymax],
    xRes=res,
    yRes=res,
    allTouched=True,  # To avoid underestimating
    targetAlignedPixels=True,
    outputType=gdal.GDT_Byte,
    creationOptions=["COMPRESS=DEFLATE"]
)
for2015 = None

# ========================================================
# Get GFC with geefcc
# ========================================================

ee.Initialize(project="deforisk",
              opt_url=("https://earthengine-highvolume."
                       "googleapis.com"))
ncpu = os.cpu_count() - 1

# GFC
out_dir = Path("data_raw", "out_geefcc", "gfc")
forest_file = out_dir / "for_2005_2015_2025_gfc.tif"
years = [2005, 2015, 2025]
geefcc.get_fcc_loss(
    aoi=os.path.join("data_raw", "gadm-newcal", "grande-terre.gpkg"),
    buff=0,
    years=years,
    source="gfc",
    perc=60,
    tile_size=1,
    output_file=forest_file,
    parallel=True,
    crop_to_aoi=True,
    ncpu=ncpu
)

# Sum raster bands
fcc_file = out_dir / "fcc_2005_2015_2025_gfc.tif"
geefcc.sum_raster_bands(
    input_file=forest_file,
    output_file=fcc_file,
    verbose=False,
)

# Compute statistics
res_df_gfc60 = geefcc.stat_fcc_loss(
    input_file=fcc_file,
    years=years,
    epsg=32758,
    output_file=out_dir / "fcc_stat_2005_2015_2025_gfc60.csv",
)

# ========================================================
# Get TMF with geefcc
# ========================================================

min_years = 1
out_dir = Path("data_raw", "out_geefcc", "tmf")
ofile = out_dir / "fcc_2005_2015_tmf.tif"
geefcc.get_fcc_loss_gain(
    aoi=os.path.join("data_raw", "gadm-newcal", "grande-terre.gpkg"),
    year1=2005,
    year2=2015,
    min_years=min_years,
    tile_size=1,
    crop_to_aoi=True,
    parallel=True,
    output_file=ofile
)


#### A FAIRE en utilisant rasterio et numpy
import rasterio
import numpy as np

# 1. Ouvrir les deux rasters d'entrée
with rasterio.open("chemin/raster1.tif") as src1, rasterio.open("chemin/raster2.tif") as src2:
    # Lire la bande 1 de chaque raster
    r1 = src1.read(1)
    r2 = src2.read(1)
    
    # Récupérer le NoData et le profil géographique du premier raster pour la sortie
    nodata1 = src1.nodata
    nodata2 = src2.nodata
    meta_out = src1.meta.copy()

# ---- OPTIONS DE CALCUL CONDITIONNEL ----

# Option A : Condition simple (si r1 > 10 ET r2 < 5 alors 1, sinon 0)
# Note : Utilisez & pour 'ET', | pour 'OU' (et entourez chaque bloc de parenthèses)
condition = (r1 > 10) & (r2 < 5)
raster_out = np.where(condition, 1, 0)


# Option B : Conditions multiples (Équivalent de if/else if/else)
# Exemple : si r1 == r2 -> 10 | si r1 > r2 -> 20 | sinon -> 30
conditions = [
    r1 == r2,
    r1 > r2
]
choix = [10, 20]
raster_out = np.select(conditions, choix, default=30)


# Option C : Gestion stricte des NoData (Style GRASS GIS)
# Si l'un des deux rasters est NoData, la sortie devient NoData
is_nodata = (r1 == nodata1) | (r2 == nodata2)
raster_out = np.where(is_nodata, nodata1, raster_out)

# ----------------------------------------

# 3. Écrire le raster final
meta_out.update(dtype=raster_out.dtype)
with rasterio.open("chemin/raster_resultat.tif", "w", **meta_out) as dst:
    dst.write(raster_out, 1)

