"""Comparison with Birnbaum forest cover map in 2015."""

import os
from pathlib import Path

import zipfile
import requests
from osgeo import gdal, ogr
import geefcc
import ee
import rasterio
from rasterio.enums import Resampling
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

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
out_dir = Path("outputs", "out_forest2015")
gpkg_file = out_dir / "gadm-newcal" / "grande-terre.gpkg"
gpkg_file_proj = out_dir / "gadm-newcal" / "grande-terre_utm58s.gpkg"
gdal.VectorTranslate(
    str(gpkg_file_proj),   # fichier de sortie
    str(gpkg_file),          # fichier d'entrée
    dstSRS="EPSG:32758",    # UTM 58S
    format="GPKG"
)

# Rasterize borders of North and South Province
raster_file = out_dir / "gadm-newcal" / "grande-terre.tif"
land = gdal.Rasterize(
    str(raster_file),
    str(gpkg_file_proj),
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
destfile = out_dir / "amap_carto_3k_20240715.zip"

response = requests.get(url, stream=True)
with open(destfile, "wb") as f:
    for chunk in response.iter_content(chunk_size=8192):
        f.write(chunk)

# Uncompress
amap_dir = out_dir / "amap_carto_3k_20240715"
os.makedirs(amap_dir, exist_ok=True)
with zipfile.ZipFile(destfile, "r") as zip_ref:
    zip_ref.extractall(amap_dir)

# Get shapefile
shp_file = None
for root, dirs, files in os.walk(amap_dir):
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
gpkg_file = amap_dir / "amap_carto_3k_20240715.gpkg"
gdal.VectorTranslate(
    str(gpkg_file),
    str(shp_file),
    format="GPKG",
    layerName="forest_nc",
    dstSRS="EPSG:32758",
    reproject=True
)

# Rasterize
raster_file = amap_dir / "forest_birn.tif"
for2015 = gdal.Rasterize(
    str(raster_file),
    str(gpkg_file),
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

# Combining forest with land
raster1 = out_dir / "gadm-newcal" / "grande-terre.tif"
raster2 = amap_dir / "forest_birn.tif"
result = amap_dir / "land_forest_birn.tif"
with rasterio.open(raster1) as src1, rasterio.open(raster2) as src2:
    profile = src1.profile
    profile.update(dtype="uint8", nodata=255, compress="deflate", predictor=2,
                   zlevel=6, tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(result, "w", **profile) as dst:
        for _, window in src1.block_windows(1):
            # Data
            a = src1.read(1, window=window)
            b = src2.read(1, window=window)
            # Result
            res = np.ones_like(a) * 255
            res[np.where((a == 1) & (b == 0))] = 0
            res[np.where((a == 1) & (b == 1))] = 1
            dst.write(res.astype("uint8"), 1, window=window)

# Areas
land_area = 0
forest_area = 0
pixel_area = 30 * 30 / 10000
file = amap_dir / "land_forest_birn.tif"
with rasterio.open(file, "r") as src:
    for _, window in src.block_windows(1):
        data = src.read(1, window=window)
        land_area += np.sum(np.isin(data, [0, 1])) * pixel_area
        forest_area += np.sum(data == 1) * pixel_area
land_area_birn = np.round(land_area).astype(int)
forest_area_birn = np.round(forest_area).astype(int)
forest_perc_birn = np.round(100 * forest_area / land_area).astype(int)

# ========================================================
# Get GFC with geefcc
# ========================================================

ee.Initialize(project="deforisk",
              opt_url=("https://earthengine-highvolume."
                       "googleapis.com"))
ncpu = os.cpu_count() - 1

# GFC
gfc_dir = out_dir / "gfc"
forest_file = gfc_dir / "for_2005_2015_2025_gfc.tif"
years = [2005, 2015, 2025]
geefcc.get_fcc_loss(
    aoi=out_dir / "gadm-newcal" / "grande-terre.gpkg",
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
fcc_file = gfc_dir / "fcc_2005_2015_2025_gfc.tif"
geefcc.sum_raster_bands(
    input_file=forest_file,
    output_file=fcc_file,
    verbose=False,
)

# Compute statistics (and reproject)
ofile = gfc_dir / "fcc_stat_2005_2015_2025_gfc60.csv"
res_df_gfc60 = geefcc.stat_fcc_loss(
    input_file=fcc_file,
    years=years,
    epsg=32758,
    output_file=ofile,
)

# Reproject on region's grid
input_file = gfc_dir / "fcc_2005_2015_2025_gfc.tif"
proj_file = gfc_dir / "fcc_2005_2015_2025_gfc_proj.tif"
ds = gdal.Warp(
    str(proj_file), str(input_file),
    outputBounds=[xmin, ymin, xmax, ymax],
    xRes=30, yRes=30,
    dstSRS="EPSG:32758",
    resampleAlg="near",
    targetAlignedPixels=True,
    creationOptions=["COMPRESS=DEFLATE"],
)
ds = None

# Combining forest with land
raster1 = out_dir / "gadm-newcal" / "grande-terre.tif"
raster2 = gfc_dir / "fcc_2005_2015_2025_gfc_proj.tif"
result = gfc_dir / "land_forest_gfc.tif"
with rasterio.open(raster1) as src1, rasterio.open(raster2) as src2:
    profile = src1.profile
    profile.update(dtype="uint8", nodata=255, compress="deflate", predictor=2,
                   zlevel=6, tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(result, "w", **profile) as dst:
        for _, window in src1.block_windows(1):
            # Data
            a = src1.read(1, window=window)
            b = src2.read(1, window=window)
            # Result
            res = np.ones_like(a) * 255
            res[np.where((a == 1) & (b <= 1))] = 0
            res[np.where((a == 1) & ((b == 2) | (b == 3)))] = 1
            dst.write(res.astype("uint8"), 1, window=window)

# Areas
land_area = 0
forest_area = 0
pixel_area = 30 * 30 / 10000
file = gfc_dir / "land_forest_gfc.tif"
with rasterio.open(file, "r") as src:
    for _, window in src.block_windows(1):
        data = src.read(1, window=window)
        land_area += np.sum(np.isin(data, [0, 1])) * pixel_area
        forest_area += np.sum(data == 1) * pixel_area
land_area_gfc = np.round(land_area).astype(int)
forest_area_gfc = np.round(forest_area).astype(int)
forest_perc_gfc = np.round(100 * forest_area / land_area).astype(int)

# ========================================================
# Get TMF with geefcc
# ========================================================

min_years = 1
tmf_dir = out_dir / "tmf"
ofile = tmf_dir / "fcc_2005_2015_tmf.tif"
geefcc.get_fcc_loss_gain(
    aoi=out_dir / "gadm-newcal" / "grande-terre.gpkg",
    year1=2005,
    year2=2015,
    min_years=min_years,
    tile_size=1,
    crop_to_aoi=True,
    parallel=True,
    output_file=ofile
)

# Reproject on region's grid
input_file = tmf_dir / "fcc_2005_2015_tmf.tif"
proj_file = tmf_dir / "fcc_2005_2015_tmf_proj.tif"
ds = gdal.Warp(
    str(proj_file), str(input_file),
    outputBounds=[xmin, ymin, xmax, ymax],
    xRes=30, yRes=30,
    dstSRS="EPSG:32758",
    resampleAlg="near",
    targetAlignedPixels=True,
    creationOptions=["COMPRESS=DEFLATE"],
)
ds = None

# Combining forest with land
raster1 = out_dir / "gadm-newcal" / "grande-terre.tif"
raster2 = tmf_dir / "fcc_2005_2015_tmf_proj.tif"
result = tmf_dir / "land_forest_tmf.tif"
with rasterio.open(raster1) as src1, rasterio.open(raster2) as src2:
    profile = src1.profile
    profile.update(dtype="uint8", nodata=255, compress="deflate", predictor=2,
                   zlevel=6, tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(result, "w", **profile) as dst:
        for _, window in src1.block_windows(1):
            # Data
            a = src1.read(1, window=window)
            b = src2.read(1, window=window)
            # Result
            res = np.ones_like(a) * 255
            res[np.where((a == 1) & (b != 1))] = 0
            res[np.where((a == 1) & (b == 1))] = 1
            dst.write(res.astype("uint8"), 1, window=window)

# Areas
land_area = 0
forest_area = 0
pixel_area = 30 * 30 / 10000
file = tmf_dir / "land_forest_tmf.tif"
with rasterio.open(file, "r") as src:
    for _, window in src.block_windows(1):
        data = src.read(1, window=window)
        land_area += np.sum(np.isin(data, [0, 1])) * pixel_area
        forest_area += np.sum(data == 1) * pixel_area
land_area_tmf = np.round(land_area).astype(int)
forest_area_tmf = np.round(forest_area).astype(int)
forest_perc_tmf = np.round(100 * forest_area / land_area).astype(int)

# ========================================================
# Comparison
# ========================================================

# Confusion matrix Birn / GFC
rasterA = amap_dir / "land_forest_birn.tif"
rasterB = gfc_dir / "land_forest_gfc.tif"
otif = gfc_dir / "comp_birn_gfc.tif"
ocsv = comp = gfc_dir / "comp_birn_gfc.csv"
(ff, fn, nf, nn) = (0, 0, 0, 0)
with rasterio.open(rasterA) as srcA, rasterio.open(rasterB) as srcB:
    profile = srcA.profile
    profile.update(dtype="uint8", nodata=255, compress="deflate", predictor=2,
                   zlevel=6, tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(otif, "w", **profile) as dst:
        for _, window in srcA.block_windows(1):
            # Data
            A = srcA.read(1, window=window)
            B = srcB.read(1, window=window)
            # Map
            res = np.ones_like(A) * 255
            res[np.where((A == 0) & (B == 0))] = 0
            res[np.where((A == 1) & (B == 1))] = 1
            res[np.where((A == 1) & (B == 0))] = 2
            res[np.where((A == 0) & (B == 1))] = 3
            dst.write(res.astype("uint8"), 1, window=window)
            # Confusion matrix
            ff += np.sum((A == 1) & (B == 1))
            fn += np.sum((A == 1) & (B == 0))
            nf += np.sum((A == 0) & (B == 1))
            nn += np.sum((A == 0) & (B == 0))
arr = np.round(np.array([ff, fn, nf, nn]) * pixel_area).astype(int)
arr = arr.reshape(2, 2, order="F")
df = pd.DataFrame(arr)
df = df.rename(columns={0: "f_ref", 1: "n_ref"},
               index={0: "f_map", 1: "n_map"}).to_csv(ocsv)
# Accuracy
total_pixels = arr.sum()
diagonal_sum = np.diag(arr).sum()
overall_accuracy = diagonal_sum / total_pixels
# Kappa
row_totals = arr.sum(axis=1)
col_totals = arr.sum(axis=0)
pe = np.sum(row_totals * col_totals) / (total_pixels ** 2)
kappa = (overall_accuracy - pe) / (1 - pe)
(OA_gfc, K_gfc) = np.round((overall_accuracy, kappa), 2)

# Confusion matrix Birn / TMF
rasterA = amap_dir / "land_forest_birn.tif"
rasterB = tmf_dir / "land_forest_tmf.tif"
otif = tmf_dir / "comp_birn_tmf.tif"
ocsv = comp = tmf_dir / "comp_birn_tmf.csv"
(ff, fn, nf, nn) = (0, 0, 0, 0)
with rasterio.open(rasterA) as srcA, rasterio.open(rasterB) as srcB:
    profile = srcA.profile
    profile.update(dtype="uint8", nodata=255, compress="deflate", predictor=2,
                   zlevel=6, tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(otif, "w", **profile) as dst:
        for _, window in srcA.block_windows(1):
            # Data
            A = srcA.read(1, window=window)
            B = srcB.read(1, window=window)
            # Map
            res = np.ones_like(A) * 255
            res[np.where((A == 0) & (B == 0))] = 0
            res[np.where((A == 1) & (B == 1))] = 1
            res[np.where((A == 1) & (B == 0))] = 2
            res[np.where((A == 0) & (B == 1))] = 3
            dst.write(res.astype("uint8"), 1, window=window)
            # Confusion matrix
            ff += np.sum((A == 1) & (B == 1))
            fn += np.sum((A == 1) & (B == 0))
            nf += np.sum((A == 0) & (B == 1))
            nn += np.sum((A == 0) & (B == 0))
arr = np.round(np.array([ff, fn, nf, nn]) * pixel_area).astype(int)
arr = arr.reshape(2, 2, order="F")
df = pd.DataFrame(arr)
df = df.rename(columns={0: "f_ref", 1: "n_ref"},
               index={0: "f_map", 1: "n_map"}).to_csv(ocsv)
# Accuracy
total_pixels = arr.sum()
diagonal_sum = np.diag(arr).sum()
overall_accuracy = diagonal_sum / total_pixels
# Kappa
row_totals = arr.sum(axis=1)
col_totals = arr.sum(axis=0)
pe = np.sum(row_totals * col_totals) / (total_pixels ** 2)
kappa = (overall_accuracy - pe) / (1 - pe)
(OA_tmf, K_tmf) = np.round((overall_accuracy, kappa), 2)

# Table synthesizing the results
data = {
    "source": ["birn", "gfc", "tmf"],
    "year": [2015, 2015, 2015],
    "tc": ["", 60, ""],
    "forest_area": [forest_area_birn, forest_area_gfc, forest_area_tmf],
    "perc": [forest_perc_birn, forest_perc_gfc, forest_perc_tmf],
    "OA": ["", OA_gfc, OA_tmf],
    "Kappa": ["", K_gfc, K_tmf]
}
ofile = Path("outputs", "comp_fc2021_birn_gfc_tmf.csv")
df = pd.DataFrame(data).to_csv(ofile, index=False)

# ========================================================
# Plots
# ========================================================

# TMF vs. Birn
# Colors
cols = [(0, 0, 0, 255),      # black for 0
        (34, 139, 34, 255),  # green for 1
        (10, 10, 150, 255),  # blue for 2
        (200, 200, 0, 255),  # orange for 3
        (255, 255, 255, 0)]  # transparent white for 255
cmax = 255.0
colors = []
for col in cols:
    col_class = tuple([i / cmax for i in col])
    colors.append(col_class)
color_map = ListedColormap(colors)

# Labels
labels = {0: "non-forest tmf / non-forest amap", 1: "forest tmf / forest amap",
          2: "non-forest tmf / forest amap", 3: "forest tmf / non-forest amap"}
patches = [mpatches.Patch(
    facecolor=col, edgecolor="black",
    label=labels[i]) for (i, col) in enumerate(colors[:-1])]

# 1. File paths
ifile = tmf_dir / "comp_birn_tmf.tif"
opng = tmf_dir / "comp_birn_tmf.png"

# 2. Generate Overviews
factors = [2, 4, 8, 16]
with rasterio.open(ifile, "r+") as dst:
    dst.build_overviews(
        factors,
        resampling=Resampling.nearest
    )

# 3. Choose the overview zoom level for plotting
# index 0 = factor 2, 1 = factor 4, 2 = factor 8, 3 = factor 16
overview_idx = 2

# 4. Read data and Plot
with rasterio.open(ifile) as ds:
    # Get the downsampling factor from the list (e.g., factors[2] = 8)
    factor = ds.overviews(1)[overview_idx]

    # Calculate target dimensions by dividing original sizes by the factor
    ov_height = int(ds.height / factor)
    ov_width = int(ds.width / factor)

    # Read only the downsampled pixels from the selected overview
    data_resampled = ds.read(
        1, 
        out_shape=(ov_height, ov_width),
        resampling=Resampling.nearest
    )

    # Get geographic boundaries (the spatial extent remains identical)
    bounds = ds.bounds
    extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]

    # Initialize the plot layout
    fig = plt.figure()
    ax = plt.subplot(111)

    # Plot the lightweight data array (geographic extent maps correctly)
    ax.imshow(data_resampled, cmap=color_map, extent=extent, interpolation="nearest")

    ax.set_aspect("equal")
    plt.title("Difference between TMF and AMAP 2015 forest cover maps")
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    # Save the output image and free up memory
    fig.savefig(opng, bbox_inches="tight", dpi=100)
    plt.close(fig)

# End
