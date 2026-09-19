"""Comparison with AMAP forest cover map in 2015."""

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
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
# Get AMAP map (~2015)
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
amap_dir = out_dir / "amap"
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
raster_file = amap_dir / "forest_amap.tif"
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
raster2 = amap_dir / "forest_amap.tif"
result = amap_dir / "land_forest_amap.tif"
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
file = amap_dir / "land_forest_amap.tif"
with rasterio.open(file, "r") as src:
    for _, window in src.block_windows(1):
        data = src.read(1, window=window)
        land_area += np.sum(np.isin(data, [0, 1])) * pixel_area
        forest_area += np.sum(data == 1) * pixel_area
land_area_amap = np.round(land_area).astype(int)
forest_area_amap = np.round(forest_area).astype(int)
forest_perc_amap = np.round(100 * forest_area / land_area).astype(int)

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

# Confusion matrix Amap / GFC
rasterA = amap_dir / "land_forest_amap.tif"
rasterB = gfc_dir / "land_forest_gfc.tif"
otif = gfc_dir / "comp_amap_gfc.tif"
ocsv = comp = gfc_dir / "comp_amap_gfc.csv"
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

# Confusion matrix Amap / TMF
rasterA = amap_dir / "land_forest_amap.tif"
rasterB = tmf_dir / "land_forest_tmf.tif"
otif = tmf_dir / "comp_amap_tmf.tif"
ocsv = comp = tmf_dir / "comp_amap_tmf.csv"
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
    "source": ["amap", "gfc", "tmf"],
    "year": [2015, 2015, 2015],
    "tc": ["", 60, ""],
    "forest_area": [forest_area_amap, forest_area_gfc, forest_area_tmf],
    "perc": [forest_perc_amap, forest_perc_gfc, forest_perc_tmf],
    "OA": ["", OA_gfc, OA_tmf],
    "Kappa": ["", K_gfc, K_tmf]
}
ofile = Path("outputs", "comp_fc2021_amap_gfc_tmf.csv")
df = pd.DataFrame(data).to_csv(ofile, index=False)


# ========================================================
# Plots
# ========================================================

# TMF vs. Amap
# Colors
cols = [
    (211, 211, 211, 255),  # 0: light grey
    (34, 139, 34, 255),    # 1: green
    (10, 10, 150, 255),    # 2: blue
    (200, 200, 0, 255),    # 3: orange
    (255, 255, 255, 0)     # 255: transparent
]
colors = [tuple(i / 255.0 for i in col) for col in cols]
color_map = ListedColormap(colors)

# Map 255 for the last color
bounds_norm = [0, 1, 2, 3, 4, 256]
norm = BoundaryNorm(bounds_norm, color_map.N)

# Labels
labels = {
    0: "non-forest TMF / non-forest AMAP",
    1: "forest TMF / forest AMAP",
    2: "non-forest TMF / forest AMAP",
    3: "forest TMF / non-forest AMAP"
}
patches = [
    mpatches.Patch(facecolor=colors[i], edgecolor="black", label=labels[i])
    for i in range(4)
]

# File paths
ifile = tmf_dir / "comp_amap_tmf.tif"
opng = tmf_dir / "comp_amap_tmf.png"

# Generate overviews
factors = [2, 4, 8, 16]
with rasterio.open(ifile, "r+") as dst:
    dst.build_overviews(factors, resampling=Resampling.nearest)
# Select overview
OV_ID = 1  # factor 4

# Defining zoom
xmin, xmax = 525000, 575000
ymin, ymax = 7610000, 7660000

with rasterio.open(ifile, overview_level=OV_ID) as ds:
    data = ds.read(1)

    bounds = ds.bounds
    extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Main map
    im = ax.imshow(
        data,
        cmap=color_map,
        norm=norm,
        extent=extent,
        interpolation="nearest"
    )
    ax.set_aspect("equal")
    plt.xlim(bounds.left - 10000, bounds.right)
    plt.ylim(bounds.bottom, bounds.top + 10000)
    plt.title("Difference between TMF and AMAP 2015 forest cover maps")

    # Legend positioned above the inset
    leg = ax.legend(
        handles=patches,
        bbox_to_anchor=(1.05, 1.0),
        loc="upper left",
        borderaxespad=0.
    )

    # Rectangle for zoom
    WIDTH = xmax - xmin
    HEIGHT = ymax - ymin
    rect = mpatches.Rectangle(
        (xmin, ymin), WIDTH, HEIGHT,
        linewidth=2, edgecolor="black", facecolor="none"
    )
    ax.add_patch(rect)

    # Create the inset axis (positioned below the legend on the right)
    # bbox_to_anchor=(X, Y, width, height) relative to the main axis
    ax_inset = inset_axes(
        ax,
        width="100%",
        height="100%",
        # bbox_to_anchor=(1.0, 0.25, 0.45, 0.45),
        bbox_to_anchor=(0.84, 0.15, 0.65, 0.65),
        bbox_transform=ax.transAxes,
        loc="upper left"
    )

    # Read raw data for the inset
    with rasterio.open(ifile) as ds_full:
        window = rasterio.windows.from_bounds(
            xmin, ymin, xmax, ymax,
            transform=ds_full.transform
        )
        data_inset = ds_full.read(1, window=window)
        extent_inset = [xmin, xmax, ymin, ymax]

    # Display data in the inset and crop
    ax_inset.imshow(
        data_inset,
        cmap=color_map,
        norm=norm,
        extent=extent_inset,
        interpolation="nearest"
    )
    ax_inset.set_xlim(xmin, xmax)
    ax_inset.set_ylim(ymin, ymax)
    ax_inset.set_aspect("equal")

    # Optional: hide inset tick marks to lighten the visual
    ax_inset.set_xticks([])
    ax_inset.set_yticks([])

    # Save
    fig.savefig(opng, bbox_inches="tight", dpi=100)
    plt.close(fig)

# End
