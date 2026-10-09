"""Forest fragmentation."""

from pathlib import Path
import sys
import subprocess

import numpy as np
import rasterio
from osgeo import gdal
import matplotlib.pyplot as plt

# Project directory
proj_dir = Path.cwd().parent

# fcc file
fcc_file = proj_dir / ("outputs/fcc-map-2015-2025-newcal"
                       "/out_tmf_10yr/fcc_tmf_epsg32758.tif")

# Forest cover in 2025 (without considering afforestation)
out_dir = proj_dir / "outputs" / "fragmentation"
for2025 = out_dir / "for2025.tif"
with rasterio.open(fcc_file) as src:
    profile = src.profile
    profile.update(dtype="uint8", nodata=255, compress="deflate", predictor=2,
                   zlevel=6, tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(for2025, "w", **profile) as dst:
        for _, window in src.block_windows(1):
            # Data
            a = src.read(1, window=window)
            # Result
            res = np.ones_like(a) * 255
            res[np.isin(a, [1, 3, 4, 5])] = 1
            dst.write(res.astype("uint8"), 1, window=window)

# =======================
# Distance to forest edge
# =======================

input_file = for2025
dist_file = out_dir / "dist_edge.tif"
input_nodata = False
verbose = True
values = 255
nodata = 0
src_ds = gdal.Open(input_file)
srcband = src_ds.GetRasterBand(1)

# Create raster of distance
drv = gdal.GetDriverByName("GTiff")
dst_ds = drv.Create(
    dist_file,
    src_ds.RasterXSize,
    src_ds.RasterYSize,
    1,
    gdal.GDT_UInt32,
    ["COMPRESS=DEFLATE", "PREDICTOR=2", "BIGTIFF=YES"],
)
dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
dst_ds.SetProjection(src_ds.GetProjectionRef())
dstband = dst_ds.GetRasterBand(1)

# Use_input_nodata
ui_nodata = "YES" if input_nodata else "NO"

# Compute distance
val = "VALUES=" + str(values)
use_input_nodata = "USE_INPUT_NODATA=" + ui_nodata
cb = gdal.TermProgress_nocb if verbose else 0
gdal.ComputeProximity(
    srcband,
    dstband,
    [val, use_input_nodata, "DISTUNITS=GEO"],
    callback=cb
)

# Set nodata value
dstband.SetNoDataValue(nodata)

# Delete objects
srcband = None
dstband = None
del src_ds, dst_ds

# Histogram of distances
ofile = out_dir / "distance_edge_histogram.png"
bins = [1, 100, 250, 500, 1000, 5000]
hist = np.zeros(len(bins) - 1, dtype=int)
with rasterio.open(dist_file) as src:
    for _, window in src.block_windows(1):
        # Data
        a = src.read(1, window=window)
        # Result
        res, _ = np.histogram(a, bins=bins, density=False)
        hist += res
perc = np.round(100 * hist / hist.sum(), 1)

# Labels for x-axis
labels = [f"{bins[i]}-{bins[i+1]}" for i in range(len(bins)-1)]
labels[0] = "0-100"
labels[-1] = "> 1000"

fig, ax = plt.subplots(figsize=(10, 8))

bars = ax.bar(range(len(hist)), hist, tick_label=labels,
              color='steelblue', edgecolor='white')

# Values above bars
for bar, count, p in zip(bars, hist, perc):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height(),
        f"{p}%",
        ha='center', va='bottom', fontsize=16, fontweight='bold'
    )

ax.set_xlabel("Distance to forest edge (m)", fontsize=18, fontweight='bold')
ax.set_ylabel("Forest pixel count", fontsize=18, fontweight='bold')
ax.tick_params(axis='both', labelsize=14)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(ofile, dpi=300, bbox_inches='tight')
plt.close()

# =======================
# Forest patches
# =======================

# Install GRASS and l'API Python
# sudo apt install grass grass-dev

# Append GRASS to the python system path
sys.path.append(
    subprocess.check_output(["grass", "--config", "python_path"],
                            text=True).strip()
)
try:
    import grass.script as gs
    from grass.tools import Tools
    print("GRASS has been imported")
except ImportError as e:
    print(f"Erreur : {e}")

# Create a new project
project = "/tmp/grassproject_epsg_32758"
gs.create_project(project, epsg="32758")
print("GRASS GIS session:")
print(gs.parse_command("g.gisenv", flags="s"))

# Run GRASS tools
session = gs.setup.init(project)
tools = Tools(session=session)

# Import the raster
for2025 = out_dir / "for2025.tif"
tools.r_in_gdal(
    input=for2025,
    output="forest"
)
tools.g_region(raster="forest")

# Create a mask for the forest (value 1)
tools.r_mapcalc(
    expression="forest_mask = if(forest == 1, 1, null())",
    overwrite=True
)

# Find the patches (connected components)
# r.clump: groups adjacent cells into patches
tools.r_clump(
    input="forest_mask",
    output="patches",
    overwrite=True
)

# Calculate statistics per patch
# r.report (statistics per category)
rapport = tools.r_report(map="patches", units="meters", flags="h")
print(rapport)

# Export to an attribute table
tools.r_to_vect(input="patches", output="patches_vect", type="area")
tools.v_to_db(map="patches_vect", option="area", columns="area_ha",
              unit="hectares")
tools.v_to_db(map="patches_vect", option="compact", columns="compact")
tools.v_to_db(map="patches_vect", option="fd", columns="fd")
tools.v_to_db(map="patches_vect", option="perimeter", columns="perimeter",
              unit="kilometers")
ofile = out_dir / "forest_patches.csv"
tools.db_out_ogr(input="patches_vect", output=ofile, format="CSV")

# End
