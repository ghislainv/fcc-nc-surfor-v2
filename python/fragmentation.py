"""Forest fragmentation."""

from pathlib import Path

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

# End
