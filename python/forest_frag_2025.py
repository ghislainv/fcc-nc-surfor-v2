"""Forest fragmentation in 2025."""

from pathlib import Path
import sys
import subprocess

import numpy as np
import rasterio
from osgeo import gdal
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# Project directory
proj_dir = Path.cwd().parent

# fcc file
fcc_file = proj_dir / ("outputs/7_fcc_map_2015_2025"
                       "/out_tmf_10yr/fcc_tmf_epsg32758.tif")

# Forest cover in 2025 (without considering afforestation)
out_dir = proj_dir / "outputs" / "8_forest_fragmentation_2025"
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
ofile = out_dir / "forest_frag_dist_edge_histogram.png"
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
# By default, the resulting clumps are connected only by their four direct neighbors
# (left, right, top, bottom). The -d flag activates also diagonal clump tracing.
for i in ["4n", "8n"]:
    diag = "d" if i == "8n" else None
    tools.r_clump(
        input="forest_mask",
        output=f"patches_{i}",
        flags=diag,
        overwrite=True,
    )

    # Calculate statistics per patch
    # r.report (statistics per category)
    rapport = tools.r_report(map=f"patches_{i}", units="meters", flags="h")

    # Export to an attribute table
    tools.r_to_vect(input=f"patches_{i}", output=f"patches_vect_{i}", type="area")
    tools.v_to_db(map=f"patches_vect_{i}", option="area", columns="area_ha",
                  unit="hectares")
    tools.v_to_db(map=f"patches_vect_{i}", option="compact", columns="compact")
    tools.v_to_db(map=f"patches_vect_{i}", option="fd", columns="fd")
    tools.v_to_db(map=f"patches_vect_{i}", option="perimeter", columns="perimeter",
                  unit="kilometers")
    ofile = out_dir / f"forest_frag_patches_{i}.csv"
    tools.db_out_ogr(input=f"patches_vect_{i}", output=ofile, format="CSV", overwrite=True)

    # Note: with 8n, the number of features (polygons) is the same,
    # but the values (categories) are changing.

# =======================
# GeoTIFF and plot
# =======================

min_ha = 0.5
tif_all = out_dir / "patches_4n.tif"
tif_min = out_dir / f"patches_4n_min{min_ha}ha.tif"

# Export the full patch raster (before filtering)
tools.r_out_gdal(
    input="patches_4n",
    output=tif_all,
    format="GTiff",
    type="Int32",
    createopt="COMPRESS=LZW",
    nodata=0,
    overwrite=True,
)

# Remove patches smaller than 0.5 ha
tools.r_reclass_area(
    input="patches_4n",
    output="patches_4n_min",
    value=min_ha,
    flags="c",  # input raster already clumped
    mode="greater",
    overwrite=True,
)
# r.reclass.area creates a virtual reclass raster: materialize it
tools.r_mapcalc(expression="patches_4n_clean = patches_4n_min", overwrite=True)

# Export the filtered raster (used for the plot)
tools.r_out_gdal(
    input="patches_4n_clean",
    output=tif_min,
    format="GTiff",
    type="Int32",
    createopt="COMPRESS=LZW",
    nodata=0,
    overwrite=True,
)

# Read the filtered raster
with rasterio.open(tif_min) as src:
    data = src.read(1)
    extent = (src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top)
    raster_crs = src.crs
    pixel_area_ha = abs(src.res[0] * src.res[1]) / 10_000  # assumes a projected CRS in meters

# Patch areas (hectares) from pixel counts; index = patch ID
counts = np.bincount(data.ravel())
counts[0] = 0  # ignore background / nodata
ids = np.nonzero(counts)[0]
areas_ha = counts[ids] * pixel_area_ha

# Select the 10 largest patches
order = np.argsort(areas_ha)[::-1]
top_ids = ids[order[:10]]
top_areas = areas_ha[order[:10]]

# Colors: pastel random colors for all patches (kept light so they do not compete)
rng = np.random.default_rng(42)
colors = 0.6 + 0.4 * rng.random((len(ids), 3))

# Strongly distinct colors for the 10 largest patches
top_colors = plt.get_cmap("tab10")(np.arange(10))[:, :3]
id_to_index = {pid: k for k, pid in enumerate(ids)}
for pid, col in zip(top_ids, top_colors):
    colors[id_to_index[pid]] = col

# Remap patch IDs to consecutive indices 1..N (0 reserved for background)
lut = np.zeros(ids.max() + 1, dtype=np.int32)
lut[ids] = np.arange(1, len(ids) + 1)
idx = lut[data]

# Colormap: white background + patch colors
cmap = ListedColormap(np.vstack([[1, 1, 1], colors]))

# New Caledonia borders, reprojected to the raster CRS
ifile = proj_dir / "gisdata/ncl_borders/borders_NCL.gpkg"
borders = gpd.read_file(ifile).to_crs(raster_crs)

# Subsample for display only (stats and top-10 selection use the full raster)
step = 4
idx_plot = idx[::step, ::step]

# Plot
fig, ax = plt.subplots(figsize=(10, 10))
ax.imshow(idx_plot, cmap=cmap, interpolation="nearest", extent=extent)
borders.boundary.plot(ax=ax, color="black", linewidth=0.8)

# Keep the view on the raster extent
ax.set_xlim(extent[0], extent[1])
ax.set_ylim(extent[2], extent[3])

# Legend for the 10 largest patches
top_pid = [i + 1 for i in range(10)]
handles = [
    Patch(facecolor=col, edgecolor="none", label=f"Patch {pid}: {area:,.0f} ha")
    for pid, col, area in zip(top_pid, top_colors, top_areas)
]
ax.legend(handles=handles, title="10 largest patches", loc="lower left",
          fontsize=8, frameon=True)

ax.set_title(f"Forest patches (4-neighbors, ≥ {min_ha} ha): {len(ids)} patches")
ax.set_axis_off()
plt.savefig(out_dir / "patches_4n.png", dpi=200, bbox_inches="tight")
plt.close()
    
# End
