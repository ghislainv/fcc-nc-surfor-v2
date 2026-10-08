"""Get fcc map with geefcc."""

from pathlib import Path
import time

import ee
import geefcc
import numpy as np
import pandas as pd

# Initialize GEE
ee.Initialize(
    project="deforisk",
    opt_url="https://earthengine-highvolume.googleapis.com")

# Working directory
wd = Path.cwd().parent / "outputs/fcc-map-2015-2025-newcal/"

# Download data from GEE
min_years = 10
out_dir = wd / f"out_tmf_{min_years}yr"
ofile = out_dir / "fcc_tmf.tif"
start_time = end_time = time.time()
if not ofile.is_file():
    start_time = time.time()
    geefcc.get_fcc_loss_gain(
        aoi=(163.5, -23, 168.15, -19.51),
        year1=2015,
        year2=2025,
        min_years=min_years,
        tile_size=1,
        crop_to_aoi=True,
        parallel=True,
        output_file=ofile,
    )
    end_time = time.time()

elapsed_time = (end_time - start_time) / 60
print('Execution time:', round(elapsed_time, 2), 'minutes')

# Plot
vfile = wd / "data" / "borders_NCL.gpkg"
geefcc.plot_fcc_loss_gain(
    input_file=ofile,
    output_file=out_dir / f"fcc_tmf_{min_years}yr.png",
    title="Forest cover change 2015\u20132025, TMF",
    dpi=200,
    borders=vfile,
    grid=out_dir / "grid.gpkg",
    xlim=(163, 169),
    ylim=(-23.25, -18.75),
)

# Stat
ifile = out_dir / "fcc_tmf.tif"
ofile = out_dir / f"fcc_statistics_{min_years}yr.csv"
res_df = geefcc.stat_fcc_loss_gain(
    input_file=ifile,
    epsg=32758,
    output_file=ofile)

# Summary
forest_t1 = res_df.loc[[1, 2, 4, 5, 6, 7, 8], "area_ha"].sum()
forest_t2 = res_df.loc[[1, 3, 4, 5, 7, 9], "area_ha"].sum()
gross_loss = - res_df.loc[[2, 4, 6, 8], "area_ha"].sum()
gross_gain = res_df.loc[[3, 4, 9], "area_ha"].sum()
lossgain_df = pd.DataFrame({
    "label": ["forest_t1", "forest_t2", "gross loss", "gross gain", "net change"],
    "area_ha": [forest_t1, forest_t2, gross_loss, gross_gain, gross_gain + gross_loss],
})
time = 2025 - 2015
lossgain_df["annual_change_ha"] = (lossgain_df["area_ha"] / time).round().astype(int)
ratio = lossgain_df["area_ha"] / forest_t1
lossgain_df["annual_change_perc"] = round(100 * (1 - pow((1 - ratio), 1 / time)), 2)
lossgain_df.iloc[:2, 2:4] = np.nan
# Export
ofile = out_dir / f"loss_gain_statistics_{min_years}yr.csv"
lossgain_df.to_csv(ofile, index=False)

# Summary removing afforestation
forest_t1 = res_df.loc[[1, 2, 4, 5, 6, 7], "area_ha"].sum()
forest_t2 = res_df.loc[[1, 3, 4, 5, 7], "area_ha"].sum()
gross_loss = - res_df.loc[[2, 4, 6], "area_ha"].sum()
gross_gain = res_df.loc[[3, 4], "area_ha"].sum()
lossgain_df = pd.DataFrame({
    "label": ["forest_t1", "forest_t2", "gross loss", "gross gain", "net change"],
    "area_ha": [forest_t1, forest_t2, gross_loss, gross_gain, gross_gain + gross_loss],
})
time = 2025 - 2015
lossgain_df["annual_change_ha"] = (lossgain_df["area_ha"] / time).round().astype(int)
ratio = lossgain_df["area_ha"] / forest_t1
lossgain_df["annual_change_perc"] = round(100 * (1 - pow((1 - ratio), 1 / time)), 2)
lossgain_df.iloc[:2, 2:4] = np.nan
# Export
ofile = out_dir / f"loss_gain_statistics_{min_years}yr_noaff.csv"
lossgain_df.to_csv(ofile, index=False)

# End
