"""
Google Earth Engine Data Extraction for a Points File (lat/lon)
================================================================

This script extracts, for each point in a CSV file:
  - Hansen Global Forest Change: treecover2000, loss, lossyear, gain, datamask
  - JRC Tropical Moist Forest (TMF) - Annual Changes: Dec1999, Dec2000,
    Dec2007, Dec2008, Dec2020, Dec2021

Points are processed in chunks to avoid timeouts / payload limits associated
with getInfo() on large collections.

Prerequisites:
    pip install earthengine-api pandas

Author: generated with Claude, proofread and corrected by Ghislain Vieilledent
"""

import time

import ee
import pandas as pd

# ---------------------------------------------------------------------------
# 1. CONFIGURATION - adapt as needed
# ---------------------------------------------------------------------------

PROJECT_ID = "deforisk"                     # your GCP project ID linked to Earth Engine
INPUT_CSV = "../data/df_surfor_for_gee.csv"  # input file, must contain 'lat' and 'lon' columns
OUTPUT_CSV = "extract_gee_gfc_tmf.csv"   # output file
LAT_COL = "lat"
LON_COL = "lon"
CHUNK_SIZE = 500                            # points per chunk (reduce if you encounter timeouts)
SCALE = 30                                  # resolution in meters (native Hansen/TMF)
MAX_RETRIES = 3                             # retries for transient errors

# Hansen bands to extract
HANSEN_ASSET = "UMD/hansen/global_forest_change_2025_v1_13"
HANSEN_BANDS = ["treecover2000", "loss", "lossyear", "gain", "datamask"]

# TMF Annual Changes bands to extract
TMF_ASSET = "projects/JRC/TMF/v1_2025/AnnualChanges"
TMF_BANDS = ["Dec2000", "Dec2003", "Dec2004", "Dec2007",
             "Dec2008", "Dec2020", "Dec2021"]


# ---------------------------------------------------------------------------
# 2. EARTH ENGINE INITIALISATION
# ---------------------------------------------------------------------------

def initialise_ee(project_id):
    """Authenticate (if needed) and initialise Earth Engine."""
    try:
        ee.Initialize(project=project_id)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=project_id)


# ---------------------------------------------------------------------------
# 3. BUILD COMBINED IMAGE (Hansen + TMF)
# ---------------------------------------------------------------------------

def build_combined_image():
    """Create a single image that contains both Hansen and TMF bands."""
    hansen = ee.Image(HANSEN_ASSET).select(HANSEN_BANDS)

    # TMF Annual Changes is a tiled ImageCollection (global coverage split
    # into tiles). We need to mosaic it into a single image before merging
    # with Hansen.
    tmf = ee.ImageCollection(TMF_ASSET).select(TMF_BANDS).mosaic()

    return hansen.addBands(tmf)


# ---------------------------------------------------------------------------
# 4. CHUNKED EXTRACTION
# ---------------------------------------------------------------------------

def extract_chunk(df_chunk, image, lat_col, lon_col, scale):
    """Extract values from `image` for a subset of points."""
    features = [
        ee.Feature(
            ee.Geometry.Point([row[lon_col], row[lat_col]]),
            # keep the original DataFrame index for later join
            {"id": idx, lat_col: row[lat_col], lon_col: row[lon_col]}
        )
        for idx, row in df_chunk.iterrows()
    ]
    fc = ee.FeatureCollection(features)

    results = image.reduceRegions(
        collection=fc,
        reducer=ee.Reducer.first(),
        scale=scale
    )

    info = results.getInfo()
    return [feat["properties"] for feat in info["features"]]


def extract_all(df, image, lat_col, lon_col, chunk_size, scale, max_retries):
    """Loop over all chunks with error handling and retries."""
    all_results = []
    n_chunks = (len(df) // chunk_size) + (1 if len(df) % chunk_size else 0)

    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i + chunk_size]
        chunk_number = i // chunk_size + 1
        attempt = 0

        while attempt < max_retries:
            try:
                rows = extract_chunk(chunk, image, lat_col, lon_col, scale)
                all_results.extend(rows)
                print(f"Chunk {chunk_number}/{n_chunks} OK ({len(rows)} points)")
                break
            except Exception as e:
                attempt += 1
                print(f"Error in chunk {chunk_number}, attempt {attempt}/{max_retries}: {e}")
                time.sleep(5)
        else:
            print(f"⚠️  Chunk {chunk_number} failed after {max_retries} attempts – points skipped")

    return all_results


# ---------------------------------------------------------------------------
# 5. MAIN SCRIPT
# ---------------------------------------------------------------------------

def main():
    print("Initialising Earth Engine...")
    initialise_ee(PROJECT_ID)

    print(f"Reading input file: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    df = df.reset_index(drop=True)
    print(f"{len(df)} points loaded.")

    print("Building combined image (Hansen + TMF)...")
    combined_image = build_combined_image()

    print("Starting chunked extraction...")
    results = extract_all(
        df,
        combined_image,
        LAT_COL,
        LON_COL,
        CHUNK_SIZE,
        SCALE,
        MAX_RETRIES
    )

    result_df = pd.DataFrame(results)
    print(f"\nExtraction completed: {result_df.shape[0]} rows retrieved "
          f"out of {len(df)} original points.")

    # Join with the original DataFrame (lat/lon + any other attributes)
    if "id" in result_df.columns:
        result_df = result_df.set_index("id")
        final_df = df.join(result_df, how="left")
    else:
        final_df = result_df

    final_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Results saved to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------
# ALTERNATIVE: export via Google Drive (more robust for very large volumes,
# does not depend on local connection stability).
# ---------------------------------------------------------------------

initialise_ee(PROJECT_ID)
df = pd.read_csv(INPUT_CSV).reset_index(drop=True)
combined_image = build_combined_image()

features = [
    ee.Feature(
        ee.Geometry.Point([row[LON_COL], row[LAT_COL]]),
        {"id": idx, LAT_COL: row[LAT_COL], LON_COL: row[LON_COL]}
    )
    for idx, row in df.iterrows()
]
fc = ee.FeatureCollection(features)

results = combined_image.reduceRegions(fc, ee.Reducer.first(), scale=SCALE)

task = ee.batch.Export.table.toDrive(
    collection=results,
    description="extract_gee_gfc_tmf",
    fileFormat="CSV"
)
task.start()

while task.active():
    print("Running...", task.status()["state"])
    time.sleep(30)
print("Finished:", task.status())

# End
