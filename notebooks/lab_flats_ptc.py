import os
import re
from glob import glob
from pathlib import Path

import numpy as np
from astropy.io import fits

from winterdrp_offline.image_operations import (
    median_combine_images,
    merge_channels_into_data,
    normalize_and_median_combine_images,
    split_data_into_channels,
)
from winterdrp_offline.labflats import (
    filter_exposures_with_both_light_and_dark,
    get_dark_frames_for_exposure,
    get_light_frames_for_exposure,
    get_unique_exposures_in_ascending_order,
    load_fits_directory,
    load_light_and_dark_frames,
    parse_filename,
)
from winterdrp_offline.utils import get_table_from_ldac, plot_image, write_image

# plan:
# 1. read in the light frames at each exposure time, and save them as multi-extension fits files
# 2. make median combined light frames at each exposure time and save them as fits images
# 3. read in the dark frames at each exposure time, and save them as multi-extension fits files
# 4. make median combined dark frames at each exposure time and save them as fits images
# 5. make a photon transfer curve by plotting the median light frame against the variance of the light frames

# 1. read in the light frames at each exposure time, and save them as multi-extension fits files
# read in the light frames

base_dir = Path(
    os.path.join(
        os.getenv("HOME"), "data", "calibration_hackathon", "hackathon_datasets"
    )
)


flat_dir = base_dir.joinpath(f"flats/lab_flats/pc_2_60nA_post_mods_ptc-coarse")
flat_light_dir = flat_dir.joinpath("flat/raw")
flat_dark_dir = flat_dir.joinpath("dark/raw")

df_all = load_light_and_dark_frames(flat_light_dir, flat_dark_dir)
print("Combined (light & dark) frames:")
print(df_all)
df_filtered = filter_exposures_with_both_light_and_dark(df_all)
print("\nFiltered to exposures with both light and dark frames:")
print(df_filtered)

exposures_sorted = get_unique_exposures_in_ascending_order(df_filtered)
print("\nUnique exposures in ascending order (that have both light and dark):")
print(exposures_sorted)

# now go through each exposure time and get the light and dark frames, and make the median combined light and dark frames
for exp_time in exposures_sorted:
    # Example of retrieving specific frames:
    print(f"\nLight frames at exposure_time={exp_time}:")
    light_frames = get_light_frames_for_exposure(df_filtered, exp_time)
    dark_frames = get_dark_frames_for_exposure(df_filtered, exp_time)

    median_light_frame = median_combine_images(light_frames["filepath"])
    median_dark_frame = median_combine_images(dark_frames["filepath"])

    # plot the median light and dark frames
    ax = plot_image(median_light_frame)
    ax.set_title(f"Median light frame at {exp_time} seconds")
    ax = plot_image(median_dark_frame)
    ax.set_title(f"Median dark frame at {exp_time} seconds")


"""filename = "pc_2_60nA_exp_0.06_light_7.66nA_test_20240801-173954-563.fits"
print(parse_filename(filename))"""


"""flat_lights = glob(str(flat_dir.joinpath("flat/raw/*exp_0.06*fits")))
print(f"Found {len(flat_lights)} flat lights in {flat_dir}")
flat_darks = glob(str(flat_dir.joinpath("dark/raw/*exp_0.06*fits")))
print(f"Found {len(flat_darks)} flat darks in {flat_dir}")
"""
