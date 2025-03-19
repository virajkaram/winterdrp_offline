import os
import re

import pandas as pd


def parse_filename(fname):
    """
    Example filename:
      pc_2_60nA_exp_0.06_light_7.66nA_test_20240801-173954-563.fits

    We want to parse out:
      address = "pc"
      dark_photocurrent = 2.60  (from "2_60" -> "2.60")
      exposure_time = 0.06
      frame_type = "light"
      light_photocurrent = 7.66
      datetime_local = "20240801-173954-563"

    Returns a dict of parsed information or None if it doesn't match.
    """
    pattern = re.compile(
        r"^(?P<address>[a-z]{2})_"  # 2-char address, e.g. pc
        r"(?P<dc_str>\d+_\d+)nA_exp_"  # e.g. 2_60nA_exp_
        r"(?P<exposure_time>[\d.]+)_"  # e.g. 0.06_
        r"(?P<frame_type>(light|dark))_"  # light or dark
        r"(?P<light_photocurrent>[\d.]+)nA_"  # e.g. 7.66nA_
        r".*_"  # skip intermediate text
        r"(?P<datetime_local>[^.]+)\.fits$"  # up to the .fits
    )
    match = pattern.match(fname)
    if not match:
        return None

    # Convert the dark photocurrent by replacing the underscore with a decimal
    dc_str = match["dc_str"]  # e.g. "2_60" -> "2.60"
    dark_photocurrent = float(dc_str.replace("_", "."))

    return {
        "address": match["address"],
        "dark_photocurrent": dark_photocurrent,
        "exposure_time": float(match["exposure_time"]),
        "frame_type": match["frame_type"],
        "light_photocurrent": float(match["light_photocurrent"]),
        "datetime_local": match["datetime_local"],
        "filename": fname,  # the bare filename
    }


def load_fits_directory(dirname):
    """
    Given a directory of FITS files, parse each using parse_filename().
    Returns a DataFrame of the parsed fields + 'filepath' column
    (the full path to each file).
    """
    rows = []
    for fname in os.listdir(dirname):
        if not fname.lower().endswith(".fits"):
            continue
        parsed = parse_filename(fname)
        if parsed is not None:
            # Attach the full file path for convenience
            parsed["filepath"] = os.path.join(dirname, fname)
            rows.append(parsed)
    return pd.DataFrame(rows)


def load_light_and_dark_frames(light_dir, dark_dir):
    """
    Loads FITS files from two separate directories (light_dir, dark_dir).
    Returns a single DataFrame containing both sets, parsed via parse_filename(),
    including a 'filepath' column for each file.
    """
    df_light = load_fits_directory(light_dir)
    df_dark = load_fits_directory(dark_dir)
    return pd.concat([df_light, df_dark], ignore_index=True)


def filter_exposures_with_both_light_and_dark(df):
    """
    Keep only rows whose exposure_times have both 'light' and 'dark' frames present.
    """
    # For each exposure_time, figure out the set of frame_types
    grouped = df.groupby("exposure_time")["frame_type"].apply(set)

    # We only want those exposure_times that contain both 'light' and 'dark'
    valid_exposures = grouped[
        grouped.apply(lambda s: {"light", "dark"}.issubset(s))
    ].index

    # Filter the original df to keep only those exposure times
    return df[df["exposure_time"].isin(valid_exposures)].copy()


def get_light_frames_for_exposure(df, exposure):
    """
    Return all rows (as a DataFrame) from df for a given exposure_time
    with frame_type == 'light'. Includes 'filepath' for opening the file.
    """
    return df[(df["exposure_time"] == exposure) & (df["frame_type"] == "light")]


def get_dark_frames_for_exposure(df, exposure):
    """
    Return all rows (as a DataFrame) from df for a given exposure_time
    with frame_type == 'dark'. Includes 'filepath' for opening the file.
    """
    return df[(df["exposure_time"] == exposure) & (df["frame_type"] == "dark")]


def get_unique_exposures_in_ascending_order(df):
    """
    Return a list of all unique exposure times in ascending order.
    """
    exposures = [float(x) for x in sorted(df["exposure_time"].unique())]
    return exposures
