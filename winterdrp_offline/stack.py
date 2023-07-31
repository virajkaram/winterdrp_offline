import os
from pathlib import Path
from typing import Optional

import numpy as np

import subprocess


def run_swarp(
    stack_paths_list: list[str] | list[Path],
    stack_paths_weight_list: list[str] | list[Path],
    swarp_config_path: str | Path,
    out_path: str | Path,
    weight_out_path: Optional[str | Path] = None,
    pixscale: Optional[float] = None,
    x_imgpixsize: Optional[float] = None,
    y_imgpixsize: Optional[float] = None,
    propogate_headerlist: Optional[list] = None,
    center_ra: Optional[float] = None,
    center_dec: Optional[float] = None,
    combine: bool = True,
    gain: Optional[float] = None,
    subtract_bkg: bool = False,
    flux_scaling_keyword: str = None,
    cache: bool = False,
    center_type: str = "ALL",
    timeout: int = 60
):
    """
    Wrapper to resample and stack images with swarp

    Parameters
    ----------
    stack_list_path : string
        Name of file containing the names of files to be stacked
        One file name per line
    weight_list_path : string
        Name of file containing the names of weight files to be stacked
        One file name per line
    swarp_config_path: str
        Path of Swarp config file
    out_path : string
        Path of stacked output file
    weight_out_path: str
        Path of output weight image
    pixscale: float
        Pixelscale in degrees
    x_imgpixsize: float
        X-dimension in pixels
    y_imgpixsize: float
        Y-dimension in pixels
    propogate_headerlist: list
        Headerlist to propagate from header
    center_ra: float
        Central RA
    center_dec: float
        Central Dec
    combine: bool
        Combine and coadd all images? For reasons internal to Swarp, it is strongly
        advised to always set this to True (even if you are running Swarp on only
        one image).
    gain: float
        Gain
    subtract_bkg: bool
        Background subtraction
    flux_scaling_keyword: str
        What flux scaling keyword do you want to use? If None, the default value in
        the config will be used
    """

    print(f"Stacking {len(stack_paths_list)} images.")
    assert len(stack_paths_list) == len(stack_paths_weight_list)
    swarp_command = (
        f"swarp -c {swarp_config_path} "
        f"{','.join(stack_paths_list)} "
        f"-IMAGEOUT_NAME {out_path} "
        f"-RESAMPLE Y -RESAMPLE_DIR {os.path.dirname(out_path)} "
    )

    if subtract_bkg:
        swarp_command += "-SUBTRACT_BACK Y "
    else:
        swarp_command += "-SUBTRACT_BACK N "
    if combine:
        swarp_command += "-COMBINE Y -COMBINE_TYPE MEDIAN "
    else:
        swarp_command += "-COMBINE N "

    if stack_paths_weight_list is not None:
        swarp_command += f" -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE " \
                         f"{','.join(stack_paths_weight_list)} "

    if weight_out_path is not None:
        swarp_command += f" -WEIGHTOUT_NAME {weight_out_path}"

    if pixscale is not None:
        swarp_command += f" -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE {pixscale}"
    else:
        swarp_command += f" -PIXELSCALE_TYPE MEDIAN"

    if propogate_headerlist is not None:
        swarp_command += " -COPY_KEYWORDS "
        for keyword in propogate_headerlist:
            swarp_command += f"{keyword},"

        # remove final comma
        swarp_command = swarp_command[:-1]

    if np.logical_and(center_ra is not None, center_dec is not None):
        swarp_command += f" -CENTER_TYPE MANUAL -CENTER {center_ra},{center_dec}"

    else:
        swarp_command += f" -CENTER_TYPE {center_type}"

    if x_imgpixsize is not None:
        swarp_command += f" -IMAGE_SIZE {x_imgpixsize}"
        if y_imgpixsize is not None:
            swarp_command += f",{y_imgpixsize}"

    if gain is not None:
        swarp_command += f" -GAIN {gain}"

    if flux_scaling_keyword is not None:
        swarp_command += f" -FSCALE_KEYWORD {flux_scaling_keyword}"

    if not cache:
        swarp_command += " -DELETE_TMPFILES Y"
    else:
        swarp_command += " -DELETE_TMPFILES N"

    print(swarp_command)
    rval = subprocess.run(swarp_command, check=True, capture_output=True, shell=True,
                          timeout=timeout)
    msg = "Successfully executed command. "

    if rval.stdout.decode() != "":
        msg += f"Found the following output: {rval.stdout.decode()}"
    print(msg)

    return
