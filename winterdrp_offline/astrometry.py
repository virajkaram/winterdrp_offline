"""
Script containing functions to run astrometry.net locally
"""
import subprocess
import os
from pathlib import Path
from typing import Optional

from astropy.io import fits


def run_astrometry_net(imagelist: str | list, output_dir: str,
                       masklist: str | list = None,
                       *args, **kwargs):
    """
    function to execute `run_astrometry_net_single` on several images in batch
    """
    if not isinstance(imagelist, list):
        imagelist = [imagelist]

    # make output directory if it doesn't exist

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    solved_images, failed_images = [], []
    for ind, img in enumerate(imagelist):
        mask_path = None
        if masklist is not None:
            mask_path = masklist[ind]
        return_imgname = run_astrometry_net_single(img, Path(output_dir),
                                                   mask_path=mask_path,
                                                   *args, **kwargs)
        if return_imgname is not None:
            solved_images.append(return_imgname)
        else:
            failed_images.append(return_imgname)
    return solved_images, failed_images


def run_astrometry_net_single(
        img_path: str | Path,
        output_dir: str | Path,
        mask_path: str | Path = None,
        scale_bounds: Optional[tuple | list] = [15, 23],
        # limits on scale (lower, upper)
        scale_units: Optional[str] = 'amw',  # scale units ('degw', 'amw')
        downsample: Optional[float | int] = None,  # downsample by factor of __
        timeout: Optional[float] = 30.0,  # astrometry cmd execute timeout, in seconds
        use_sextractor: bool = False,
        sextractor_path: str = "sex",
        search_radius_deg: float = 1.0,
        parity: str = None,
        sextractor_config_path: str = None,
        x_image_key: str = "X_IMAGE",
        y_image_key: str = "Y_IMAGE",
        sort_key_name: str = "MAG_AUTO",
):
    """
    function to run astrometry.net locally on one image, with options to adjust settings
    default: solve-field <img> -D <output_dir> -N <newname> -O
    """
    # name for new file if a-net solves (otherwise a-net writes to '<img>.new')
    newname = output_dir.joinpath(Path(os.path.basename(img_path).replace(".fits",
                                                                          ".solved.fits"
                                                                          )))
    basename = os.path.basename(img_path).split(".fits")[0]

    # run a-net (solve-field)
    cmd = (
        f"solve-field {img_path} "
        f"--dir {output_dir} "
        f"--new-fits {newname} "
        f"--overwrite "
        f"--out {basename} "  # use this base name for outputs (instead of 'temp_...')
    )

    if scale_bounds is not None:
        cmd += f" --scale-high {max(scale_bounds)} "
        cmd += f" --scale-low {min(scale_bounds)} "

    if scale_units is not None:
        cmd += f"--scale-units {scale_units} "

    if downsample is not None:
        cmd += f"--downsample {downsample} "

    # cmd with a ra, dec first guess (speeds up solution)
    header = fits.open(img_path)[0].header  # pylint: disable=no-member
    ra_req, dec_req = header["RADEG"], header["DECDEG"]  # requested ra, dec
    if use_sextractor:
        if mask_path is not None:
            sextractor_path += f" -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE {mask_path} "
        cmd += f"--use-source-extractor --source-extractor-path '{sextractor_path}' "

    if sextractor_config_path is not None:
        cmd += f"--source-extractor-config {sextractor_config_path} "

    cmd += f"-X {x_image_key} -Y {y_image_key} -s {sort_key_name} --sort-ascending "

    if parity is not None:
        assert parity in ["pos", "neg"]
        cmd += f"--parity {parity} "

    cmd_loc = (
            cmd + f"--ra {ra_req} --dec {dec_req} --radius {search_radius_deg} "
    )  # radius takes on units of ra, dec

    return_imgname = None
    try:
        print(
            f"Running a-net with ra,dec guess and timeout {timeout}. \n"
            f"A-net command:\n {cmd_loc}"
        )

        rval = subprocess.run(
            cmd_loc, check=True, capture_output=True, shell=True, timeout=timeout
        )

        msg = "Successfully executed command. "

        if rval.stdout.decode() != "":
            msg += f"Found the following output: {rval.stdout.decode()}"
        print(msg)
        return_imgname = newname.as_posix()

    except Exception as err:
        print(f"Error running a-net with ra,dec guess: {err}")

    return return_imgname
