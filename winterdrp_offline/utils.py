import os
from glob import glob
from pathlib import Path
from astropy.io import fits
import numpy as np


def subtract_dark(image: np.ndarray, masterdark: np.ndarray):
    """
    Subtract a master dark from an image

    :param image: image to subtract dark from
    :param masterdark: master dark to subtract
    :return: image with dark subtracted
    """
    return image - masterdark


def combine_images(imagepaths):
    data = np.array([fits.getdata(x) for x in imagepaths])
    combined_data = np.nanmedian(data, axis=0)
    return combined_data


def get_imagelist_from_directory(dir, select_type="*.fits"):
    imglist = glob(os.path.join(dir, select_type))
    nonmask_imglist = [x for x in imglist if (("mask" not in x) & ("back" not in x))]
    return nonmask_imglist


def copy_files_to_directory(filelist, dir):
    Path(dir).mkdir(parents=True, exist_ok=True)
    for filename in filelist:
        os.system(f"cp {filename} {dir}")


def write_mask(imagelist):
    masklist = []
    for imgname in imagelist:
        data = fits.getdata(imgname)
        header = fits.getheader(imgname)
        mask = np.ones_like(data, dtype=float)
        mask[np.isnan(data)] = 0.0
        maskname = imgname.replace(".fits", "_mask.fits")
        fits.writeto(maskname, data=mask, header=header, overwrite=True)
        masklist.append(maskname)
    return masklist
