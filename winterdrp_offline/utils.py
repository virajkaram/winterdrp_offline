import os
from glob import glob
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from astropy.table import Table
from pathlib import Path


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


def write_weight_images(imagelist: str | Path | list):
    """
    Function to write weight images for a list of images. The weight image is
    a binary mask where 1 is a valid pixel and 0 is a bad pixel.
    :param imagelist:
    :return:
    """
    if isinstance(imagelist, str) or isinstance(imagelist, Path):
        imagelist = [imagelist]
    masklist = []
    for imgname in imagelist:
        data = fits.getdata(imgname)
        header = fits.getheader(imgname)
        mask = np.ones_like(data, dtype=float)
        mask[np.isnan(data)] = 0.0
        maskname = imgname.replace(".fits", "_weight.fits")
        fits.writeto(maskname, data=mask, header=header, overwrite=True)
        masklist.append(maskname)
    return masklist



def convert_hdu_to_ldac(hdu):
    """
    Convert an hdu table to a fits_ldac table (format used by astromatic suite)

    Parameters
    ----------
    hdu: `astropy.io.fits.BinTableHDU` or `astropy.io.fits.TableHDU`
        HDUList to convert to fits_ldac HDUList

    Returns
    -------
    tbl1: `astropy.io.fits.BinTableHDU`
        Header info for fits table (LDAC_IMHEAD)
    tbl2: `astropy.io.fits.BinTableHDU`
        Data table (LDAC_OBJECTS)
    """
    from astropy.io import fits
    import numpy as np
    tblhdr = np.array([hdu.header.tostring(',')])
    col1 = fits.Column(name='Field Header Card', array=tblhdr, format='13200A')
    cols = fits.ColDefs([col1])
    tbl1 = fits.BinTableHDU.from_columns(cols)
    tbl1.header['TDIM1'] = '(80, {0})'.format(len(hdu.header))
    tbl1.header['EXTNAME'] = 'LDAC_IMHEAD'
    tbl2 = fits.BinTableHDU(hdu.data)
    tbl2.header['EXTNAME'] = 'LDAC_OBJECTS'
    return (tbl1, tbl2)


def convert_table_to_ldac(tbl):
    """
    Convert an astropy table to a fits_ldac

    Parameters
    ----------
    tbl: `astropy.table.Table`
        Table to convert to ldac format
    Returns
    -------
    hdulist: `astropy.io.fits.HDUList`
        FITS_LDAC hdulist that can be read by astromatic software
    """
    from astropy.io import fits
    import tempfile
    f = tempfile.NamedTemporaryFile(suffix='.fits', mode='rb+')
    tbl.write(f, format='fits')
    f.seek(0)
    hdulist = fits.open(f, mode='update')
    tbl1, tbl2 = convert_hdu_to_ldac(hdulist[1])
    new_hdulist = [hdulist[0], tbl1, tbl2]
    new_hdulist = fits.HDUList(new_hdulist)
    return new_hdulist


def save_table_as_ldac(tbl, filename, **kwargs):
    """
    Save a table as a fits LDAC file

    Parameters
    ----------
    tbl: `astropy.table.Table`
        Table to save
    filename: str
        Filename to save table
    kwargs:
        Keyword arguments to pass to hdulist.writeto
    """
    hdulist = convert_table_to_ldac(tbl)
    hdulist.writeto(filename, **kwargs)


def plot_image(data, thresh=3.0, ax=None):
    """
    Plot an image with plt.imshow and threshold
    :param data:
    :param thresh:
    :return:
    """
    if ax is None:
        fig, ax = plt.subplots()

    _, med, std = sigma_clipped_stats(data, sigma=thresh)
    ax.imshow(data, vmin=med-3*std, vmax=med+3*std, cmap='gray', origin='lower')
    # plot colorbar
    cbar = plt.colorbar(ax.imshow(data, vmin=med-3*std, vmax=med+3*std, cmap='gray', origin='lower'), ax=ax)

    return ax


def write_image(image, filename, header=None, overwrite=True):
    """
    Write an image to a fits file
    :param image:
    :param filename:
    :param header:
    :param overwrite:
    :return:
    """
    # make parent directories if they don't exist
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    hdu = fits.PrimaryHDU(image, header=header)
    hdu.writeto(filename, overwrite=overwrite)


def get_table_from_ldac(filename: str | Path, frame: int = 1) -> Table:
    """
    Load an astropy table from a fits_ldac by frame (Since the ldac format has column
    info for odd tables, giving it twce as many tables as a regular fits BinTableHDU,
    match the frame of a table to its corresponding frame in the ldac file).

    Parameters
    ----------
    filename: str
        Name of the file to open
    frame: int
        Number of the frame in a regular fits file
    """
    if frame > 0:
        frame = frame * 2
    tbl = Table.read(filename, hdu=frame)
    return tbl
