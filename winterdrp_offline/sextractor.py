import subprocess
from pathlib import Path
from winter_utils.ldactools import get_table_from_ldac, save_table_as_ldac, \
    convert_table_to_ldac
from astropy.io import fits
import numpy as np

astrom_sex = Path(__file__).parent.joinpath("config/astrom.sex")
astrom_param = Path(__file__).parent.joinpath("config/astrom.param")
astrom_filter = Path(__file__).parent.joinpath("config/default.conv")
astrom_nnw = Path(__file__).parent.joinpath("config/default.nnw")


def run_sextractor(
        imgname: str | Path,
        pixscale: float = 1.12,
        regions: bool = True,
        weightimg: str = 'weight.fits'
):
    """
    Run sextractor on the proc image file

    Parameters
    ----------
    imgname: Name of the image to run sextractor on
    pixscale: Pixel scale of the image
    regions: boolean whether to write out regions
    weightimg: Weight image to use

    Returns
    -------
    None

    """
    image_basepath = imgname.replace('.fits', '')
    try:
        command = f'sex -c {astrom_sex} {imgname} -CATALOG_NAME {imgname}.cat ' + \
                  f'-CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME {astrom_param} ' + \
                  f'-FILTER_NAME {astrom_filter} -STARNNW_NAME {astrom_nnw} ' \
                  f'-PIXEL_SCALE {pixscale} -DETECT_THRESH 3 -ANALYSIS_THRESH 3 ' \
                  f'-SATUR_LEVEL 60000 -WEIGHT_TYPE MAP_WEIGHT ' \
                  f'-WEIGHT_IMAGE {weightimg} -CHECKIMAGE_TYPE BACKGROUND,-BACKGROUND '\
                  f'-CHECKIMAGE_NAME {image_basepath}.back.fits,' \
                  f'{image_basepath}.bkgsub.fits'

        print('Executing command : %s' % (command))
        subprocess.run(command.split(), check=True, capture_output=True)

    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with error %s.' % (err))
        return

    if regions:
        t = get_table_from_ldac(imgname + '.cat')

        with open(imgname + '.cat' + '.stats.reg', 'w') as f:
            f.write('image\n')
            for row in t:
                f.write('CIRCLE(%s,%s,%s) # text={%.2f}\n' % (
                    row['X_IMAGE'], row['Y_IMAGE'], row['FWHM_IMAGE'] / 2,
                    row['FWHM_IMAGE']))

        clean_catalog = t[(t['FLAGS'] == 0) & (t['FWHM_IMAGE'] > 0) &
                          (t['SNR_WIN'] > 0)]
        clean_hdu = convert_table_to_ldac(clean_catalog)
        with fits.open(imgname + '.cat') as hdul:
            clean_hdulist = fits.HDUList([hdul[0], hdul[1], clean_hdu[2]])
            clean_hdulist.writeto(imgname + '.clean.cat', overwrite=True)

        with open(imgname + '.cat' + '.clean.stats.reg', 'w') as f:
            f.write('image\n')
            for row in clean_catalog:
                f.write('CIRCLE(%s,%s,%s) # text={%.2f}\n' % (
                    row['X_IMAGE'], row['Y_IMAGE'], row['FWHM_IMAGE'] / 2,
                    row['FWHM_IMAGE']))


def mask_sextractor_skysub_pixels(imgname, skysub_name):
    with fits.open(skysub_name, 'update') as hdul:
        raw_data = fits.getdata(imgname)
        mask = np.isnan(raw_data)
        hdul[0].data[mask] = np.nan


def subtract_sky(imagelist: list[str],
                 write_sky_model: bool = True) -> list[str]:
    """
    Function to subtract sky from a list of images

    :param imagelist: 2D array with boardids as first index and MEF image as second
    index
    :return: list of sky subtracted image paths
    """
    if not isinstance(imagelist, np.ndarray):
        imagelist = np.array(imagelist)
    sky_subtracted_paths = []

    board_image_data = [fits.getdata(x) for x in imagelist]
    board_image_data = np.array(board_image_data)
    print(f"Median combining {len(imagelist)} images to make sky model")
    sky_model = np.nanmedian(board_image_data, axis=0)
    sky_model = np.array(sky_model, dtype=np.float32)

    if write_sky_model:
        fits.writeto(f"{imagelist[0].replace('.fits', '')}_sky_model.fits",
                     sky_model, overwrite=True)

    for count, image in enumerate(imagelist):
        sky_sub_header = fits.getheader(image)
        sky_sub_header['SATURATE'] -= np.nanmedian(sky_model)
        sky_sub_path = image.replace('.fits', '_sky.fits')
        image_data = fits.getdata(image)
        sky_sub_data = image_data - \
                       sky_model * np.nanmedian(image_data) / np.nanmedian(sky_model)
        fits.writeto(sky_sub_path, sky_sub_data, fits.getheader(image),
                     overwrite=True)
        sky_subtracted_paths.append(sky_sub_path)

    return sky_subtracted_paths
