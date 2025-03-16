from astropy.io import fits
import numpy as np


def median_combine_images(imglist: list):
    """
    Function to return median combined data
    :param imglist:
    :return:
    """
    data_list = []
    for imgname in imglist:
        with fits.open(imgname) as hdul:
            data = hdul[0].data
            data_list.append(data)

    data_array = np.array(data_list)
    return np.nanmedian(data_array, axis=0)



def normalize_and_median_combine_images(imglist: list):
    """
    Function to normalize and median combine images
    :param imglist:
    :return:
    """
    data_list = []
    for imgname in imglist:
        with fits.open(imgname) as hdul:
            data = hdul[0].data
            data_list.append(data/np.nanmedian(data))

    data_array = np.array(data_list)
    return np.nanmedian(data_array, axis=0)
