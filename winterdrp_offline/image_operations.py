from typing import Any

import numpy as np
import numpy.typing as npt
from astropy.io import fits


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
            data_list.append(data / np.nanmedian(data))

    data_array = np.array(data_list)
    return np.nanmedian(data_array, axis=0)


def split_data_into_channels(
    data: npt.NDArray[Any],
) -> npt.NDArray[Any]:
    """
    Splits a 2D array into 8 channels by sampling rows and columns
    at regular offsets.

    Each of the 8 channels is extracted by:
      1. Choosing a row offset (j::2), where j = 1 for the first four channels
         (channels 0..3) and j = 0 for the latter four (channels 4..7).
      2. Choosing a column offset ((3 - i) % 4::4), where i is the channel index.

    Parameters
    ----------
    data : numpy.typing.NDArray[Any]
        2D array of shape (height, width). Must be evenly divisible so that
        height % 2 == 0 and width % 4 == 0.

    Returns
    -------
    channels_3d : numpy.typing.NDArray[Any]
        3D array of shape (8, height//2, width//4). The first dimension indexes
        the 8 channels, and the remaining two dimensions are the downsampled rows
        and columns for each channel.
    """
    height, width = data.shape

    # Number of channels to produce
    channels = 8

    # Prepare the output array
    data_8ch = np.zeros((channels, height // 2, width // 4), dtype=data.dtype)

    # Fill each of the 8 channels
    for i in range(channels):
        # Row offset (use j=1 for channels 0..3, j=0 for channels 4..7)
        j = 1 if (3 - i) >= 0 else 0
        # Column offset is (3 - i) % 4
        data_8ch[i] = data[j::2, (3 - i) % 4 :: 4]

    return data_8ch


def merge_channels_into_data(data_8ch):
    """
    Inverts the `split_data_into_channels()` step, returning data with shape:
        (N, 2 * data_8ch.shape[2], 4 * data_8ch.shape[3])
    which corresponds to the cropped region that was originally split.

    Parameters
    ----------
    data_8ch : np.ndarray
        Array of shape (8, N, h, w) returned by `split_data_into_channels`.

    Returns
    -------
    data : np.ndarray
        Reconstructed array with shape (N, 2*h, 4*w), which is the merged
        version of the 8-channel data.
    """
    # data_8ch has shape: (8, N, h, w)
    # The final output should have shape: (N, 2*h, 4*w)
    _, N, h, w = data_8ch.shape
    out = np.zeros((N, 2 * h, 4 * w), dtype=data_8ch.dtype)

    for i in range(8):
        # This is the same indexing logic used in the original splitting:
        #   j controls which set of rows (odd/even)
        #   (3 - i) % 4 controls the column offset
        j = 1 if (3 - i) >= 0 else 0
        out[:, j::2, (3 - i) % 4 :: 4] = data_8ch[i]

    return out


def package_image_list_into_mef(imglist: list, output_mef: str):
    """
    Function to package a list of images into a multi-extension fits file
    :param imglist:
    :param output_mef:
    :return:
    """
    data_list = []
    for imgname in imglist:
        with fits.open(imgname) as hdul:
            data = hdul[0].data
            data_list.append(data)

    data_array = np.array(data_list)
    hdu = fits.PrimaryHDU(data_array)
    hdu.writeto(output_mef, overwrite=True)
    return output_mef


def package_image_data_into_mef(data_list: list, output_mef: str):
    """
    Function to package a list of images into a multi-extension fits file
    :param data_list:
    :param output_mef:
    :return:
    """
    data_array = np.array(data_list)
    hdu = fits.PrimaryHDU(data_array)
    hdu.writeto(output_mef, overwrite=True)
    return output_mef
