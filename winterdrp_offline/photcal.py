from winterdrp_offline.utils import get_table_from_ldac
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astropy.stats import sigma_clip, sigma_clipped_stats


def xmatch_img_ref_cats(img_cat, tmass_cat, xmatch_radius_arcsec: float = 2):
    """
    Crossmatch the image catalog with the 2MASS catalog
    :param img_cat:
    :param tmass_cat:
    :param xmatch_radius_arcsec:
    :return:
    """
    img_cat = img_cat[img_cat['MAG_AUTO']<99]
    cat_crds = SkyCoord(img_cat['ALPHAWIN_J2000'], img_cat['DELTAWIN_J2000'], unit=u.deg)
    tmass_crds = SkyCoord(tmass_cat['ra'], tmass_cat['dec'], unit=u.deg)
    id1, id2, d2d, d3d = cat_crds.search_around_sky(tmass_crds,
                                                    xmatch_radius_arcsec*u.arcsec)
    tm_matched = tmass_cat[id1]
    cat_matched = img_cat[id2]
    return cat_matched, tm_matched


def xmatch_img_ref_cat_filenames(
    img_cat_filename: str,
    tmass_cat_filename: str,
    xmatch_radius_arcsec: float = 2,
):
    """
    Crossmatch the image catalog with the 2MASS catalog
    :param img_cat_filename:
    :param tmass_cat_filename:
    :param xmatch_radius_arcsec:
    :return:
    """
    img_cat = get_table_from_ldac(img_cat_filename)
    tmass_cat = get_table_from_ldac(tmass_cat_filename)
    cat_matched, tm_matched = xmatch_img_ref_cats(img_cat, tmass_cat, xmatch_radius_arcsec)
    return cat_matched, tm_matched


def calculate_zeropoint_outlier_rejection(img_mags, ref_mags,
                                          outlier_rejection_thresholds = [3.0, 2.5, 2.0],
                                          num_stars_threshold = 10):
    offsets = np.ma.array(
        ref_mags - img_mags
    )
    zp_mean, zp_med, zp_std, num_stars = np.nan, np.nan, np.nan, 0
    for outlier_thresh in outlier_rejection_thresholds:
        cl_offset = sigma_clip(offsets, sigma=outlier_thresh)
        num_stars = np.sum(np.invert(cl_offset.mask))

        zp_mean, zp_med, zp_std = sigma_clipped_stats(
            offsets, sigma=outlier_thresh
        )

        if num_stars > num_stars_threshold:
            break

    return zp_mean, zp_std, num_stars