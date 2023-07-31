from astroquery.gaia import Gaia
import os
from winter_utils.ldactools import save_table_as_ldac
from sextractor import run_sextractor
import subprocess
from pathlib import Path
from astropy.io import fits
import numpy as np
from utils import write_mask

astrom_scamp = Path(__file__).parent.joinpath("config/astrom.scamp")


def make_gaia_catalog(ra, dec, tmcatname, catalog_box_size_arcmin, catalog_min_mag,
                      catalog_max_mag, writeldac=True, write_regions=True):
    print(
        "SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.tmass_best_neighbour AS "
        "tbest, gaiadr1.tmass_original_valid AS tmass WHERE g.source_id = "
        "tbest.source_id AND tbest.tmass_oid = tmass.tmass_oid AND CONTAINS(POINT("
        "'ICRS', g.ra, g.dec), CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND tmass.j_m > "
        "%.2f AND tmass.j_m < %.2f AND tbest.number_of_mates=0 AND "
        "tbest.number_of_neighbours=1;" % (
            ra, dec, catalog_box_size_arcmin / 60, catalog_min_mag, catalog_max_mag))
    try:
        job = Gaia.launch_job_async(
            "SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.tmass_best_neighbour AS "
            "tbest, gaiadr1.tmass_original_valid AS tmass WHERE g.source_id = "
            "tbest.source_id AND tbest.tmass_oid = tmass.tmass_oid AND CONTAINS("
            "POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND "
            "tmass.j_m > %.2f AND tmass.j_m < %.2f AND tbest.number_of_mates=0 AND "
            "tbest.number_of_neighbours=1;" % (
                ra, dec, catalog_box_size_arcmin / 60, catalog_min_mag,
                catalog_max_mag),
            dump_to_file=False)
        print('Yay')
        t = job.get_results()
    except:
        print('Could not query Gaia .. ')
        return 0

    # print(t.colnames)
    # t.remove_columns(['designation', 'phot_variable_flag', 'datalink_url',
    # 'epoch_photometry_url', 'original_ext_source_id', 'designation_2'])
    t.remove_columns(
        ['designation', 'phot_variable_flag', 'datalink_url', 'original_ext_source_id',
         'DESIGNATION'])
    t['ph_qual'] = t['ph_qual'].astype(str)
    t['ra_errdeg'] = t['ra_error'] / 3.6e6
    t['dec_errdeg'] = t['dec_error'] / 3.6e6
    t['FLAGS'] = 0

    if writeldac:
        if os.path.exists(tmcatname):
            os.remove(tmcatname)
        save_table_as_ldac(t, tmcatname)

    if write_regions:
        with open(tmcatname + '.reg', 'w') as f:
            f.write("wcs\n")
            for i in range(len(t)):
                f.write("circle(%.7f,%.7f,%.7f) # color=red\n" % (
                    t['ra'][i], t['dec'][i], 0.0005))


def run_scamp(filelist: list | str, output_dir: str, write_file: bool = True):
    if isinstance(filelist, str):
        filelist = [filelist]

    masklist = write_mask(filelist)
    catnames = []
    for ind, file in enumerate(filelist):
        run_sextractor(file, weightimg=masklist[ind])
        catnames.append(file + '.clean.cat')

    scamp_filename = Path(output_dir).joinpath('scamp.txt')
    with open(scamp_filename, 'w') as f:
        for catalog_name in catnames:
            f.write(catalog_name + '\n')

    center_ras = [fits.getheader(x)['CRVAL1'] for x in filelist]
    center_decs = [fits.getheader(x)['CRVAL2'] for x in filelist]
    center_ra = np.median(center_ras)
    center_dec = np.median(center_decs)

    make_gaia_catalog(ra=center_ra, dec=center_dec,
                      tmcatname=output_dir + '/gaia.cat',
                      catalog_box_size_arcmin=15,
                      catalog_min_mag=7,
                      catalog_max_mag=20)

    print('Downloaded gaia sources')
    run_scamp_command(scamp_filename, catalog_name = output_dir + '/gaia.cat')

    if write_file:
        for imgname in filelist:
            write_scamped_file(imgname, imgname + '.clean.head')


def run_scamp_command(scamp_filename, catalog_name):
    # Run scamp on the proc image file
    try:
        command = 'scamp -c ' + astrom_scamp.as_posix() + f' @{scamp_filename}' \
                  + f' -ASTREFCAT_NAME {catalog_name}' \
                  + ' -ASTREFMAG_KEY j_m'

        print('Executing command : %s' % command)
        rval = subprocess.run(command.split(), check=True, capture_output=True)
        print('Process completed')
        print(rval.stdout.decode())

    except subprocess.CalledProcessError as err:
        print('Could not run scamp with error %s.' % (err))


def write_scamped_file(imgname, headername):
    img = fits.open(imgname)
    header = img[0].header
    data = img[0].data
    img.close()

    with open(headername, 'r') as f:
        a = f.read()

    # Remove any existing astrometry keywords
    astrometry_keys = get_astrometry_keys()
    for k in astrometry_keys:
        if k in header:
            del header[k]

    h = fits.Header()
    h = h.fromstring(a, sep='\n')
    for k in h.keys():
        if k == 'HISTORY' or k == 'COMMENT':
            continue

        header[k] = h[k]
    header.add_history('%s' % (h['HISTORY']))
    procHDU = fits.PrimaryHDU(data)
    procHDU.header = header
    procHDU.writeto(imgname.replace('.fits', '') + '.scampastrom.fits', overwrite=True)


def get_astrometry_keys() -> list:
    """
    Function to get a list of common astrometric keywords that could be present in a
    fits header
    Returns:

    """
    # List for all astrometric keywords that could go in a header
    # First add basic keywords that could be present in all wcs headers
    astrometric_keywords = [
        "CTYPE1",
        "CTYPE2",
        "CRVAL1",
        "CRVAL2",
        "CRPIX1",
        "CRPIX2",
        "CD1_1",
        "CD1_2",
        "CD2_1",
        "CD2_2",
        "CDELT1",
        "CDELT2",
        "PC1_1",
        "PC1_2",
        "PC2_1",
        "PC2_2",
        "PC001001",
        "PC002001",
        "PC001002",
        "PC002002",
    ]
    # Add TPV/ZPN distortion keywords -
    # https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html
    for i in range(40):
        astrometric_keywords.append(f"PV1_{i}")
        astrometric_keywords.append(f"PV2_{i}")

    # Add SIP distortion keywords, upto order 10
    astrometric_keywords.append("A_ORDER")
    astrometric_keywords.append("B_ORDER")
    astrometric_keywords.append("AP_ORDER")
    astrometric_keywords.append("BP_ORDER")
    for i in range(10):
        for j in range(10):
            astrometric_keywords.append(f"A_{i}_{j}")
            astrometric_keywords.append(f"AP_{i}_{j}")
            astrometric_keywords.append(f"B_{i}_{j}")
            astrometric_keywords.append(f"BP_{i}_{j}")

    astrometric_keywords += [
        "LONPOLE",
        "LATPOLE",
        "CUNIT1",
        "CUNIT2",
        "IMAGEW",
        "IMAGEH",
        "WCSAXES",
        "EQUINOX",
    ]
    # Add old style WCS keywords
    for i in range(40):
        astrometric_keywords.append(f"PROJP{i}")

    # Add some SWARP-specific keywords that can come from Scamp

    astrometric_keywords.append("FLXSCALE")
    return astrometric_keywords
