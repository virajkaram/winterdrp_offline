from astropy.io import ascii
from astropy.time import Time
import os
from astrometry import run_astrometry_net
from stack import run_swarp
from pathlib import Path
from sextractor import run_sextractor, mask_sextractor_skysub_pixels, subtract_sky
from scamp import make_gaia_catalog, run_scamp_command, run_scamp
from astropy.io import fits
import numpy as np
from utils import get_imagelist_from_directory, copy_files_to_directory, write_mask, \
    combine_images, subtract_dark


if __name__ == "__main__":
    ########### User-specifiable options ##############
    night = "20230716"
    base_dir = f"/Users/viraj/winter_data/winter/{night}/" # base directory
                                                           # where the raw data is stored
    unpack_dir = os.path.join(base_dir,
                              "raw_unpacked")  # directory where the unpacked
                                               # raw images are stored
    obslog_file = os.path.join(base_dir, f"observation_log_{night}.csv")
    # path to the observation log file obtained from the winter server.

    save_dark_sub = True  # Set to True if you want to save the dark-subtracted images
    skip_making_cals = True  # Set to True if you have already made
    # the masterdark  and masterflat
    subids_to_reduce = ["0_0"]  # a combination of 0_0 and 0_1

    start_step = "astrometry"
    all_steps = ["astrometry"]  # subset of ["detrend", "astrometry", "stack"]

    exptime = 120.0
    board_id = 4
    filter = "J"
    fieldid = 7443  # field of the science observation that is being reduced.

    flat_exptime = 120  # exptime of the images to use for flat-fielding

    astrometry_type = "astrometry_net"
    bkg_sub_type = "sextractor"

    ############################### END OF USER-SPECIFIABLE OPTIONS ###################

    masterdark_dir = os.path.join(base_dir, "masterdark")
    dark_sub_dir = os.path.join(base_dir, "dark_sub")
    detrend_dir = os.path.join(base_dir, "detrend")
    astrometry_dir = os.path.join(base_dir, "astrometry")
    stack_dir = os.path.join(base_dir, "stack")

    astrom_anet_sex = os.path.join(os.path.dirname(__file__), "config/astrom_anet.sex")
    swarp_config = os.path.join(os.path.dirname(__file__), "config/config.swarp")

    obs_log = ascii.read(obslog_file)
    t_sunset = Time('2023-07-14T04:00:00') # only use dark images after sunset

    steps_to_perform = all_steps[all_steps.index(start_step):]

    ###########################################
    # Start processing #########################
    for sub_id in subids_to_reduce:
        full_pathnames = [os.path.join(unpack_dir,
                                       x['BASENAME'].split('.fits')[0]
                                       + f"_{board_id}_{sub_id}.fits")
                          for x in obs_log]

        obs_log['PATHNAME'] = full_pathnames

        if not skip_making_cals:
            # Make master darks and flats
            uniq_exptimes = np.unique(obs_log['EXPTIME'])
            obstimes = Time(obs_log['UTCISO'])

            for exp_time in uniq_exptimes:
                darks = obs_log[(obs_log['FILTERID'] == "dark") &
                                (obs_log["BOARD_ID"] == board_id) &
                                (obs_log["EXPTIME"] == exp_time) &
                                (obstimes > t_sunset)
                                ]

                if len(darks) == 0:
                    print(f"No darks found for exptime={exp_time}")
                    continue
                print(f"Using {len(darks)} images to make a dark.")

                master_dark_data = combine_images(list(darks['PATHNAME']))

                master_dark_header = fits.getheader(darks['PATHNAME'][0])
                masterdark_name = os.path.join(masterdark_dir,
                                               f"masterdark_{exp_time}_{board_id}"
                                               f"_{sub_id}.fits")
                if not os.path.exists(masterdark_dir):
                    os.makedirs(masterdark_dir)

                fits.writeto(masterdark_name,
                             data=master_dark_data, header=master_dark_header,
                             overwrite=True)
                print(f"Wrote master dark to {masterdark_name}")

            # Select files to make a flat image.
            flats = obs_log[(obs_log['FILTERID'] == filter) &
                            (obs_log["BOARD_ID"] == board_id) &
                            (obs_log["OBSTYPE"] == "SCIENCE") &
                            (obs_log["EXPTIME"] == flat_exptime)]

            print(f"Using {len(flats)} images to make a flat.")
            print(f"Applying dark correction to {len(flats)} images.")

            norm_flat_array = []
            for flat in flats:
                masterdark_name = os.path.join(masterdark_dir,
                                               f"masterdark_{flat['EXPTIME']}_"
                                               f"{board_id}_{sub_id}.fits")
                if not os.path.exists(masterdark_name):
                    continue
                dark_data = fits.getdata(masterdark_name)
                data = fits.getdata(flat['PATHNAME'])
                header = fits.getheader(flat['PATHNAME'])
                dark_sub_data = subtract_dark(data, dark_data)

                if save_dark_sub:
                    if not os.path.exists(dark_sub_dir):
                        os.makedirs(dark_sub_dir)
                    dark_subname = os.path.join(dark_sub_dir,
                                                os.path.basename(flat['PATHNAME']))
                    fits.writeto(dark_subname, data=dark_sub_data, header=header,
                                 overwrite=True)

                norm_flat_array.append(dark_sub_data / np.nanmedian(dark_sub_data))

            masterflat_data = np.nanmedian(norm_flat_array, axis=0)
            masterflat_header = header
            masterflat_header['EXPTIME'] = 1.0

            masterflatname = os.path.join(masterdark_dir,
                                          f"masterflat_{filter}_{board_id}_{sub_id}"
                                          f".fits")
            fits.writeto(masterflatname, masterflat_data, masterflat_header,
                         overwrite=True)

        field_dir_base_name = f"field_{fieldid}_{filter}_{board_id}_{sub_id}"

        # Detrend : dark subtract, flat field images and subtract background.
        if "detrend" in steps_to_perform:
            science = obs_log[(obs_log['FIELDID'] == fieldid) &
                              (obs_log['FILTERID'] == "J") &
                              (obs_log["BOARD_ID"] == board_id) &
                              (obs_log["EXPTIME"] == exptime)]

            print(f"Found {len(science)} images to reduce")

            detrend_namelist = []
            for sci in science:
                masterdark_name = os.path.join(masterdark_dir,
                                               f"masterdark_{sci['EXPTIME']}_{board_id}_"
                                               f"{sub_id}.fits")
                if not os.path.exists(masterdark_name):
                    print(f"No dark found for image{sci['BASENAME']} with "
                          f"exptime={sci['EXPTIME']}")
                    continue
                dark_data = fits.getdata(masterdark_name)
                data = fits.getdata(sci['PATHNAME'])
                header = fits.getheader(sci['PATHNAME'])
                dark_sub_data = subtract_dark(data, dark_data)

                masterflat_name = os.path.join(masterdark_dir,
                                               f"masterflat_{filter}_{board_id}_"
                                               f"{sub_id}.fits")
                flat_data = fits.getdata(masterflat_name)
                norm_flat_data = dark_sub_data / flat_data

                detrend_name = os.path.join(detrend_dir + '/' + field_dir_base_name,
                                            os.path.basename(sci['PATHNAME']))
                Path(detrend_name).parent.mkdir(parents=True, exist_ok=True)

                header["SATURATE"] = (40000 -
                                      np.nanmedian(dark_data)) / np.nanmedian(flat_data)

                fits.writeto(detrend_name, data=norm_flat_data, header=header,
                             overwrite=True)
                detrend_namelist.append(detrend_name)

            # Subtract background
            if bkg_sub_type == 'sextractor':
                detrend_mask_list = write_mask(detrend_namelist)
                skysub_list = []
                for ind, imgname in enumerate(detrend_namelist):
                    run_sextractor(imgname, weightimg=detrend_mask_list[ind])
                    skysub_name = imgname.replace('.fits', '') + '.bkgsub.fits'
                    skysub_list.append(skysub_name)
                    mask_sextractor_skysub_pixels(imgname, skysub_name)
            else:
                skysub_list = subtract_sky(detrend_namelist)

            copy_files_to_directory(skysub_list,
                                    astrometry_dir + '/' + field_dir_base_name)

        # Astrometry
        if "astrometry" in steps_to_perform:
            skysub_list = get_imagelist_from_directory(astrometry_dir + '/' +
                                                       field_dir_base_name,
                                                       select_type="*.bkgsub.fits")
            print(f"Found {len(skysub_list)} images to solve in "
                  f"{astrometry_dir + '/' + field_dir_base_name}")
            skysub_mask_list = write_mask(skysub_list)

            if astrometry_type == "astrometry_net":
                solved_images, failed_images \
                    = run_astrometry_net(skysub_list,
                                         masklist=skysub_mask_list,
                                         output_dir=astrometry_dir +
                                                    f'/{field_dir_base_name}',
                                         use_sextractor=True,
                                         sextractor_config_path=astrom_anet_sex,
                                         )
                print(f"Solved {len(solved_images)} images, "
                      f"failed to solve {len(failed_images)} images")

                run_scamp(solved_images, astrometry_dir + f'/{field_dir_base_name}')

            elif astrometry_type == "scamp":
                for ind, imgname in enumerate(skysub_list):
                    run_sextractor(imgname=imgname, weightimg=skysub_mask_list[ind],
                                   regions=True)
                    catname = imgname + ".cat"
                    regname = imgname + ".reg"

                    h = fits.getheader(imgname)
                    ra, dec = h['CRVAL1'], h['CRVAL2']
                    make_gaia_catalog(ra=ra, dec=dec,
                                      tmcatname=imgname.replace(".fits", ".gaia"),
                                      catalog_box_size_arcmin=15,
                                      catalog_min_mag=7,
                                      catalog_max_mag=20)

                    print('Downloaded gaia sources')
                    run_scamp_command(imgname)
                    run_swarp([imgname], [skysub_mask_list[ind]],
                              swarp_config_path=swarp_config,
                              out_path=imgname.replace(".fits", "_scamp.fits"))

        # Stack
        if "stack" in steps_to_perform:
            if not os.path.exists(stack_dir):
                os.makedirs(stack_dir)
            solved_images = get_imagelist_from_directory(astrometry_dir +
                                                         f'/{field_dir_base_name}',
                                                         select_type=
                                                         "*scampastrom.fits")
            solved_mask_paths = write_mask(solved_images)
            run_swarp(solved_images,
                      solved_mask_paths,
                      swarp_config_path=swarp_config,
                      out_path=f'{stack_dir}/stack.fits')
