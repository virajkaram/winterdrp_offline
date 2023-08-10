# winterdrp_offline
Some scripts to reduce winter data.

# Usage 
1. First, grab some unpacked and split images from the winter machine at caltech - `/data/loki/data/winter/<night>/raw_unpacked`. <br>
2. Then, grab an observing log from the winter machine - `/data/loki/observing_logs/observation_log_<night>.csv`. <br>
3. In `winterdrp_offline/run.py`, edit the USER-SPECIFIABLE options at the top of the file. <br>
4. Run using `python winterdrp_offline/run.py`.

You will nees the following packages - 
1. sextractor, scamp, swarp : Installable from conda
2. astrometry.net : Also installable with conda, but you will need to download the 
index files separately, from https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/
and put them in a folder named `data` inside the astrometry-net installation folder.

- utils git: https://github.com/winter-telescope/winter_utils/tree/main
    In winter_utils directory, pip install -e .
- sextractor: conda install -c conda-forge astromatic-source-extractor
- scamp: conda install -c conda-forge astromatic-scamp
- swarp: conda install -c conda-forge astromatic-swarp
