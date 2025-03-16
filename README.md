# winterdrp_offline
Some scripts to reduce winter data.

# Install
```
cd winterdrp_offline
pip install -e .
```

# Usage 
See [notebooks/example.ipynb](winterdrp_offline/notebooks/example.ipynb) for an example of usage.

You will nees the following packages - 
Required : 
1. sextractor, swarp : Installable from conda
Optional (but required for running mirar): 
2. scamp, astrometry.net : Also installable with conda, but you will need to download the 
index files separately, from https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/
and put them in a folder named `data` inside the astrometry-net installation folder.

- sextractor: conda install -c conda-forge astromatic-source-extractor
- scamp: conda install -c conda-forge astromatic-scamp
- swarp: conda install -c conda-forge astromatic-swarp
