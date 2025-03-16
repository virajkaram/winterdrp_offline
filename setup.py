import setuptools
from pathlib import Path

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='winterdrp_offline',
    version='0.0.1',
    author='Viraj Karambelkar',
    requires=['astropy', 'numpy', 'scipy', 'matplotlib',
              'astroquery', 'pathlib', 'os', 'glob', 'matplotlib',
              ],
)

