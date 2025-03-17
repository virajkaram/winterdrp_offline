import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="winterdrp_offline",
    version="0.0.1",
    author="Viraj Karambelkar",
    install_requires=[
        "astropy",
        "numpy",
        "scipy",
        "matplotlib",
        "astroquery",
    ],
    # Let setuptools find your package automatically from the current directory
    packages=setuptools.find_packages(),
    # Possibly also:
    python_requires=">=3.7",
)
