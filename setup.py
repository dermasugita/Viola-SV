"""setup script for SV_parser relying on setup tools
"""
import setuptools


VERSION = "0.1"
#PACKAGES = setuptools.find_packages( include = ["sv_parser", "sv_parser.*"] )
PACKAGES = setuptools.find_packages('src')
DEPENDENCIES = ["pandas", "PyVCF", "intervaltree", "biopython"]

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PySgt",
    version=VERSION,
    author="Itsuki Sugita",
    author_email="itsukiisogin@gmail.com",
    description="Useful SV vcf parser",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dermasugita/SV_parser",
    packages=PACKAGES,
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
    install_requires=DEPENDENCIES,
    keywords="bioinformatics",
)
