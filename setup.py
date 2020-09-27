"""setup script for SV_parser relying on setup tools
"""
import setuptools


VERSION = "0.0.18"
PACKAGES = setuptools.find_packages( include = ["sv_parser", "sv_parser.*"] )
DEPENDENCIES = ["pandas", "PyVCF"]

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SV_parser",
    version=VERSION,
    author="Itsuki Sugita",
    author_email="itsukiisogin@gmail.com",
    description="Useful SV vcf parser",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dermasugita/SV_parser",
    packages=PACKAGES,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=DEPENDENCIES,
    keywords="bioinformatics",
)
