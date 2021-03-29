"""setup script for SV_parser relying on setup tools
"""
import setuptools


PACKAGES = setuptools.find_packages('src')
DEPENDENCIES = ["pandas",
                "PyVCF",
                "intervaltree",
                "biopython",
                "scikit-learn",
                "click",]

with open("README.rst", "r") as fh:
    long_description = fh.read()

exec(open("src/viola/_version.py", "r").read())

setuptools.setup(
    name="Viola-SV",
    version=__version__,
    author="Itsuki Sugita, Shohei Matsuyama, Hiroki Dobashi",
    author_email="itsukiisogin@gmail.com",
    license='Apache License 2.0',
    description="SV signature analysis tool with custom SV type",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/dermasugita/Viola-SV",
    packages=PACKAGES,
    package_dir={'': 'src'},
    entry_points={
        "console_scripts":[
            "vcf2bedpe=viola.cli.viola:main",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    extras_require={
        'test': ['pytest']
    },
    install_requires=DEPENDENCIES,
    keywords="bioinformatics",
)
