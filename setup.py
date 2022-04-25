"""setup script for Viola-SV relying on setup tools
"""
import setuptools


PACKAGES = setuptools.find_packages('src')
DEPENDENCIES = ["pandas >= 1.1.5",
                "PyVCF >= 0.6.8",
                "intervaltree >= 3.1.0",
                "biopython >= 1.79",
                "scikit-learn >= 0.24.2",
                "scipy >= 1.5.4",
                "click >= 8.0.4",]

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
    include_package_data=True,
    entry_points={
        "console_scripts":[
            "viola=viola.cli.viola:viola",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "programming language :: Python :: 3.6",
        "programming language :: Python :: 3.7",
        "programming language :: Python :: 3.8",
        "programming language :: Python :: 3.9",
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
