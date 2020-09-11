import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SV_parser",
    version="0.0.1",
    author="Itsuki Sugita",
    author_email="itsukiisogin@gmail.com",
    description="Useful SV vcf parser",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dermasugita/SV_parser",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
