import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bigants",
    version="0.0.8",
    author="Olga Lazareva",
    author_email="olga.lazareva@tum.de",
    description="BiGAnts - a package for network-constrained biclustering of omics data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/biomedbigdata/BiGAnts-PyPI-package",
    packages=setuptools.find_packages(exclude=['script_main.py']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'pandas',
        'numpy',
        'networkx',
	'matplotlib',
	'scipy'],
	

)
