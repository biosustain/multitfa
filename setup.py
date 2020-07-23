""" Thermodynamics-based Flux Analysis

.. moduleauthor:: vishnu


"""

from setuptools import setup, find_packages

version_tag = "0.1.0"

setup(
    name="multiTFA",
    version=version_tag,
    author="vishnu",
    author_email="v.mahamkali@uq.edu.au",
    url="https://github.com/biosustain/tmfa",
    download_url="https://github.com/EPFL-LCSB/pytfa/archive/"
    + version_tag
    + ".tar.gz",
    install_requires=["cobra>0.13", "optlang", "pytest", "scipy", "numpy"],
    packages=find_packages(),
    python_requires=">=3.6, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4",
    description="multiTFA, Thermodynamics-based Flux Analysis in Python",
    keywords=["multiTFA", "tfa", "thermodynamics", "flux analysis"],
    license="Apache 2.0",
    include_package_data=True,
    package_data={"multiTFA": ["Data/*.pickle", "Data/*.map", "Data/*.npz"]},
    # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 4 - Beta",
        # Indicate who your project is intended for
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Environment :: Console",
        # Pick your license as you wish (should match "license" above)
        "License :: OSI Approved :: Apache Software License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
)

