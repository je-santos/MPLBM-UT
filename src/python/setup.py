from setuptools import setup, find_packages

VERSION = '1.0'
DESCRIPTION = 'Utilities for the MPLBM-UT lattice-Boltzmann simulator.'
LONG_DESCRIPTION = 'The mplbm_utils package contains all the pre- and post-processing utilities to run MPLBM-UT.'

setup(
      name="mplbm_utils",
    version=VERSION,
    author="Alex Gigliotti and Javier Santos",
    author_email="alex.gigliotti@utexas.edu",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['porespy', 'numpy', 'pyyaml', 'vedo', 'pyvista', 'matplotlib', 'scikit-image', 'pandas'],
    python_requires='>=3.6',
    keywords=['python', 'lattice boltzmann method', 'lattice boltzmann', 'digital rocks'],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ]
)