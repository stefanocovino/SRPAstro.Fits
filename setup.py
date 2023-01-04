
from setuptools import setup, find_packages
from os import path
import glob
import SRPFITS as FITS
import os

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


# Look for script files
lscr = glob.glob(os.path.join('Scripts', 'SRP*'))
lscrex = []
for i in lscr:
    if os.path.splitext(i)[1] == '':
        lscrex.append(i)


# Look for data files
lsdt = glob.glob(os.path.join('SRPFITS/Data/*', '*'))
lsdtex = []
for i in lsdt:
    if os.path.splitext(i)[1] == '':
        lsdtex.append(i)


setup(
    name='SRPAstro.FITS',
    version=FITS.__version__,
    description='Tools for handling FITS files under SRP',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://pypi.python.org/pypi/SRPAstro.FITS',
    author='Stefano Covino',
    author_email='stefano.covino@inaf.it',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        ],
    keywords='astronomy data analysis',
    include_package_data=True,
    packages=find_packages(),
    python_requires='>=3',
    package_data={'SRPFITS':lsdtex},
    scripts=lscrex,
    install_requires=['SRPAstro >= 4.4', 'sep', 'photutils', 'astropy', 'astLib >= 0.11.5', 'astroalign'],
    ) 

