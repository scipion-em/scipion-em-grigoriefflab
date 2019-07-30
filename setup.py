"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from grigoriefflab import PLUGIN_VERSION
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Load constants, specially PLUGIN_VERSION
exec(open(path.join(here, "spider", "constants.py")).read())

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-em-grigoriefflab',  # Required
    version=PLUGIN_VERSION,  # Required
    description='Grigoriefflab ready to use in scipion.',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-em/scipion-em-grigoriefflab',  # Optional
    author='I2PC',  # Optional
    author_email='scipion@cnb.csic.es',  # Optional
    keywords='scipion cryoem imageprocessing scipion-2.0',  # Optional
    packages=find_packages(),
    package_data={  # Optional
       'grigoriefflab': ['grigoriefflab_logo.png', 'protocols.conf'],
    }
)
