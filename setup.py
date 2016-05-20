from setuptools import setup, find_packages

# Check if the user has `pyfits` as part
# of the `astropy` distribution...
try:
  import astropy.io.fits
  pyfits = 'astropy'
except:
  pyfits = 'pyfits'

# Get the long description from the README
def readme():
  with open('README.rst') as f:
    return f.read()

# Setup!
setup(name = 'everest',
      version = '1.0',
      description = 'EPIC Variability Extraction and Removal for Exoplanet Science Targets',
      long_description = readme(),
      classifiers = [
                      'Development Status :: 3 - Alpha',
                      'License :: OSI Approved :: MIT License',
                      'Programming Language :: Python :: 3.4',
                      'Topic :: Scientific/Engineering :: Astronomy',
                    ],
      url = 'http://github.com/rodluger/everest',
      author = 'Rodrigo Luger',
      author_email = 'rodluger@uw.edu',
      license = 'MIT',
      packages = find_packages(),
      install_requires = [
                          'numpy',
                          'scipy',
                          'matplotlib',
                          'george',
                          'sklearn',
                          'astroML',
                          'six',
                          pyfits,
                         ],
      include_package_data = True,
      zip_safe = False)