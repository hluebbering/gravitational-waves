# Initializes project after it has been cloned + describes all of the metadata about project
# Required fields: name, version, and packages

from distutils.core import setup

setup(
  name = 'GravitationalWaves', # Must be unique if published to PyPI
  version = '0.1dev', # Keeps track of different project releases
  packages = ['towelstuff',], # Describes where the Python source code is within project
  license = 'Creative Commons Attribution-Noncommercial-Share Alike license', # Includes information about the license 
  long_description = open('README.md').read(), # Re-uses README file for long_description field
)
