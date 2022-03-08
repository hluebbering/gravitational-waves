# Initializes project after it has been cloned + describes all of the metadata about project
# Required fields: name, version, and packages

from distutils.core import setup

setup(
    name = 'gravitational-waves',
    version = '0.0.0', # '0.1dev'
    packages = ['gravitational-waves',], # Describes where the Python source code is within project
    license = 'MIT',
    long_description = open('README.md').read(),
    url='https://github.com/hluebbering/gravitational-waves',
)
