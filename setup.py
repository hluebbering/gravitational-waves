# Initializes project after it has been cloned + describes all of the metadata about project
# Required fields: name, version, and packages

from distutils.core import setup

setup(
    name = 'gravitational-waves',
    version = '0.0.1', 
    packages = ['gravitational-waves',], # Describes where the Python source code is within project
    license = 'MIT',
    author = "Hannah Luebbering",
    author_email = "luebhr@gmail.com",
    description = "Gravitational wave analysis in Python",
    long_description = open('README.md').read(),
    url = 'https://github.com/hluebbering/gravitational-waves',
)
