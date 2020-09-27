from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='AeroMDAO',
    version="0.0.1dev",
    description='A Python 3 package for investigations into aerodynamic shape optimization, aircraft design, multidisciplinary optimization and surrogate modeling.',
    packages=find_packages(),
    install_requires=['numpy', 'scipy', 'matplotlib'],
    python_requires='>=3.6',
)