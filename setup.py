#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from distutils.command.build import build as _build
from distutils import file_util
import os


setup(
    name='cosmoslik-uspype',
    version='0.1.0',
    author='Marius Millea',
    author_email='mmillea@ucdavis.edu',
    packages=find_packages(),
#    namespace_packages = ['cosmoslik','cosmoslik.plugins'],
    description='USPype plugins for Cosmoslik.',
)
