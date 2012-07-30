#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from distutils.command.build import build as _build

try:
    from sphinx.setup_command import BuildDoc as _BuildDoc
    has_sphinx = True
except ImportError:
    has_sphinx = False


class sdist(_sdist):
    """
    Modified sdist which first builds the sphinx documentation so we can include it.
    """
    def run(self):
        if has_sphinx: self.run_command("build_sphinx") 
        _sdist.run(self)


if has_sphinx:
    class BuildDoc(_BuildDoc):
        """
        Modified build_sphinx to build the documenation in the same 
        place "make html" would, under doc/_build.
        """
        def finalize_options(self):
            self.build_dir = "doc/_build"
            _BuildDoc.finalize_options(self)


class build(_build):
    """
    Modified build to run cosmoslik.py --build, then copy the built files 
    (which get built in-place in cosmoslik/) over to build/ so that install 
    automatically sees them and installs them. 
    """
    def run(self): 
        if has_sphinx: self.run_command("build_sphinx") 
        _build.run(self)
        #run cosmoslik --build
        #copy files to build/ directory for optional installation


cmdclass = {}
cmdclass['sdist']=sdist
cmdclass['build']=build
if has_sphinx: cmdclass['build_sphinx']=BuildDoc


setup(
    name='cosmoslik',
    version='0.1.0',
    author='Marius Millea',
    author_email='mmillea@ucdavis.edu',
    packages=find_packages(),
    url='http://pypi.python.org/pypi/cosmoslik/',
    license='LICENSE.txt',
    description='A modular cosmology likelihood sampler.',
    long_description=open('README.rst').read(),
    cmdclass=cmdclass
)
