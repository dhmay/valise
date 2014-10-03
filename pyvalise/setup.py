#!/usr/bin/env python

from setuptools import setup

setup(name='pyvalise',
      version='0.1dev',
      description='Python code for Valise',
      author='Damon May',
      author_email='damonmay@uw.edu',
      packages=['pyvalise'],
      license='None',
      install_requires=['numpy', 'scipy', 'matplotlib', 'matplotlib_venn', 'pyteomics', 'lxml', 'Bio', 'statsmodels', 'pandas', 'patsy', 'rpy2', 'pysam', 'xlsxwriter'],
      long_description=open('README.md').read(),
     )
