""" Set-up script to install PyFAST locally
"""
from setuptools import setup

setup(name='pyfast',
      version='0.1',
      description='Tools for working with wind turbine simulator FAST',
      url='https://github.com/jennirinker/PyFAST.git',
      author='Jenni Rinker',
      author_email='jennifer.m.rinker@gmail.com',
      license='GPL',
      packages=['pyfast'],
      zip_safe=False)