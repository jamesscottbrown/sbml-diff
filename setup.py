#!/usr/bin/env python

from setuptools import setup

setup(name='sbml-diff',
      version='1.0',
      description='sbml-diff is a tool for visually comparing SBML models',
      author='James Scott-Brown',
      author_email='james@jamesscottbrown.com',
      url='',
      packages=['sbml_diff'],
      scripts=['sbml-diff.py'],
      install_requires=['BeautifulSoup']
      )
