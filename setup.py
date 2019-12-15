#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='wes-pipeline',
      version='1.0',
      description='Whole Exome Analysis Utilities',
      author='Italo Faria do Valle',
      author_email='italo.fariadovalle2@unibo.it',
      url='https://bitbucket.org/BBDA-UNIBO/wes-pipeline',
      package=['wes'],
      install_requires= ['pyvcf','pyyaml','pandas'],
      dependency_links=['https://github.com/jamescasbon/PyVCF/',
                        'https://bitbucket.org/xi/pyyaml'] 
     )

