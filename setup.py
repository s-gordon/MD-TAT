#!/usr/bin/env python
# Author: Shane Gordon
# Date Created:

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()


config = {
    'description': 'MDTAT: Molecular Dynamic Trajectory Analysis Tool',
    'author': 'Shane Gordon',
    'url': 'github.com/s-gordon/mdtat.git',
    'download_url': 'Where to download it.',
    'author_email': 'se2gordon@students.latrobe.edu.au',
    'version': '0.1',
    'install_requires': required,
    'packages': ['mdtat', 'mdtat.analysis'],
    'scripts': ['mdtat/scripts/analyse.py', 'mdtat/scripts/compress.py'],
    'name': 'mdtat'
}

setup(**config)
