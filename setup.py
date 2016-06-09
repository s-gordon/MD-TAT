#!/usr/bin/env python
# Author: Shane Gordon

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()


config = {
    'description': 'MD-TAT: Molecular Dynamics Trajectory Analysis Tool',
    'author': 'Shane Gordon',
    'url': 'https://github.com/s-gordon/md-tat',
    'download_url': 'https://github.com/s-gordon/md-tat',
    'author_email': 'se2gordon@students.latrobe.edu.au',
    'maintainer': 'Shane Gordon',
    'version': '0.1',
    'install_requires': required,
    'packages': ['mdtat', 'mdtat.analysis', 'mdtat.utils'],
    'scripts': ['mdtat/scripts/analyse.py', 'mdtat/scripts/compress.py'],
    'name': 'mdtat'
}

setup(**config)
