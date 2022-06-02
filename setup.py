from os import path
from setuptools import find_packages, setup

if path.exists('README.md'):
    with open('README.md') as readme:
        long_description = readme.read()



setup(
    name="highz_qso_arxiv",
    version='0.0.dev0',
    description="Repository for High-z Quasar Search Reductions",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/enigma-igm/highz_qso_arxiv",
    author='Daming Yang and enigma folks',
    packages=find_packages(exclude=["tests"]),
    install_requires=[
        "numpy",
        "astropy",
        "matplotlib",
    ], dependency_links=[],
    scripts=[],)
