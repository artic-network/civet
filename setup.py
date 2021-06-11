from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from civet import __version__, _program


setup(name='civet',
      version=__version__,
      packages=find_packages(),
      scripts=[
            "civet/scripts/civet.smk",
            "civet/scripts/build_catchment_trees.smk"
            ],
      package_data={"civet":["data/*"]},
      install_requires=[
            "biopython>=1.70",
            "baltic",
            "mako>=1.1"
        ],
      description='Cluster Investivation & Virus Epidemiology Tool',
      url='https://github.com/artic-network/civet',
      author='Aine OToole, Verity Hill, Rachel Colquhoun',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = civet.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
