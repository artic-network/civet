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
            "civet/scripts/build_catchment_trees.smk",
            "civet/scripts/snipit_runner.smk",
            "civet/scripts/scorpio_runner.smk",
            "civet/scripts/generate_background_data.smk",
            "civet/scripts/global_snipit.smk"
            ],
      package_data={"civet":["data/*","data/report_modules/*","data/map_data/*","tests/action_test_data/*.fa"]},
      install_requires=[
            "biopython>=1.70",
            "mako>=1.1",
            "tabulate==0.8.9",
            "snipit"
        ],
      description='Cluster Investivation & Virus Epidemiology Tool',
      url='https://github.com/artic-network/civet',
      author='Aine OToole, Verity Hill & Rambaut Group',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = civet.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
