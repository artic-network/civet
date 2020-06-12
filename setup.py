from setuptools import setup, find_packages
import glob
import os
import pkg_resources
# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from civet import __version__, _program

setup(name='civet',
      version=__version__,
      packages=find_packages(),
      scripts=["civet/scripts/Snakefile",
      "civet/scripts/parse_paf.py",
      "civet/scripts/find_closest_cog.smk",
      "civet/scripts/assess_input_file.smk",
      "civet/scripts/process_catchment_trees.smk",
      "civet/scripts/make_report.py",
      "civet/scripts/data_parsing.py",
      "civet/scripts/utils/baltic.py",
      "civet/scripts/civet_template.pmd",
      "civet/scripts/make_tree_figures.py"],
      package_data={"civet":["data/HelveticaNeue.ttf","data/reference.fasta","data/polytomies.png"]},
      install_requires=[
            "biopython>=1.70",
            "dendropy>=4.4.0",
            "pytools>=2020.1",
            'pandas>=1.0.1',
            'pysam>=0.15.4',
            "matplotlib>=3.2.1",
            "pweave>=0.30.3",
            "scipy>=1.4.1",
            "numpy>=1.13.3"
        ],
      description='Cluster Investivation & Virus Epidemiology Tool',
      url='https://github.com/aineniamh/civet',
      author='Aine OToole, Verity Hill',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = civet.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
