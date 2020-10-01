from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources
# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from civet import __version__, _program

class NPMInstall(build_py):
    def run(self):
        self.run_command('npm install -g markdown-pdf --ignore-scripts')
        build_py.run(self)

setup(name='civet',
      version=__version__,
      packages=find_packages(),
      scripts=["civet/scripts/Snakefile",
      "civet/scripts/assess_input_file.smk",
      "civet/scripts/find_closest_cog.smk",
      "civet/scripts/cluster_civet.smk",
      "civet/scripts/process_collapsed_trees.smk",
      "civet/scripts/just_collapse_trees.smk",
      "civet/scripts/find_snps.smk",
      "civet/scripts/make_report.py",
      "civet/scripts/local_scale_analysis.py",
      "civet/scripts/check_cog_db.py",
      "civet/scripts/civet_template.pmd",
      "civet/scripts/COG_template.pmd",
      "civet/scripts/parse_closest.py",
      "civet/scripts/mapping.py",
      "civet/scripts/civetfunks.py"],
      package_data={"civet":["data/*","data/headers/*","data/mapping_files/*","data/vega_templates/*"]},
      install_requires=[
            "biopython>=1.70",
            "dendropy>=4.4.0",
            "pytools>=2020.1",
            'pandas>=1.0.1',
            'pysam>=0.15.4',
            "matplotlib>=3.2.1",
            "pweave>=0.30.3",
            "scipy>=1.4.1",
            "numpy>=1.13.3",
            "geopandas>=0.7.0",
            "snipit"
        ],
      cmdclass={
        'npm_install': NPMInstall
      },
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
