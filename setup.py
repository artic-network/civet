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
      scripts=[
                ],
      install_requires=[
            "biopython>=1.70",
            "dendropy>=4.4.0",
            "pytools>=2020.1",
            'pandas>=1.0.1',
            'pysam>=0.15.4'
        ],
      description='CoV Introductions & Variant Epi Tool',
      url='https://github.com/aineniamh/civet',
      author='Aine OToole, Verity Hill & Stefan Rooke',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = civet.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)