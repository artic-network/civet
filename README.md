# civet


## Quick links

  * [Requirements](#requirements)
  * [Install civet](#install-civet)
  * [Check the install worked](#check-the-install-worked)
  * [Updating civet](#updating-civet)
  * [Usage](#usage)
  * [Analysis pipeline](#analysis-pipeline)
  * [Output](#output)
  * [Source data](#source-data)
  * [Authors](#authors)
  * [Acknowledgements](#acknowledgements)
  * [References](#references)


### Requirements

civet runs on MacOS and Linux. The conda environment recipe may not build on Windows and is not supported but can be run using the Windows subsystem for Linux.

1. Some version of conda, we use Miniconda3. Can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html)
2. Access to CLIMB
3. Your input csv with minimally a column with `name` as a header
4. Optional fasta file with sequences you know will not yet be on CLIMB

### Install civet

1. Clone this repository and ``cd civet``
2. ``conda env create -f environment.yml``
3. ``conda activate civet``
4. ``python setup.py install`` or ``pip install .``

> Note: we recommend using pangolin in the conda environment specified in the ``environment.yml`` file as per the instructions above. If you can't use conda for some reason, dependency details can be found in the ``environment.yml`` file.


### Check the install worked

Type (in the civet environment):

```
civet -v
```
and you should see the version of civet printed.


### Updating civet

> Note: Even if you have previously installed ``civet``, as it is being worked on intensively, we recommend you check for updates before running.

To update:

1. ``conda activate civet``
2. ``git pull`` \
pulls the latest changes from github
3. ``python setup.py install`` \
re-installs civet
4. ``conda env update -f environment.yml`` \
updates the conda environment 


### Usage

1. Activate the environment ``conda activate civet``
2. Run ``civet <query> [options]``

Example usage:
> ``civet civet/tests/test.csv --fasta civet/tests/test.fasta --datadir path/to/data --remote -uun <your-user-name>``, where `<your-user-name>` represents your unique climb identifier.

Full usage:
```
positional arguments:
  query                 Input csv file with minimally "name" as a column header.
                        Can include additional fields to be incorporated into the
                        analysis, e.g. sample_date

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA         Optional fasta query.
  --remote              Remotely access lineage trees from CLIMB (note read
                        access required)
  -uun UUN, --your-user-name UUN
                        Your CLIMB COG-UK username. Required if running with
                        --remote flag
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: current working directory
  --datadir DATADIR     Data directory. Default with --remote flag will
                        rsync COG_UK data from CLIMB.
  -n, --dry-run         Go through the motions but don't actually run
  -f, --force           Overwrite all output
  -t THREADS, --threads THREADS
                        Number of threads
  --verbose             Print lots of stuff to screen
  --max-ambig MAXAMBIG  Maximum proportion of Ns allowed to attempt analysis.
                        Default: 0.5
  --min-length MINLEN   Minimum query length allowed to attempt analysis.
                        Default: 10000
  -v, --version         show program's version number and exit
```

### Analysis pipeline

- From the input csv (`<query>`), `civet` attempts to match the ids with COG-UK ids in the up-to-date metadata database.

- If the id matches with a record in COG-UK, the corresponding metadata is pulled out.

- If the id doesn't match with a record in COG-UK and a fasta sequence of that id has been provided, it's passed into a workflow to identify the closest sequence in COG-UK. In brief, this search consists of quality control steps that maps the sequence against a reference (`MN908947.3`), pads any indels relative to the reference and masks non-coding regions. civet then runs a `minimap2` search against the COG-UK database and finds the best hit to the query sequence.

- The metadata for the closest sequences are also pulled out of the large COG-UK database.

- For each set of sequences in the large COG-UK tree, local subtrees are pruned out and vizualised.

- A report summarising the query sequences is generated, providing information about global and UK lineages.

### Output

Your output will be a markdown report, with summaries of lineages and genetic diversity present in your query. Trees are visualised in the report and compared to the diversity of lineages present in the community.

### Source data

Data associated with COG-UK is pulled from CLIMB, which is why access to CLIMB is required to run remotely. 

Data files:

- `/cephfs/covid/bham/artifacts/published/latest/phylogenetics/alignments/cog_2020-06-05_metadata.csv`

- `/cephfs/covid/bham/artifacts/published/latest/phylogenetics/trees/cog_global_2020-06-05_tree.newick`

- `/cephfs/covid/bham/artifacts/published/latest/phylogenetics/alignments/cog_2020-06-05_all_alignment.fasta`


### Authors

`civet` was created by Áine O'Toole, Verity Hill & Stefan Rooke, on behalf of the COG_UK consortium.

### Acknowledgements

We acknowledge the hard work from all members of the COG-UK consortium that has gone into generating the data used by `civet`.

`civet` makes use of [`datafunk`](https://github.com/cov-ert/datafunk) and [`clusterfunk`](https://github.com/cov-ert/clusterfunk) functions which have been written by members of the Rambaut Lab, specificially Rachel Colquhoun, JT McCrone, Ben Jackson and Shawn Yu.

[`baltic`](https://github.com/evogytis/baltic/tree/master/baltic) by Gytis Dudas is used to visualize the trees.

### References

[`minimap2`](https://github.com/lh3/minimap2) 

Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191
