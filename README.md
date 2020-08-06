# civet

**C**luster **I**nvestigation & **V**irus **E**pidemiology **T**ool

<img src="./docs/civet_logo.png" width="450">

civet is a tool developed with 'real-time' genomics in mind. With the large phylogeny available through the COG-UK infrastructure on CLIMB, civet will generate a report for a set of sequences of interest i.e. an outbreak investigation. If the sequences are already on CLIMB and part of the large tree, civet will pull out the local context of those sequences, merging the smaller local trees as appropriate. If sequences haven't yet been uploaded to CLIMB, for instance if they have just been sequenced, civet will find the closest sequence in the COG-UK database on climb, pull the local tree of that sequence out and add your sequence in. The local trees then get collapsed to display in detail only the sequences of interest so as not to inform investigations beyond what was suggested by epidemiological data. 

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
  * [Software versions](#software-versions)


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

> Note: we recommend using civet in the conda environment specified in the ``environment.yml`` file as per the instructions above. If you can't use conda for some reason, dependency details can be found in the ``environment.yml`` file.


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
> ``civet civet/tests/test.csv --fasta civet/tests/test.fasta --remote -uun <your-user-name>``, where `<your-user-name>` represents your unique CLIMB identifier.

Full usage:
```
usage: civet <query> [options]

civet: Cluster Investivation & Virus Epidemiology Tool

positional arguments:
  query                 Input csv file with minimally `name` as a column
                        header. Can include additional fields to be
                        incorporated into the analysis, e.g. `sample_date`

optional arguments:
  -h, --help            show this help message and exit
  -i, --id-string       Indicates the input is a comma-separated id string
                        with one or more query ids. Example:
                        `EDB3588,EDB3589`.
  --fasta FASTA         Optional fasta query.
  -sc SEQUENCING_CENTRE, --sequencing-centre SEQUENCING_CENTRE
                        Customise report with logos from sequencing centre.
  --CLIMB               Indicates you're running CIVET from within CLIMB, uses
                        default paths in CLIMB to access data
  --cog-report          Run summary cog report. Default: outbreak
                        investigation
  -r, --remote-sync     Remotely access lineage trees from CLIMB, need to also
                        supply -uun,--your-user-name
  -uun UUN, --your-user-name UUN
                        Your CLIMB COG-UK username. Required if running with
                        --remote-sync flag
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: current working directory
  -b, --launch-browser  Optionally launch md viewer in the browser using grip
  --datadir DATADIR     Local directory that contains the data files
  --fields FIELDS       Comma separated string of fields to colour by in the
                        report. Default: country
  --label-fields LABEL_FIELDS
                        Comma separated string of fields to add to tree report
                        labels.
  --search-field SEARCH_FIELD
                        Option to search COG database for a different id type.
                        Default: COG-UK ID
  --distance DISTANCE   Extraction from large tree radius. Default: 2
  -g, --global          Search globally.
  -n, --dry-run         Go through the motions but don't actually run
  --tempdir TEMPDIR     Specify where you want the temp stuff to go. Default:
                        $TMPDIR
  --no-temp             Output all intermediate files, for dev purposes.
  --collapse-threshold THRESHOLD
                        Minimum number of nodes to collapse on. Default: 1
  -t THREADS, --threads THREADS
                        Number of threads
  --verbose             Print lots of stuff to screen
  --max-ambig MAXAMBIG  Maximum proportion of Ns allowed to attempt analysis.
                        Default: 0.5
  --min-length MINLEN   Minimum query length allowed to attempt analysis.
                        Default: 10000
  --local-lineages      Contextualise the cluster lineages at local regional
                        scale. Requires at least one adm2 value in query csv.
  --date-restriction    Chose whether to date-restrict comparative sequences
                        at regional-scale.
  --date-range-start DATE_RANGE_START
                        Define the start date from which sequences will COG
                        sequences will be used for local context. YYYY-MM-DD
                        format required.
  --date-range-end DATE_RANGE_END
                        Define the end date from which sequences will COG
                        sequences will be used for local context. YYYY-MM-DD
                        format required.
  --date-window DATE_WINDOW
                        Define the window +- either side of cluster sample
                        collection date-range. Default is 7 days.
  --map-sequences       Maps coordinates of sequences. Requires arguments x-col,
                        y-col and input-crs.
  --x-col X_COL         Name of column in input csv containing x coordinates of sequences for mapping.
  --y-col Y_COL         Name of column in input csv containing x coordinates of sequences for mapping.
  --input-crs           Coordinate system of sequence coordinates in EPSG numbers eg World Mercator (WGS84) is EPSG:4326, 
                        and ordnance survey is EPSG:27700.
                        For more information see https://geopandas.org/projections.html and https://spatialreference.org/ref/epsg/
  --mapping-trait MAPPING_TRAIT
                        Trait to colour by on the map. Must match a column header in the input csv.
  -v, --version         Show program's version number and exit

```

### Analysis pipeline

Detailed description of civet found [here](./docs/civet_detailed_description.md).

Overview:

<img src="./docs/workflow_diagram.png" width="700">

- From the input csv (`<query>`), `civet` attempts to match the ids with COG-UK ids in the up-to-date metadata database.

- If the id matches with a record in COG-UK, the corresponding metadata is pulled out.

- If the id doesn't match with a record in COG-UK and a fasta sequence of that id has been provided, it's passed into a workflow to identify the closest sequence in COG-UK. In brief, this search consists of quality control steps that maps the sequence against a reference (`MN908947.3`), pads any indels relative to the reference and masks non-coding regions. civet then runs a `minimap2` search against the COG-UK database and finds the best hit to the query sequence.

- The metadata for the closest sequences are also pulled out of the large COG-UK database.

- Combining the metadata from the COG-UK records of the closest hit and the exact matching records found in COG-UK, `civet` queries the large global phylogeny (also from COG-UK database)containing all COG-UK and all GISAID sequences. The local trees around the relevant tips are pruned out of the large phylogeny, merging overlapping local phylogenys as needed.

- If these local trees contain "closest-matching" tips, the sequence records for the tips on the tree and the sequences of the relevant queries are added into an alignment. Any peripheral sequences coming off of a polytomy are collapsed to a single node and summaries of the tip's contents are output.

- After collapsing the nodes, civet runs `iqtree` on the new alignment, now with query sequences in. Optionally, the `--delay-tree-collapse` argument will wait to collapse nodes until after `iqtree` has added the new query sequences in, but be wary as some of these local trees can be very large and may take a number of hours to run. 

- `civet` then generates a report summarising the query sequences, providing information about global and UK lineages.

### Output

Your output will be a markdown report, with summaries of lineages and genetic diversity present in your query. Trees are visualised in the report and compared to the diversity of lineages present in the community.

**An example report can be found [here](https://github.com/COG-UK/civet/blob/master/docs/civet_report_example.md).**

The default output is a markdown file, which can then be converted to a file format of your choice. In addition to this, if you provide the `--launch-browser` option in the command line, an html document will be outputted using grip (https://github.com/joeyespo/grip), which will also appear in your browser. You can then save this as a pdf using your browser.

### Source data

Data associated with COG-UK is pulled from CLIMB, which is why access to CLIMB is required to run remotely. 

Data files:

- `/cephfs/covid/bham/civet-cat/cog_global_tree.nexus`

- `/cephfs/covid/bham/civet-cat/cog_metadata.csv`

- `/cephfs/covid/bham/civet-cat/cog_metadata_all.csv`

- `/cephfs/covid/bham/civet-cat/cog_global_alignment.fasta`

- `/cephfs/covid/bham/civet-cat/cog_alignment.fasta`

- `/cephfs/covid/bham/civet-cat/cog_alignment_all.fasta`

### Authors

`civet` was created by Áine O'Toole & Verity Hill & Rambaut Group on behalf of the COG_UK consortium.

### Acknowledgements

We acknowledge the hard work from all members of the COG-UK consortium that has gone into generating the data used by `civet`.

`civet` makes use of [`datafunk`](https://github.com/cov-ert/datafunk) and [`clusterfunk`](https://github.com/cov-ert/clusterfunk) functions which have been written by members of the Rambaut Lab, specificially Rachel Colquhoun, JT McCrone, Ben Jackson and Shawn Yu.

[`baltic`](https://github.com/evogytis/baltic/tree/master/baltic) by Gytis Dudas is used to visualize the trees.

### References

[`minimap2`](https://github.com/lh3/minimap2) 

Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191

[iqtree](http://www.iqtree.org/#download)

L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300

D.T. Hoang, O. Chernomor, A. von Haeseler, B.Q. Minh, L.S. Vinh (2018) UFBoot2: Improving the ultrafast bootstrap approximation. Mol. Biol. Evol., 35:518–522. https://doi.org/10.1093/molbev/msx281

Stéphane Guindon, Jean-François Dufayard, Vincent Lefort, Maria Anisimova, Wim Hordijk, Olivier Gascuel, New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0, Systematic Biology, Volume 59, Issue 3, May 2010, Pages 307–321, https://doi.org/10.1093/sysbio/syq010

[snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

### Software versions

    - python=3.6
    - snakemake-minimal=5.13 
    - iqtree=1.6.12
    - minimap2=2.17-r941
    - pandas==1.0.1
    - pytools=2020.1
    - dendropy=4.4.0
    - tabulate=0.8.7
    
