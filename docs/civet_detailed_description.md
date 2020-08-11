## Detailed description of civet

### 1) Finding appropriate local files and output directories

- Find master Snakefile
- If specified, check if input fasta file exists, otherwise exit with informative message 
- Output directory:
	if given, get path and create directory if necessary
	if not given, create a timestamped output directory
- Set up temporary directory, can be specified. If `--no-temp`, set it to the output directory
- Check if input query is a file:
	if not and if the `-ids` flag is given, parse the input comma-separated id string
	if not and the -ids flag isn't given, exit with informative message
- Check input query file has appropriate headers (i.e. at least `name`)
- If no input `--fields`, default to adm1

- Initialise config dictionary to pass through to snakemake

### 2) Finding COG-UK data

- If `--CLIMB`, check path to civet directory exists, else exit
- If `--datadir`, if not `--remote`, check if it has the appropriate files, else exit
- If not `--datadir`, default to current working directory
- If `--remote` and `-uun` have been specified, rsync data from CLIMB
	If rsync fails exit with error showing username
- If `-uun` hasn't been specified alongside `--remote` exit with error
- If neither `--remote` nor `--datadir` nor `--CLIMB` have been given, exit with error

### 3) QC on input fasta file
- Check minimum length (Default: 10,000)
- Check maximum ambiguity content (Default: 0.5)
- Write out failed query ids and why they failed (to pass to report)
- Write out sequences that passed QC to pass into pipeline

### 4) Access installed package data
- Find reference fasta file for downstream datafunk command
- Find polytomies figure for the report
- Find the report template
- If `--sequencing-centre` specified, check it's a valid centre name and find the appropriate header for the report. If not a real sequencing centre, exit and display valid names

### 5) Miscellaneous arguments
- Check if `--distance` an integer, else exit
- Add `--verbose` to config 
- --global boolean to config
- --delay-tree-collapse to config
- --dry-run to snakemake
- --force to snakemake and config

### 6) Call master snakemake
Passing through the configuration dictionary and setting the working directory to the temporary directory (if --no-temp, this inherits the output directory).  

### 7) Within master Snakefile
- Create sequencing centre flag string for report generator
- If `--global` specified, sequence database is all global sequences and COG-UK, if not specified the default is to just search within COG-UK
- Reformat miscellaneous arguments
- Call target rule all and include `assess_input_file.smk`

### 8) assess_input_file.smk
1) check_cog_db
Check against the COG-UK phylogenetics database for your query id, will check search field, which by default is `central_sample_id`, this finds all sequences that will already be in the big tree
2) check_cog_all
Check against the whole COG-UK database, including the sequences that failed the original phylogenetics QC, will find any record that has been uploaded to COG-UK and passes Sam's QC step. If your query was uploaded to CLIMB but didn't pass Sam's QC step, it will not be found and you should supply it separately, but also question whether your sequence is of good enough quality to give a reliable result.
3) get_closest_cog
If you have supplied an additional fasta file with sequences that haven't yet been uploaded to the COG-UK database, or if at least one sequence in your query hadn't passed the phylogenetics QC thresholds, this sub-pipeline will be spawned off. In short, it finds the closest matching COG-UK sequence, or closest --global sequence if the global flag is used, details in Section 9. This workflow returns a csv of the closest hits within the COG-UK database for each query not in the tree and a csv of sequences that were not present in the database and didn't have a sequence supplied. This information gets added to the report. It also returns the input query sequences aligned, with the UTRs masked out. 
4) combine_metadata
Metadata from the closest COG sequence is combined with the metadata found directly within the database.
5) prune_out_catchments
Using the names of the closest COG sequences and of the sequences already in the tree, local trees with a radius of `--distance` number of nodes are pulled out of the large global tree. 
6) process_catchments
Another sub-workflow that gets spawned off with the set of catchment trees produced from the prune_out_catchments rule. Full details of this step are found in Section 10, but in brief there are three alternative snakefiles can can be called at this step depending on whether --delay-tree-collapse has been called and depending on whether new sequences need to be added into the tree or not. The output of this step is a directory of "local_trees" that have the nodes we have not queried collapsed to the nearest polytomy, with an output summary file of the tips that had been in that sub-tree prior to collapse.
7) make_report

### 9) Making the report

Parses metadata CSVs and any input csv provided, along with tree files generated previously. There are a number of different options for figures, some of which must be specified when calling CIVET. 

NB There are two templates available for use, but for the vast majority of CIVET runs, the default option will be sufficient. The other option (call `--cog-report` ) is a slimmed down report designed to explore new sequences that were added to the database, rather than a cluster investigation.

Core report:

1) Custom header - use the `--sc` argument when calling CIVET, followed by one option from BIRM, CAMB, EDIN, EXET, GLAS, GSTT, LIVE, LOND, NORT, NORW, NOTT, OXON, PHEC, PHWC, PORT, SANG or SHEF. This will provide logos for the appropriate sequencing centre at the top of the report. If not specified, the default header will be inserted, which is simply the CIVET and COG logos.
2) A table is produced showing all the sequences that have been provided in the input csv. Any additional fields provided in the input csv and specified on the command line (see below) will be included in this table, along with the default columns of Query ID, sequence name in the tree, sample date, closest sequence in the tree (if a fasta is provided and the sequence is not in the COG database), UK lineage, Global lineage, Phylotype, and Tree number. 
At this point, any sequences which were in the input csv but not in the fasta (if provided) or the COG database will be flagged, as well as any which did not meet quality control requirements.
3) The trees are rendered. If they are too large (ie more than 1000 sequences), they will not render due to pixel restrictions in the python packages used. To colour by a specific trait, use the `--fields` argument and provide the column titles in the input csv of the desired trait. It is possible to colour by more than one trait by providing a comma separated string at the command line, and different coloured dots will be placed next to the tree tips. For clarity of visualisation, we recommend not colouring by more than a few traits so that there is space to clearly see the different dots.
If there are many factors in a trait, it will be difficult to differentiate between colours. In this case, traits can also be added to the tip labels of the trees. To do this, use the `--label-fields` argument, again specifying the relevant column titles in a comma separated string. 
The default colouring scheme is Adm1 (eg England, Wales, Scotland and Northern Ireland) if available, and the base tip label is "ID|County|sample_date".
4) Tree background - bar charts displaying more detailed country information for collapsed nodes. They only show the ten largest countries in each collapsed node.

Optional figures:

1) Putting sequences on a map using coordinates: if the coordinates of the sequences are present, they can be plotted on a map using `--map-sequences`. There are three compulsory command line options after this: `--x-col` and `--y-col` which must contain the column headers of the x and y coordinates in the spreadsheet; and `--input-crs` which is the coordinate reference system that these coordinates are in. For example, latitute and longitude (WGS84) is "EPSG:4326". For more information see https://geopandas.org/projections.html. Whatever the input is, it will be converted to a flattened projection. If the map looks very odd (eg sequence data points not being anywhere near the land), you have likely provided the wrong input CRS.
Sequences will be placed on the map of the UK which will then be filtered to only show the region of interest, along with an overlay of urban areas and labels of counties and major cities.
An optional additional argument `--mapping-trait` may be provided, which will colour the dots on the map by a specific trait, again sepcified by providing the  appropriate column header.

2) Local background diversity of UK lineages: Provides a map of the diversity of lineages in the adm2 that the sequences in the query dataset are present in, and in neighbouring regions by using `--local-lineages` flag. 
Options: `--date-restriction` flag allows the diversity to be examined for a date range, specified either by using `--date-range-start` and `--date-range-end`, both in YYYY-MM-DD format, or by using `--date-window` which finds the first and last sample date in the query dataset and adds the window onto either end. Default is 7 days.

### 10) find_closest_cog.smk
1) non_cog_minimap2_to_reference
Map the query fasta file containing sequences from CLIMB that didn't meet the phylogenetic cut offs and sequences from the input fasta file against the Wuhan-Hu-1 reference sequence (Genbank ID: MN908947.3). This mapping runs `minimap2` with the present option `asm5` that is configured for mapping a long, closely related assembly against a reference. 
2) non_cog_remove_insertions_and_trim_and_pad
This step runs the [`datafunk`](https://github.com/cov-ert/datafunk) command `sam_2_fasta`, which parses the cigar strings produced by `minimap2` and scaffolds the query sequence against the reference by removing insertions relative to the reference. The insertions are logged in a text file for reference, and can be retrieved by running `civet` with the `--no-temp` flag. This step also masks the UTR regions with N's as the terminal ends of the genome. 
3) minimap2_against_cog
The padded, mapped set of query sequences then get mapped against the entire COG-UK database, or if the `--global` flag is used against all COG-UK and all GISAID, and finds the top hit within the database for each query. Again, `minimap2` is run with the `asm5` presets. 
4) parse_paf
The output from minimap2 is parsed and the relevant COG-UK metadata and fasta sequences are pulled out of the big database.

### 11) process_catchments

#### just_collapse_trees.smk
Runs if no new sequences need to be added to the phylogeny.
1) summarise_polytomies
`clusterfunk focus` is run, which collapses and summarises non-query or non-closest sequences to the parent polytomy and outputs a summary file of the tips that make up the collapsed node.
2) remove_str_for_baltic
Processes the tip names to make them readable for the baltic tree parser. 
3) to_nexus
Converts each tree from newick to nexus format.
4) summarise_processing
Combines all the information about collapsed tips into a single file.

#### process_collapsed_trees.smk
Runs by default if new trees need to be estimated with query sequences added in.
1) summarise_polytomies
`clusterfunk focus` is run, which collapses and summarises non-query or non-closest sequences to the parent polytomy and outputs a summary file of the tips that make up the collapsed node. 
2) get_collapsed_representative
For each collapsed node, the sequences of the tips that comprise that collapsed node are pulled out of the global alignment and the sequence with the least amount of ambiguity is selected to represent that node in the phylogeny.
3) extract_taxa
All taxon labels are extracted from the collapsed phylogeny.
4) gather_fasta_seqs
From the the respective sequences from each tip label in the tree are extracted from either the global phylogeny,the 'get_collapsed_representative' step or from the aligned input query sequences.
5) hash_for_iqtree
A hash for taxon labels is created for the input to `iqtree` so it doesn't corrupt the sequence names.
6) hash_tax_labels
`clusterfunk relabel_tips` is run that relabels the tips of the tree with the taxon hash.
7) iqtree_catchment
`iqtree` is run on the new alignment with added sequences, with `-m HKY` and `-au` options.
8) restore_tip_names
The hashed tip names are restored with the same `clusterfunk relabel_tips` command.
9) remove_str_for_baltic
Processes the tip names to make them readable for the baltic tree parser.
10) to_nexus
Converts each tree from newick to nexus format.
11) summarise_processing
Combines all the information about collapsed tips into a single file.
