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
This step summarises the metadata and visualises the tree in a markdown report. 
