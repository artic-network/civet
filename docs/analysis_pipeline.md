

<section id="banner">
    <div class="content">
      <header>
        <h2>civet pipeline</h2>
        <p>Detailed description of civet analysis pipeline</p>
      </header>
    </div>
    <span class="image object">
        <img src="./figures/civet_logo.png" alt="" style="max-width:150px"/>
        </span>
</section>

<img src="./figures/workflow_diagram.png" width="700">


**An example report can be found [here](https://github.com/COG-UK/civet/blob/master/docs/civet_report_example.md).**

The default output is a markdown file, which can then be converted to a file format of your choice. In addition to this, if you provide the `--launch-browser` option in the command line, an html document will be outputted using grip (https://github.com/joeyespo/grip), which will also appear in your browser. You can then save this as a pdf using your browser.


- From the input csv (`<query>`), `civet` attempts to match the ids with COG-UK ids in the up-to-date metadata database.

- If the id matches with a record in COG-UK, the corresponding metadata is pulled out.

- If the id doesn't match with a record in COG-UK and a fasta sequence of that id has been provided, it's passed into a workflow to identify the closest sequence in COG-UK. In brief, this search consists of quality control steps that maps the sequence against a reference (`MN908947.3`), pads any indels relative to the reference and masks non-coding regions. civet then runs a `minimap2` search against the COG-UK database and finds the best hit to the query sequence.

- The metadata for the closest sequences are also pulled out of the large COG-UK database.

- Combining the metadata from the COG-UK records of the closest hit and the exact matching records found in COG-UK, `civet` queries the large global phylogeny (also from COG-UK database)containing all COG-UK and all GISAID sequences. The local trees around the relevant tips are pruned out of the large phylogeny, merging overlapping local phylogenys as needed.

- If these local trees contain "closest-matching" tips, the sequence records for the tips on the tree and the sequences of the relevant queries are added into an alignment. Any peripheral sequences coming off of a polytomy are collapsed to a single node and summaries of the tip's contents are output.

- After collapsing the nodes, civet runs `iqtree` on the new alignment, now with query sequences in. Optionally, the `--delay-tree-collapse` argument will wait to collapse nodes until after `iqtree` has added the new query sequences in, but be wary as some of these local trees can be very large and may take a number of hours to run. 

- `civet` then generates a report summarising the query sequences, providing information about global and UK lineages.


### 1) Initialising directories & detect input type

- Initialise config dictionary to pass through to snakemake, load defaults
- Detect the input type (``-i / --input``), either csv, yaml/yml or a string of IDS
- If a string of IDs, they're written to a temporary 'input.csv' file
- If a yaml file is detected, all options specified are loaded and added to config dict
- Output directory:
	if given, get path and create directory if necessary
	if not given, create a timestamped output directory
- Set up temporary directory, can be configured but by default is $TEMPDIR. If `--no-temp`, set temporary directory to the output directory

### 2) Finding COG-UK data

- If `--CLIMB`, check path to civet directory exists, which verifies the user is on CLIMB
- If `--datadir` is specified and the `--remote` flag is not used, check if it has the appropriate files
- If `--datadir` not given, default to look for a directory called 'civet-cat' in the current working directory
- If `--remote` has been specified, rsync data from CLIMB
	If rsync fails exit with informative error
- If SSH keys not configured, specify `--remote` alongside `-uun` 
- Check background metadata has `data_column` (default: central_sample_id)

### 3) ``-fm / --from-metadata``

- Exit if used in conjunction with input.csv, ID string or fasta options (not compatible yet)
- Find all records in the background metadata that match search criteria and pass on into civet

### 4) QC on input csv
- Check the file exists
- Check the file has `input_column` (default: name)

### 5) QC on input fasta
- If specified, check if input fasta file exists, otherwise exit with informative message 
- Check `--min-length` (Default: 10,000)
- Check `--max-ambiguity` content (Default: 0.5)
- Check if sequence name in the input.csv file
- Write out failed query ids and why they failed (to pass to report)
- Write out sequences that passed QC to pass into pipeline
- If no sequences pass (or no fasta is given), check if any queries in the input.csv match with the database (if not, exit)

### 6) Access installed package data
- Find reference fasta file for downstream datafunk command
- Find polytomies figure for the report
- Find the report template
- Find map files
- If `-sc / --sequencing-centre` specified, check it's a valid centre name and find the appropriate header for the report. If not a real sequencing centre, exit and display valid names

### 7) QC of mapping and report configuration

- Load the report options
- Check if the columns specified exist in the metadata

### 8) Configure local tree size and collapsing

- Check `--up-distance` and `--down-distance`, if either are not specified, use `--distance` (default: 2)
- Check if the columns specified exist in the metadata
- Add `--collapse-threshold` input

### 9) Load the free text specified for the report

- Make the report title
- Parse free text specified in the config file

### 10) Miscellaneous arguments
- Any option not specified as a cmd line argument or in a config file defaults to default values
- If `--verbose`, print shell commands set to True, quiet mode off and custom logger disabled
- If `--launch-browser`, launch grip at end of report generation
- `--threads` added for number of parallel jobs to run
- If `--generate-config`, dump these config options to a file and exit

### 11) Pass options to snakemake
Passing through the configuration dictionary and setting the working directory to the temporary directory (if `--no-temp`, this inherits the output directory).  

### 12) Main analysis pipeline
1) check_cog_db
Check against the COG-UK phylogenetics database for your query id, will check search field, which by default is `central_sample_id`, this finds all sequences that will already be in the big tree. Use `data_column` to configure a custom search field. 
3) get_closest_cog
If you have supplied an additional fasta file with sequences that haven't yet been uploaded to the COG-UK database, this side pipeline finds the closest matching sequence, details in Section 13. This workflow returns a csv of the closest hits within the background phylogeny for each query not already in the tree. This information gets added to the report. It also returns the input query sequences aligned, with the UTRs masked out. 
4) combine_metadata
Metadata from the closest COG sequence is combined with the metadata found directly within the database.
5) prune_out_catchments
Using the names of the closest COG sequences and of the sequences already in the tree, local trees with a radius of `--distance` number of nodes are pulled out of the large global tree. 
6) process_catchments
Another sub-workflow that gets spawned off with the set of catchment trees produced from the prune_out_catchments rule. Full details of this step are found in Section 10, but in brief there are three alternative snakefiles can can be called at this step depending on whether --delay-tree-collapse has been called and depending on whether new sequences need to be added into the tree or not. The output of this step is a directory of "local_trees" that have the nodes we have not queried collapsed to the nearest polytomy, with an output summary file of the tips that had been in that sub-tree prior to collapse.
7) find_snps
A side pipeline that runs `snipit` on your set of query sequences, producing a figure that summarises all the snps relative to the reference sequence. Turn these figures off with `--no-snipit`
8) regional_mapping
A side pipeline that will render maps summarising the variation in the local area surrounding your query sequences of interest. 
9) regional_map_rendering
Render the figures showing local maps

### 9) Making the report

Parses metadata CSVs and any input csv provided, along with tree files generated previously. There are a number of different options for figures, some of which must be specified when calling CIVET. 

The report broadly consists of:

- A list of those sequences which did not meet quality control requirements
- Two summary tables containing information about sequences which were already in the database, and those which were not, but had sequence data provided in fasta file. 
- Trees showing the phylogenetic context of each of the query sequences which passed quality control.

There also a number of optional figures, including different SNP data tables, maps and bar charts describing collpased nodes.

For information on these optional figures and how to configure the core report, see the [report documentation](https://github.com/COG-UK/civet/docs/report_docs.md) and [map documentation](https://github.com/COG-UK/civet/docs/map_option_docs.md) files.


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
