![](./doc_figures/website_header.png)


## Output

Description of output civet files

<strong>civet_report.md</strong>

A markdown report is generated, with summaries of lineages and genetic diversity present in your query. Trees are visualised in the report and compared to the diversity of lineages present in the community. Optional figures include summary barcharts of collapsed-node composition, genome 'snipit' graphs and maps of the local community diversity. Full report options described [here](./report_docs.md).

<strong>combined_metadata.csv</strong>

A metadata file with the background metadata pulled out for the query sequences of interest, including the metadata for the 'closest' sequences that were found for the optional input fasta file.

<strong>config.final.yaml</strong>

The final config file that gets passed to the report (can be used as a template to run additional analyses).

<strong>query.failed_qc.csv</strong>

If a fasta file is input, this file shows details of any sequences that didn't pass the quality control cut offs and why. For example:

| name | reason_for_failure |
| --- | --- |
| seq1  | fail=seq_len:5052 |
| seq2 | fail=N_content:0.98 |
| seq3 | fail=N_content:1.0 | 
| seq4 | fail=not_in_query_csv | 

<strong>query.post_qc.fasta</strong>

A file of input fasta sequences that passed the quality control cut offs.

<strong>local_trees</strong>

A directory of the collapsed local trees, with corresponding text files of the content of each collapsed node. 

<strong>catchment_trees</strong>

A directory of the catchment trees, with corresponding text files of the content of each collapsed node (will not have the additional fasta input added in at this stage, these are the trees clipped directly out of the global tree). 


<strong>figures</strong>

A directory of figures that get included in the markdown report.


<strong>--no-temp</strong>

With the ``--no-temp`` flag, many more intermediate files are output, including the local tree files prior to collapse. 