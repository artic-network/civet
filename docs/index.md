# civet
**C**luster **I**nvestigation & **V**irus **E**pidemiology **T**ool


<strong>civet</strong> is a tool developed with 'real-time' genomics in mind. 

Using a background phylogeny, such as the large phylogeny available through the COG-UK infrastructure on CLIMB, <strong>civet</strong> will generate a report for a set of sequences of interest i.e. an outbreak investigation. 

If the sequences are already on CLIMB and part of the large tree, <strong>civet</strong> will pull out the local context of those sequences, merging the smaller local trees as appropriate. If sequences haven't yet been incorporated into the large phylogeny, for instance if they have just been sequenced, <strong>civet</strong> will find the closest sequence in the large tree, pull the local tree of that sequence out and add your sequence in. The local trees then get collapsed to display in detail only the sequences of interest so as not to inform investigations beyond what was suggested by epidemiological data. 

A fully customisable report is generated, summarising information about the sequences of interest. The tips of these trees can be coloured by any categorical trait present in the input csv, and additional fields added to the tip labels. Optional figures may be added to describe the local background of UK lineages and to map the query sequences using coordinates, again colourable by a custom trait. 


### civet documentation
  * [Data anonymisation](./safety_level.md)
  * [Install and update civet](./installation.md)
  * [Input options](./input_options.md)
  * [Background data](./background_data.md)
  * [Usage](./usage.md)
  * [Analysis pipeline](./analysis_pipeline.md)
  * [Output](./output.md)
  * [Report options and descriptions](./report_docs.md)
  * [Mapping options and usage instructions](./map_option_docs.md)
  * [Administrative level 2 regions](./geographic_data.md)
  * [Example report](./civet_report_example.md)
  * [Interpretation guide](./interpretation.md)
  * [Contributors & acknowledgements](./acknowledgements.md)
  * [Software versions](./acknowledgements.md)



<strong>civet</strong> was created by Áine O'Toole & Verity Hill, Rambaut Group, Edinburgh University

<img src="./doc_figures/workflow_diagram.png">