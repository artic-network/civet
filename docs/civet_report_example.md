











![](./figures/DEFAULT.png)


# Cluster investigation



This investigation was started on 2020-10-05

**Authors**: Verity, Áine, Andrew










The following sequences have failed QC:

 - This_seq_is_too_short Sequence too short: only 5052 bases.
 - This_seq_has_lots_of_Ns Sequence has too many Ns: 98.0\% of bases
 - This_seq_is_literally_just_N Sequence has too many Ns: 100.0\% of bases




Showing adm1,HCW_status on the tree.
adm1is shown graphically on the tree using viridisas the colour scheme.
Summarising nodes by country
Mapping sequences using columns adm2 for outer postocdes.


7 queries after QC (5 matched to background database).
2 additional sequences found in CLIMB or provided in a fasta file.
Time fields provided: sample_date,date_2,date_3
Earliest date: 2020-02-27
Latest date: 2020-06-17





**Table 1** | Queries found in COG-UK database.

| Query ID   | Sequence name in tree   | sample_date   | uk_lineage   | lineage   | phylotype   | Tree           |
|:-----------|:------------------------|:--------------|:-------------|:----------|:------------|:---------------|
| EDB3588    | Scotland/EDB3588/2020   | 2020-04-20    | UK1437       | B.1.93    | UK1437_1.10 | tree_subtree_3 |
| EDB2533    | Scotland/EDB2533/2020   | 2020-04-29    | UK5          | B.1.1     | UK5_1.1     | tree_subtree_4 |
| PHEC-1A65C | England/PHEC-1A65C/2020 | 2020-05-04    | UK120        | B         | UK120_1     | tree_subtree_1 |
| PHEC-1AD2A | England/PHEC-1AD2A/2020 | 2020-03-02    | UK120        | B         | UK120_1     | tree_subtree_1 |
| PHEC-1A917 | England/PHEC-1A917/2020 | 2020-02-27    | UK120        | B         | UK120_1     | tree_subtree_1 |


**Table 2** | Queries matched to closest COG-UK sequence using input sequences

| Query ID        | Closest sequence in tree   |   Distance to closest sequence | SNPs         | sample_date   | uk_lineage   | lineage   | phylotype   | Tree           |
|:----------------|:---------------------------|-------------------------------:|:-------------|:--------------|:-------------|:----------|:------------|:---------------|
| EDB129_closest  | Scotland/EDB6546/2020      |                              2 | 666GA;1313TA | 2020-06-05    | NA           | NA        | NA          | tree_subtree_2 |
| EDB129_closestb | Scotland/EDB6546/2020      |                              0 |              | 2020-03-30    | NA           | NA        | NA          | tree_subtree_2 |













> **Tree 1** | 
3 sequences of interest
   
![](./figures/civet_report_make_trees_1.png)
<img src="./figures/civet_report_make_legend_1.png" alt="drawing" style="width:100%;"/>


![](./figures/civet_report_time_plot_1.png)
![](./figures/genome_graph_tree_subtree_1.png)
> **Tree 2** | 
2 sequences of interest
   
![](./figures/civet_report_make_trees_2.png)
<img src="./figures/civet_report_make_legend_1.png" alt="drawing" style="width:100%;"/>


![](./figures/civet_report_time_plot_2.png)
![](./figures/genome_graph_tree_subtree_2.png)
> **Tree 3** | 
1 sequence of interest
   
![](./figures/civet_report_make_trees_3.png)
<img src="./figures/civet_report_make_legend_1.png" alt="drawing" style="width:100%;"/>


![](./figures/genome_graph_tree_subtree_3.png)
> **Tree 4** | 
1 sequence of interest
   
![](./figures/civet_report_make_trees_4.png)
<img src="./figures/civet_report_make_legend_1.png" alt="drawing" style="width:100%;"/>


![](./figures/genome_graph_tree_subtree_4.png)




### Tree background

The following plots describe the data in the collapsed nodes in more detail.
If more than one country was present, the bar chart describes the number of sequences present in each country. 
Where there were 10 options or more, the largest 10 have been taken.
If a UK sequence is present in the collapsed node, it is always shown in the plot.



![](./figures/civet_report_tree_background_1.png)
![](./figures/civet_report_tree_background_2.png)
![](./figures/civet_report_tree_background_3.png)




## Plotting sequences
There are sequences from 3 admin2 regions
This is divided into:
28.57% (2) in EDINBURGH
42.86% (3) in FIFE
28.57% (2) in ABERDEEN



![](./figures/civet_report_map_sequences_1.png)




## Regional-scale background UK lineage mapping
These figures show the background diversity of lineages in the local area to aid with identifying uncommon lineages.
Based on the sample density for submitted sequences with adm2 metadata, **Lothian** was determined to be the focal NHS Health-board.

The below figure visualises the relative proportion of assigned UK-Lineages for samples sampled in **Lothian** for the defined time-frame.
![](./figures/central_map_ukLin.png)

The below figure visualises the relative proportions of assigned UK-Lineages for samples collected in the whole region for the defined time-frame. Plot-size demonstrates relative numbers of sequences across given NHS healthboards.
![](./figures/region_map_ukLin.png)

Tabulated lineage data for the **central** health-board region:

### Lothian

| UK Lineage   |   Count | Global Lineage   | Count   |
|:-------------|--------:|:-----------------|:--------|
| UK109        |     275 | B.1              | 396.0   |
| UK199        |     109 | B.1.5            | 144.0   |
| UK5          |     105 | B.1.1            | 131.0   |
| UK1611       |      54 | B                | 63.0    |
| UK2464       |      50 | B.1.1.4          | 53.0    |
| UK1531       |      49 | B.1.1.14         | 48.0    |
| UK2913       |      43 | B.1.11           | 42.0    |
| UK21         |      38 | B.1.40           | 30.0    |
| UK72         |      28 | B.1.1.12         | 23.0    |
| UK44         |      27 | B.2              | 23.0    |
| UK5676       |      19 |                  |         |
| UK945        |      19 |                  |         |
| UK512        |      18 |                  |         |
| UK1437       |      17 |                  |         |
| UK107        |      16 |                  |         |
| UK740        |      15 |                  |         |
| UK43         |      15 |                  |         |
| UK939        |      14 |                  |         |
| UK175        |      14 |                  |         |
| UK436        |      12 |                  |         |

Tabulated lineage data for the **neighbouring** health-board regions:

## Borders

| UK Lineage   |   Count | Global Lineage   | Count   |
|:-------------|--------:|:-----------------|:--------|
| UK175        |      83 | B.1.89           | 66.0    |
| UK199        |      31 | B.1.5.10         | 28.0    |
| UK5          |       7 | B.1.71           | 16.0    |
| UK133        |       4 | B.1              | 13.0    |
| UK5676       |       3 | B.1.1            | 8.0     |
| UK1437       |       3 | B.1.5            | 5.0     |
| UK2464       |       2 | B.2              | 4.0     |
| UK1254       |       1 | B.1.1.1          | 4.0     |
| UK21         |       1 | B.1.93           | 3.0     |
| UK889        |       1 | B.2.6            | 1.0     |
| UK285        |       1 |                  |         |
| UK1238       |       1 |                  |         |
| UK786        |       1 |                  |         |
| UK2916       |       1 |                  |         |
| UK190        |       1 |                  |         |
| UK969        |       1 |                  |         |
| UK1126       |       1 |                  |         |
| UK109        |       1 |                  |         |
| UK2913       |       1 |                  |         |
| UK5098       |       1 |                  |         |

### Fife

| UK Lineage   |   Count | Global Lineage   | Count   |
|:-------------|--------:|:-----------------|:--------|
| UK5          |      15 | B.1              | 13.0    |
| UK1437       |      10 | B.1.93           | 10.0    |
| UK175        |       5 | B.1.1.1          | 10.0    |
| UK107        |       3 | B.1.1            | 6.0     |
| UK2464       |       2 | B.2              | 4.0     |
| UK183        |       2 | B.2.1            | 3.0     |
| UK2916       |       2 | B                | 3.0     |
| UK247        |       1 | B.1.5            | 2.0     |
| UK199        |       1 | B.1.1.12         | 2.0     |
| UK958        |       1 | B.1.1.14         | 1.0     |
| UK258        |       1 |                  |         |
| UK667        |       1 |                  |         |
| UK763        |       1 |                  |         |
| UK600        |       1 |                  |         |
| UK2913       |       1 |                  |         |
| UK5676       |       1 |                  |         |
| UK191        |       1 |                  |         |
| UK740        |       1 |                  |         |
| UK5498       |       1 |                  |         |
| UK55         |       1 |                  |         |

### Forth_Valley

| UK Lineage   |   Count | Global Lineage   | Count   |
|:-------------|--------:|:-----------------|:--------|
| UK5098       |      21 | B.1              | 44.0    |
| UK5676       |      19 | B.2              | 19.0    |
| UK1437       |      11 | B                | 16.0    |
| UK2464       |      10 | B.1.1            | 13.0    |
| UK5          |       9 | B.1.93           | 11.0    |
| UK167        |       5 | B.1.5            | 6.0     |
| UK120        |       5 | B.1.1.14         | 6.0     |
| UK44         |       5 | B.1.1.20         | 3.0     |
| UK14         |       4 | A.1              | 3.0     |
| UK939        |       3 | A.2              | 2.0     |
| UK1289       |       3 |                  |         |
| UK552        |       3 |                  |         |
| UK1040       |       3 |                  |         |
| UK271        |       3 |                  |         |
| UK109        |       3 |                  |         |
| UK2913       |       2 |                  |         |
| UK317        |       2 |                  |         |
| UK39         |       2 |                  |         |
| UK40         |       2 |                  |         |
| UK1074       |       2 |                  |         |

### Lanarkshire

| UK Lineage   |   Count | Global Lineage   | Count   |
|:-------------|--------:|:-----------------|:--------|
| UK1437       |      65 | B.1.93           | 64.0    |
| UK5098       |      40 | B.1              | 55.0    |
| UK5          |      20 | B.1.1            | 24.0    |
| UK1165       |      15 | B.1.1.10         | 15.0    |
| UK100        |       8 | B.1.1.1          | 11.0    |
| UK2913       |       8 | B.1.11           | 8.0     |
| UK191        |       7 | B.1.101          | 8.0     |
| UK175        |       7 | B.1.77           | 7.0     |
| UK669        |       6 | B.2              | 7.0     |
| UK5676       |       6 | B.1.5            | 4.0     |
| UK2464       |       5 |                  |         |
| UK939        |       3 |                  |         |
| UK1230       |       3 |                  |         |
| UK39         |       3 |                  |         |
| UK1187       |       2 |                  |         |
| UK1531       |       2 |                  |         |
| UK261        |       2 |                  |         |
| UK14         |       2 |                  |         |
| UK944        |       1 |                  |         |
| UK199        |       1 |                  |         |







## Appendix

This report summarises the information provided by whole genome sequencing of SARS-COV-2 generated by the COG consortium. 
It is intended to provide an additional layer of analysis for infection control efforts, and to aid in the investigation of outbreak clusters.

For each query sequence, CIVET either finds them in the COG database, or matches them as closely as possible to a sequence in the COG database, and puts them into a UK lineage.

Key points for interpreting this information:

 - This type of analysis is not able to infer direct transmission between two samples. Even identical sequences may be unrelated as SARS-COV2 is relatively slow evolving for an RNA virus. Previous analysis has shown that samples taken over 100 days apart can be identical. 
 - UK lineage and UK phylotype designations are not yet stable, so they can change with each build of the COG-UK phylogeny.
 - If sequences have different global or UK lineage designations, within the same analysis/report, we can rule out close epidemiological linkage.
 - If sequences have different phylotypes, within the same analysis/report, it’s very unlikely that they are direct transmissions. 
 - If sequences share the same lineage and the same phylotype, within the same analysis/report, transmission cannot be ruled out and also cannot be confirmed.


The figure below shows the distribution of time differences that two sequences can be sampled and still be identical. 
It is to illustrate that identical sequences does not confirm linked cases.



![](./figures/polytomies.png)



### Definitions

*Phylotype* 

Each lineage phylogeny is labelled with phylotypes that describe shared mutations in the tree. If two sequences have the same phylotype it means the share mutations. They may also have additional, unique mutations. So having the same phylotype doesn't mean the seqeunces are identical. If sequences have different phylotypes however it means they are present on distinct parts of the phylogenetic tree.

*UK lineage* 

UK lineages are an approximation to distinct introductions of SARS-CoV-2 to the UK based on the phylogenetic tree.

*Global lineage* 

Assigned using the pangolin software, these are phylogenetic lineages. More information can be found at https://github.com/hCoV-2019/lineages
### Software versions

This report was made using:



Python 3.6.11
Matplotlib version 3.2.1
Pandas version 1.1.0
Tabulate version 0.8.7
CSV version 1.0
Numpy version 1.18.5
Scipy version 1.5.1
No version number for Baltic
COG data is now submitted every day, so the background data was updated yesterday
CIVET version is 0.1


## Acknowledgements

This report was generated by CIVET, made primarily by Áine O'Toole and Verity Hill, using code from Rambaut Lab members.

The background data from the UK was generated by the COG consortium (https://www.cogconsortium.uk/), a national, multi-centre consortium for the sequencing and analysis of SARS-CoV-2 genomes for Public Health.

We also use some background data from GISAID (https://www.gisaid.org/) in the phylogenies. We thank everyone involved in the global sequencing effort for making their data available. 

Tree data was visualised using baltic (https://github.com/evogytis/baltic)

Mapping data was downloaded from the Global Administrative Database (https://gadm.org/) and Natural Earth (https://www.naturalearthdata.com/)



![](./figures/footer.png)

