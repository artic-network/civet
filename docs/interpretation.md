![](./doc_figures/website_header.png)

# Interpretation


For each query sequence, CIVET either finds them in the COG database, or matches them as closely as possible to a sequence in the COG database, and puts them into a UK lineage.

### Frequently asked questions


 <strong>What does it mean if two samples have different global or UK lineages?</strong>
 - If sequences have different global or UK lineage designations, within the same analysis/report, we can rule out close epidemiological linkage.
 - Bear in mind that these lineage systems are heirarchical, so if the differencce is a lineage nested within the parent there may be very little genetic distance between the sequences. For example: B.1.24 and B.1.24.1 may in fact be quite closely linked in the phylogenetic tree (which will be reflected in the figure shown in the report), whereas A.1 and B.1.1.1.32 are very distant and transmission can be ruled out in this case.

<strong>What does it mean if two samples have a different phylotype?</strong>
 - If sequences have very different phylotypes, it’s unlikely that they are direct transmissions. 
 - However, phylotypes are assigned heirarchically so nested phylotypes may only be one mutation away from one another. Look at the snipit graph to inform your interpretation of different but similar phylotypes.

 <strong>What is they have the same phylotype?</strong>
 - If sequences share the same lineage and the same phylotype, within the same analysis/report, transmission cannot be ruled out and also cannot be confirmed.


 <strong>Why have the lineage designations changed?</strong>
 - UK lineage and UK phylotype designations are not stable, so they can change with each build of the COG-UK phylogeny. 
 - The phylogenetic tree changes each day with new data and there is an associated uncertainty with every tree produced. 

 <strong>Is this a transmission event?</strong>
 - This type of analysis is not able to infer direct transmission between two samples. Even identical sequences may be unrelated as SARS-COV2 is relatively slow evolving for an RNA virus. Previous analysis has shown that samples taken over 100 days apart can be identical (see Figure 1). 
 - If there is convincing genetic evidence that the viruses are very different, we can rule out transmission however. 

![](./doc_figures/polytomies.png)
<strong>Figure 1: Polytomies in SARS-CoV-2 phylogeny. </strong> Shows the distribution of difference in sampling time between two sequences that are still identical.

### Glossary of terms

*SNP or Single Nucleotide Polymorphism*

This term describes a point mutation e.g. A -> T or G -> C. SNPs are often reported with the position and the change: T1234G or 1234TG describes a mutation from the ancestral variant T to a G at position 1234.   

*Polytomy*

An internal node in the phylogenetic tree that has more than two descendants. The global SARS-CoV-2 phylogeny has many polytomies, meaning many similar or identicial sequences exist in the tree. This is a result of the intensity of sequencing over a relatively short period of time and the mutation rate of SARS-CoV-2.  

*Phylotype* 

Each lineage phylogeny is labelled with phylotypes that describe shared mutations in the tree. If two sequences have the same phylotype it means the share mutations. They may also have additional, unique mutations. So having the same phylotype doesn't mean the seqeunces are identical. If sequences have different phylotypes however it means they are present on distinct parts of the phylogenetic tree.

*UK lineage* 

UK lineages are an approximation to distinct introductions of SARS-CoV-2 to the UK based on the phylogenetic tree.

*Global lineage* 

Assigned using the pangolin software, these are phylogenetic lineages. More information can be found at https://github.com/hCoV-2019/lineages