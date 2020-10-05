
<section id="banner">
    <div class="content">
      <header>
        <h2>Contributors and acknowledgements</h2>
        <p>Contributions, References and Data</p>
      </header>
    </div>
    <span class="image object">
        <img src="./figures/civet_logo.png" alt="" style="max-width:150px"/>
        </span>
</section>


### Authors

`civet` was created by Áine O'Toole & Verity Hill

### Acknowledgements

This project was part-funded by the ARTIC Network 

We acknowledge the hard work from all members of the COG-UK consortium that has gone into generating the data used by `civet`.

`civet` makes use of [`datafunk`](https://github.com/cov-ert/datafunk),  `gofasta`, [`jclusterfunk`](https://github.com/cov-ert/jclusterfunk), [`clusterfunk`](https://github.com/cov-ert/clusterfunk) functions which have been written by members of the Rambaut Lab, specificially Andrew Rambaut, Rachel Colquhoun, JT McCrone, Ben Jackson and Shawn Yu. The local lineages analysis in `civet` was written by Stefan Rooke.

[`baltic`](https://github.com/evogytis/baltic/tree/master/baltic) by Gytis Dudas is used to visualize the trees.

We acknowledge the hard work and open-science of the individual research labs and public health bodies that have made their genome data accessible on GISAID. A detailed table of acknowledgements for the SARS-CoV-2 genome sequences used to build this lineage system is hosted [here](https://github.com/cov-ert/clusterfunk))

### References

[minimap2](https://github.com/lh3/minimap2) 

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
    - pandas==1.1.0
    - pytools=2020.1
    - dendropy=4.4.0
    - tabulate=0.8.7
    
