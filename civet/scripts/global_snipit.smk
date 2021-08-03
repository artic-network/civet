
from Bio import Phylo
from Bio import SeqIO
import csv
import collections

rule all:
    input:
        os.path.join(config["outdir"],"report","figures","genome_graph_global_focal.png")

rule snipit_all:
    input:
        aln = config["focal_alignment"]
    params:
        out_path = os.path.join(config["outdir"],"report","figures","genome_graph_global_focal")
    output:
        os.path.join(config["outdir"],"report","figures","genome_graph_global_focal.png")
    shell:
        """
        snipit {input.aln:q} -r "MN908947" -o {params.out_path}
        """

