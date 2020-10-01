from Bio import Phylo
from Bio import SeqIO
import csv
import collections

rule all:
    input:
        os.path.join(config["outdir"],"cluster_update.csv")

rule get_node_info:
    input:
        tree = config["background_tree"],
        query = config["query"]
    output:
        csv = os.path.join(config["tempdir"], "cluster_civet","node_information.csv")
    shell:
        """
        jclusterfunk 
        """

rule get_common_ancestor:
    input:
        tree = config["background_tree"],
        query = config["query"]
    output:
        ancestor = os.path.join(config["tempdir"], "cluster_civet","common_ancestor.csv")
    run:
        """
        jclusterfunk 
        """

rule find_new_seqs:
    input:
    output:
    run:

rule combine_metadata:
    input:
    output:
    run:
