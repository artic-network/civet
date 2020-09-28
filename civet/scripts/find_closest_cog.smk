
import os
from Bio import SeqIO
import csv

rule all:
    input:
        os.path.join(config["tempdir"],"closest_cog.csv"),
        os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")


rule non_cog_minimap2_to_reference:
    input:
        fasta = config["to_find_closest"],
        reference = config["reference_fasta"]
    output:
        sam = os.path.join(config["tempdir"],"post_qc_query.reference_mapped.sam")
    shell:
        """
        minimap2 -a -x asm5 {input.reference:q} {input.fasta:q} > {output.sam:q}
        """

rule non_cog_remove_insertions_and_trim_and_pad:
    input:
        sam = rules.non_cog_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        insertions = os.path.join(config["tempdir"],"post_qc_query.insertions.txt")
    output:
        fasta = os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam:q} \
          -r {input.reference:q} \
          -o {output.fasta:q} \
          -t [{config[trim_start]}:{config[trim_end]}] \
          --pad \
          --log-inserts 
        """

rule gofasta_against_cog:
    input:
        query_seqs = rules.non_cog_remove_insertions_and_trim_and_pad.output.fasta,
        cog_seqs = config["seq_db"]
    output:
        csv = os.path.join(config["tempdir"],"closest_gofasta.csv")
    shell:
        # produces query,closest,SNPdistance,SNPs
        """
        gofasta closest --target {input.cog_seqs:q} \
        --query {input.query_seqs:q} \
        --outfile {output.csv:q} \
        -t {workflow.cores}
        """

rule parse_closest:
    input:
        csv = rules.gofasta_against_cog.output.csv,
        metadata = config["cog_metadata"],
        fasta = config["seq_db"]
    output:
        fasta = os.path.join(config["tempdir"],"closest_cog.fasta"),
        csv = os.path.join(config["tempdir"],"closest_cog.csv")
    shell:
        """
        parse_closest.py \
        --csv {input.paf:q} \
        --metadata {input.metadata:q} \
        --seqs {input.fasta:q} \
        --csv-out {output.csv:q} \
        --seqs-out {output.fasta:q} \
        --search-field {config[search_field]}
        """
