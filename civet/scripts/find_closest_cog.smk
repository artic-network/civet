"""
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force} "
                        "{params.quiet_mode} "
                        "--directory {params.tempdir:q} "
                        "--config "
                        "outdir={params.outdir:q} "
                        # "tempdir={params.tempdir:q} "
                        # "seq_db={input.seq_db:q} "
                        # "to_find_closest={output.combined_query} "
                        "in_all_cog_metadata={input.in_all_cog_metadata:q} "
                        # "search_field={params.search_field} "
                        "cog_seqs={input.cog_seqs:q} "
                        # "trim_start={params.trim_start} "
                        # "trim_end={params.trim_end} "
                        # "reference_fasta={input.reference_fasta:q} "
                        # "cog_metadata={input.cog_metadata:q} "
                        "--cores {params.cores}")
"""

import os
from Bio import SeqIO
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
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = os.path.join(config["tempdir"],"post_qc_query.insertions.txt")
    output:
        fasta = os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam:q} \
          -r {input.reference:q} \
          -o {output.fasta:q} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts &
        mv insertions.txt {params.insertions:q}
        """

rule minimap2_against_cog:
    input:
        query_seqs = rules.non_cog_remove_insertions_and_trim_and_pad.output.fasta,
        cog_seqs = config["seq_db"]
    threads: workflow.cores
    output:
        paf = os.path.join(config["tempdir"],"post_qc_query.cog_mapped.paf")
    shell:
        """
        minimap2 -x asm5  -t {threads} --secondary=no --paf-no-hit {input.cog_seqs:q} {input.query_seqs:q} > {output.paf:q}
        """

rule parse_paf:
    input:
        paf = rules.minimap2_against_cog.output.paf,
        metadata = config["cog_metadata"],
        fasta = config["seq_db"]
    params:
        search_field = config["search_field"]
    output:
        fasta = os.path.join(config["tempdir"],"closest_cog.fasta"),
        csv = os.path.join(config["tempdir"],"closest_cog.csv")
    shell:
        """
        parse_paf.py \
        --paf {input.paf:q} \
        --metadata {input.metadata:q} \
        --seqs {input.fasta:q} \
        --csv-out {output.csv:q} \
        --seqs-out {output.fasta:q} \
        --search-field {params.search_field}
        """
