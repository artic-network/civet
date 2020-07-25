"""
Passed into config:

catchment_str=tree_1,tree_2,...,tree_X
"snakemake --nolock --snakefile {input.snakefile_collapse_before:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            # "--directory {params.tempdir:q} "
                            "--config "
                            f"catchment_str={catchment_str} "
                            "outdir={params.outdir:q} "
                            # "tempdir={params.tempdir:q} "
                            "not_cog_csv={input.not_cog_csv:q} "
                            "post_qc_query={input.not_cog_query_seqs:q} "
                            "all_cog_seqs={input.all_cog_seqs:q} "
                            "combined_metadata={input.combined_metadata:q} "
                            "--cores {params.cores}"
"""
from Bio import Phylo
from Bio import SeqIO
import csv
import collections

config["tree_stems"] = config["catchment_str"].split(",")

rule all:
    input:
        expand(os.path.join(config["tempdir"], "collapsed_trees","{tree}.tree"), tree = config["tree_stems"]),
        os.path.join(config["outdir"],"local_trees","collapse_report.txt"),
        expand(os.path.join(config["outdir"],"local_trees","{tree}.tree"), tree = config["tree_stems"])

rule summarise_polytomies:
    input:
        tree = os.path.join(config["tempdir"], "catchment_trees","{tree}.newick"),
        metadata = config["combined_metadata"]
    params:
        tree_dir = os.path.join(config["tempdir"],"catchment_trees")
    output:
        collapsed_tree = os.path.join(config["tempdir"],"collapsed_trees","{tree}.tree"),
        collapsed_information = os.path.join(config["outdir"],"local_trees","{tree}.txt")
    shell:
        """
        clusterfunk focus -i {input.tree:q} \
        -o {output.collapsed_tree:q} \
        --metadata {input.metadata:q} \
        --index-column closest \
        --in-format newick \
        --out-format newick \
        --output-tsv {output.collapsed_information:q}
        """

rule remove_str_for_baltic:
    input:
        tree = os.path.join(config["tempdir"],"collapsed_trees","{tree}.tree")
    output:
        tree = os.path.join(config["tempdir"],"collapsed_trees","{tree}.newick")
    run:
        with open(output.tree,"w") as fw:
            with open(input.tree, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    l = l.replace("'","")
                    fw.write(l)

rule to_nexus:
    input:
        tree = os.path.join(config["tempdir"],"collapsed_trees","{tree}.newick")
    output:
        tree = os.path.join(config["outdir"],"local_trees","{tree}.tree")
    run:
        Phylo.convert(input[0], 'newick', output[0], 'nexus')


rule summarise_processing:
    input:
        collapse_reports = expand(os.path.join(config["outdir"],"local_trees","{tree}.txt"), tree=config["tree_stems"])
    output:
        report = os.path.join(config["outdir"],"local_trees","collapse_report.txt")
    run:
        with open(output.report, "w") as fw:
            for report in input.collapse_reports:
                fn = os.path.basename(report)
                with open(report, "r") as f:
                    for l in f:
                        l = l.rstrip("\n")
                        new_l = f"{fn}\t{l}\n"
                        fw.write(new_l)
