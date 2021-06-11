from civet.utils.log_colours import green,cyan,red
from civet.analysis_functions import catchment_parsing
from civet.utils import misc

import collections
import sys
import yaml
from Bio import SeqIO

"""
To do: 

- input fasta, map to ref trim and pad- take note of insertions
account for ones that dont map

- *done* merge and hash aligned input fasta and extracted matched fasta

- *done* pull out local catchments from hashed sequences (configure distance and shape)

- *done* find a sensible way to merge local catchments

-qs - are we building a tree?
    - define analysis/ report content
    if yes - downsample local catchments (by what trait, enrich for things, protect things)
           - add an outgroup
           - build the tree with outgroup and downsample

- summarise non-downsampled catchments in the data

"""
rule all:
    input:
        os.path.join(config["data_outdir"],"catchments","query_metadata.catchments.csv"),
        os.path.join(config["tempdir"],"catchments","tree.txt")

rule align_to_reference:
    input:
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        sam = os.path.join(config["tempdir"],"mapped.sam")
    output:
        fasta = os.path.join(config["tempdir"],"query.aln.fasta")
    log:
        os.path.join(config["tempdir"], "logs/minimap2_sam.log")
    run:
        if config["query_fasta"]:
            print(green("Aligning supplied fasta to reference."))
            shell("""
                    minimap2 -a -x asm5 --sam-hit-only --secondary=no -t  {workflow.cores} {input.reference:q} '{config[fasta]}' -o {params.sam:q} &> {log:q} 
                    gofasta sam toMultiAlign \
                        -s {params.sam:q} \
                        -t {workflow.cores} \
                        --reference {input.reference:q} \
                        --trimstart {params.trim_start} \
                        --trimend {params.trim_end} \
                        --trim \
                        --pad > '{output.fasta}'
                    """)
        else:
            shell("touch {output.fasta:q}")

rule seq_brownie:
    input:
        query_fasta = os.path.join(config["tempdir"],"query.aln.fasta")
    output:
        fasta = os.path.join(config["tempdir"],"hashed.aln.fasta"),
        csv = os.path.join(config["tempdir"],"metadata.seq_brownie.master.csv")
    run:
        records = 0
        
        seq_map = {}
        hash_map = collections.defaultdict(list)
        hash_map_for_metadata = {}

        if config["matched_fasta"]:
            records = catchment_parsing.add_to_hash(config["matched_fasta"],seq_map,hash_map,records)

        if config["query_fasta"]:
            records = catchment_parsing.add_to_hash(input.query_fasta,seq_map,hash_map,records)
        
        with open(output.fasta,"w") as fseqs:
            for key in seq_map:
                fseqs.write(f">{key}\n{seq_map[key]}\n")

            for key in hash_map:
                for seq in hash_map[key]:
                    hash_map_for_metadata[seq] = key
        
        misc.add_col_to_metadata("hash", hash_map_for_metadata, config["query_metadata"], output.csv, config["fasta_column"], config)

        config["query_metadata"] = output.csv
        print(green("Query sequences collapsed from ") + f"{records}" +green(" to ") + f"{len(seq_map)}" + green(" unique sequences."))

rule find_catchment:
    input:
        fasta = rules.seq_brownie.output.fasta
    output:
        catchments = os.path.join(config["tempdir"],"catchments.csv")
    shell:
        """
        gofasta updown topranking \
        -q {input.fasta:q} \
        -t '{config[background_search_file]}' \
        -o {output.catchments:q} \
        --reference '{config[reference_fasta]}' \
        --size-total {config[catchment_size]}
        """

rule merge_catchments:
    input:
        catchments = rules.find_catchment.output.catchments,
        csv = rules.seq_brownie.output.csv
    params:
        catchment_dir = os.path.join(config["data_outdir"],"catchments")
    output:
        yaml = os.path.join(config["tempdir"],"catchments","config.yaml"),
        merged = os.path.join(config["tempdir"],"catchments","catchments_merged.csv"),
        csv = os.path.join(config["tempdir"],"query_metadata.key.csv"),
        catchment_csv = os.path.join(config["data_outdir"],"catchments","query_metadata.catchments.csv")
    run:
        catchment_dict, catchment_key, catchment_count = catchment_parsing.get_merged_catchments(input.catchments,output.merged,config)
        config["catchment_count"] = catchment_count

        if config["catchment_count"] == 0:
            sys.stderr.write(cyan(f"Error: no catchments matched in background data file.\n"))
            sys.exit(-1)

        misc.add_col_to_metadata("catchment", catchment_key, input.csv, output.csv, "hash", config)

        config["query_metadata"] = output.csv

        print(green("Merged into ")+f'{catchment_count}' + green(" catchments."))

        catchment_parsing.add_catchments_to_metadata(config["background_csv"],output.csv,output.catchment_csv,catchment_dict,config)
        
        if '3' in config["report_content"]:
            with open(output.yaml, 'w') as fw:
                yaml.dump(config, fw) #so at the moment, every config option gets passed to the report
            print("Writing catchment fasta files.")
            catchment_parsing.write_catchment_fasta(catchment_dict,params.catchment_dir,config)

"""
rule downsampling:
    input:

    output:

    run:
        optional downsampling of catchments, with a protection/ enrichment metric 
        that can prevent certain sequences being removed
"""

rule tree_building:
    input:
        csv = rules.merge_catchments.output.csv,
        yaml = rules.merge_catchments.output.yaml,
        snakefile = os.path.join(workflow.current_basedir,"build_catchment_trees.smk")
    output:
        txt = os.path.join(config["tempdir"],"catchments","tree.txt")
    run:
        if '3' in config["report_content"]:
            print(green("Running tree building pipeline."))
            # spawn off a side snakemake with catchment wildcard that builds an iqtree for each one
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "{config[log_string]} "
                    "--directory {config[tempdir]:q} "
                    "--configfile {input.yaml:q} "
                    "--cores {workflow.cores} && touch {output.txt}")

"""
rule render_report:
    input:

    output:

    run:
        take processed metadata
        take treefiles
        render report using mako
"""