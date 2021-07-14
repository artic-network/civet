from civet.utils.log_colours import green,cyan,red
from civet.analysis_functions import catchment_parsing
from civet.utils import misc
from civet.report_functions import report
import collections
import sys
import yaml
from Bio import SeqIO
import csv

"""
To do: 

- input fasta, map to ref trim and pad- take note of insertions
account for ones that dont map

- *done* merge and hash aligned input fasta and extracted matched fasta

- *done* pull out local catchments from hashed sequences (configure distance and shape)

- *done* find a sensible way to merge local catchments

-qs - *done* are we building a tree?
    - *done* define analysis/ report content
    if yes - *done* downsample local catchments (by what trait, enrich for things, protect things)
           - *done* add an outgroup
           - *done* build the tree with outgroup 
           - *done* prune out outgroup
    if no, pass to next step: report making
- summarise non-downsampled catchments in the data

"""
rule all:
    input:
        os.path.join(config["outdir"],"master_metadata.csv"),
        os.path.join(config["outdir"],config["output_reports"][0])

rule align_to_reference:
    input:
        reference = config["reference_sequence"]
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
            print(green("Aligning supplied sequences to reference."))
            shell("""
                    minimap2 -a -x asm5 --sam-hit-only --secondary=no -t  {workflow.cores} {input.reference:q} '{config[query_fasta]}' -o {params.sam:q} &> {log:q} 
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
        
        misc.add_col_to_metadata("hash", hash_map_for_metadata, config["query_metadata"], output.csv, config["sequence_id_column"], config)

        config["query_metadata"] = output.csv
        print(green("Query sequences collapsed from ") + f"{records}" +green(" to ") + f"{len(seq_map)}" + green(" unique sequences."))


"""
check_if_int("snp_distance_up",config) 
        check_if_int("snp_distance_down",config)
        check_if_int("snp_distance_side",config)
"""
rule find_catchment:
    input:
        fasta = rules.seq_brownie.output.fasta
    log: os.path.join(config["tempdir"],"logs","updown_top_ranking.txt")
    output:
        txt = os.path.join(config["tempdir"],"updown_ignore.txt"),
        catchments = os.path.join(config["tempdir"],"catchments.csv")
    run:
        with open(output.txt,"w") as fw:
            for i in config["ids"]:
                fw.write(f"{i}\n")
        shell("""gofasta updown topranking \
        -q {input.fasta:q} \
        -t '{config[background_search_file]}' \
        -o {output.catchments:q} \
        --reference '{config[reference_sequence]}' \
        --dist-push \
        --dist-up {config[snp_distance_up]} \
        --dist-down {config[snp_distance_down]} \
        --dist-side {config[snp_distance_side]} \
        --ignore {output.txt:q} &> {log:q}
        """)

rule merge_catchments:
    input:
        catchments = rules.find_catchment.output.catchments,
        csv = rules.seq_brownie.output.csv
    output:
        yaml = os.path.join(config["outdir"],"config.yaml"),
        merged = os.path.join(config["tempdir"],"catchments","catchments_merged.csv"),
        csv = os.path.join(config["tempdir"],"query_metadata.key.csv"),
        catchment_csv = os.path.join(config["tempdir"],"catchments","query_metadata.catchments.csv")
    run:
        catchment_dict, catchment_key, catchment_count = catchment_parsing.get_merged_catchments(input.catchments,output.merged,config)
        config["catchment_count"] = catchment_count

        if config["catchment_count"] == 0:
            sys.stderr.write(cyan(f"Error: no catchments matched in background data file.\n"))
            sys.exit(-1)

        misc.add_col_to_metadata("catchment", catchment_key, input.csv, output.csv, "hash", config)

        config["query_metadata"] = output.csv

        print(green("Merged into ")+f'{catchment_count}' + green(" catchments."))

        catchment_parsing.add_catchments_to_metadata(config["background_metadata"],output.csv,output.catchment_csv,catchment_dict,config)
        with open(output.yaml, 'w') as fw:
            yaml.dump(config, fw) 
        if config["verbose"]:
            print(red("\n**** CONFIG UPDATED ****"))
            for k in sorted(config):
                print(green(f" - {k}: ") + f"{config[k]}")

rule downsample_catchments:
    input:
        fasta  = rules.seq_brownie.output.fasta,
        csv= rules.merge_catchments.output.catchment_csv
    params:
        catchment_dir = os.path.join(config["data_outdir"],"catchments")
    output:
        csv = os.path.join(config["outdir"],"master_metadata.csv")
    run:
        catchment_parsing.downsample_if_building_trees(input.csv, output.csv, config)
        if '3' in config["report_content"]:
            print(green("Writing catchment fasta files."))
            if not os.path.exists(params.catchment_dir):
                os.mkdir(params.catchment_dir)
            catchment_parsing.write_catchment_fasta(output.csv,input.fasta,params.catchment_dir,config)
        
rule tree_building:
    input:
        yaml = rules.merge_catchments.output.yaml,
        csv = rules.downsample_catchments.output.csv,
        snakefile = os.path.join(workflow.current_basedir,"build_catchment_trees.smk")
    output:
        txt = os.path.join(config["tempdir"],"catchments","prompt.txt")
    run:
        if '3' in config["report_content"]:
            print(green("Running tree building pipeline."))
            # spawn off a side snakemake with catchment wildcard that builds an iqtree for each one
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "{config[log_string]} "
                    "--directory {config[tempdir]:q} "
                    "--configfile {input.yaml:q} "
                    "--config csv={input.csv:q} "
                    "--cores {workflow.cores} && touch {output.txt:q}")
        else:
            shell("touch {output.txt:q}")

rule snipit:
    input:
        yaml = rules.merge_catchments.output.yaml,
        csv = rules.merge_catchments.output.catchment_csv,
        fasta = rules.seq_brownie.output.fasta,
        snakefile = os.path.join(workflow.current_basedir,"snipit_runner.smk")
    output:
        txt = os.path.join(config["tempdir"],"snipit","prompt.txt")
    run:
        if '4' in config["report_content"]:
            print(green("Running snipit pipeline."))
            # spawn off a side snakemake with catchment wildcard that builds an iqtree for each one
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "{config[log_string]} "
                    "--directory {config[tempdir]:q} "
                    "--configfile {input.yaml:q} "
                    "--config fasta={input.fasta:q} csv={input.csv:q} "
                    "--cores {workflow.cores} && touch {output.txt:q}")
        else:
            shell("touch {output.txt:q}")


rule render_report:
    input:
        csv = rules.merge_catchments.output.catchment_csv,
        yaml = rules.merge_catchments.output.yaml,
        snipit = rules.snipit.output.txt,
        trees = rules.tree_building.output.txt
    output:
        html = os.path.join(config["outdir"],config["output_reports"][0])
    run:
        with open(input.yaml, 'r') as f:
            config_loaded = yaml.safe_load(f)

        for report_to_generate in config["output_reports"]:
            report_path = os.path.join(config["outdir"],report_to_generate)
            report.make_report(input.csv,report_path,config_loaded)
        
