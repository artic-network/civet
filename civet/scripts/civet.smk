from civet.utils.log_colours import green,cyan,red
from civet.analysis_functions import catchment_parsing
from civet.utils import misc
from civet.report_functions import report
import collections
import sys
import yaml
from Bio import SeqIO
import csv

from civet.utils.config import *


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
        os.path.join(config[KEY_OUTDIR],"master_metadata.csv"),
        os.path.join(config[KEY_OUTDIR],config[KEY_OUTPUT_REPORTS][0])

rule align_to_reference:
    input:
        reference = config[KEY_REFERENCE_SEQUENCE]
    params:
        trim_start = config[KEY_TRIM_START],
        trim_end = config[KEY_TRIM_END],
        sam = os.path.join(config[KEY_TEMPDIR],"mapped.sam")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"query.aln.fasta")
    log: os.path.join(config[KEY_TEMPDIR], "logs/minimap2_sam.log")
    run:
        if config[KEY_QUERY_FASTA]:
            print(green("Aligning supplied sequences to reference."))
            shell("""
                    minimap2 -a -x asm20 --score-N=0 --sam-hit-only --secondary=no -t  {workflow.cores} {input.reference:q} '{config[query_fasta]}' -o {params.sam:q} &> {log:q} 
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
        query_fasta = os.path.join(config[KEY_TEMPDIR],"query.aln.fasta")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"hashed.aln.fasta"),
        csv = os.path.join(config[KEY_TEMPDIR],"metadata.seq_brownie.master.csv")
    run:
        records = 0
        
        seq_map = {}
        hash_map = collections.defaultdict(list)
        hash_map_for_metadata = {}

        if config[KEY_MATCHED_FASTA]:
            records = catchment_parsing.add_to_hash(config[KEY_MATCHED_FASTA],seq_map,hash_map,records)

        if config[KEY_QUERY_FASTA]:
            records = catchment_parsing.add_to_hash(input.query_fasta,seq_map,hash_map,records)
        
        with open(output.fasta,"w") as fseqs:
            for key in seq_map:
                fseqs.write(f">{key}\n{seq_map[key]}\n")

        for hash_str in hash_map:
            for record_id in hash_map[hash_str]:
                hash_map_for_metadata[record_id] = hash_str
                
        if config[KEY_QUERY_FASTA]:
            misc.add_col_to_metadata(KEY_HASH, hash_map_for_metadata, config[KEY_QUERY_METADATA], output.csv, config["input_id_column"], config)
        elif config[KEY_MATCHED_FASTA]:
            misc.add_col_to_metadata(KEY_HASH, hash_map_for_metadata, config[KEY_QUERY_METADATA], output.csv, config["sequence_id_column"], config)

        config[KEY_QUERY_METADATA] = output.csv
        print(green("Query sequences collapsed from ") + f"{records}" +green(" to ") + f"{len(seq_map)}" + green(" unique sequences."))
            

"""
check_if_int("snp_distance_up",config) 
        check_if_int("snp_distance_down",config)
        check_if_int("snp_distance_side",config)
"""
rule find_catchment:
    input:
        fasta = rules.seq_brownie.output.fasta
    log: os.path.join(config[KEY_TEMPDIR],"logs","updown_top_ranking.txt")
    output:
        txt = os.path.join(config[KEY_TEMPDIR],"updown_ignore.txt"),
        catchments = os.path.join(config[KEY_TEMPDIR],"catchments.csv")
    run:
        with open(output.txt,"w") as fw:
            for i in config[KEY_IDS]:
                fw.write(f"{i}\n")
        shell("""gofasta updown topranking \
        -q {input.fasta:q} \
        -t '{config[background_search_file]}' \
        -o {output.catchments:q} \
        --reference '{config[reference_sequence]}' \
        --dist-push {config[push_distance]} \
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
        yaml = os.path.join(config[KEY_OUTDIR],"config.yaml"),
        merged = os.path.join(config[KEY_TEMPDIR],"catchments","catchments_merged.csv"),
        csv = os.path.join(config[KEY_TEMPDIR],"query_metadata.key.csv"),
        catchment_csv = os.path.join(config[KEY_TEMPDIR],"catchments","query_metadata.catchments.csv")
    run:
        catchment_dict, catchment_key, catchment_count = catchment_parsing.get_merged_catchments(input.catchments,output.merged,config)
        config[KEY_CATCHMENT_COUNT] = catchment_count

        if config[KEY_CATCHMENT_COUNT] == 0:
            sys.stderr.write(cyan(f"Error: no catchments matched in background data file.\n"))
            sys.exit(-1)

        misc.add_col_to_metadata("catchment", catchment_key, input.csv, output.csv, KEY_HASH, config)

        config[KEY_QUERY_METADATA] = output.csv

        print(green("Merged into ")+f'{catchment_count}' + green(" catchments."))

        catchment_parsing.add_catchments_to_metadata(config[KEY_BACKGROUND_METADATA],output.csv,output.catchment_csv,catchment_dict,config)
        
        catchment_parsing.which_catchments_too_large(output.catchment_csv,config)

        with open(output.yaml, 'w') as fw:
            yaml.dump(config, fw) 

        if config[KEY_VERBOSE]:
            print(red("\n**** CONFIG UPDATED ****"))
            for k in sorted(config):
                print(green(f" - {k}: ") + f"{config[k]}")

rule downsample_catchments:
    input:
        fasta  = rules.seq_brownie.output.fasta,
        csv= rules.merge_catchments.output.catchment_csv,
        yaml = os.path.join(config[KEY_OUTDIR],"config.yaml")
    params:
        catchment_dir = os.path.join(config["data_outdir"],"catchments")
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"catchments","catchment_metadata.csv")
    run:
        with open(input.yaml, 'r') as f:
            config_loaded = yaml.safe_load(f)

        catchment_parsing.downsample_if_building_trees(input.csv, output.csv, config_loaded)
        if '3' in config[KEY_REPORT_CONTENT]:
            print(green("Writing catchment fasta files."))
            if not os.path.exists(params.catchment_dir):
                os.mkdir(params.catchment_dir)
            catchment_parsing.write_catchment_fasta(output.csv,input.fasta,params.catchment_dir,config_loaded)


rule scorpio_type:
    input:
        snakefile = os.path.join(workflow.current_basedir,"scorpio_runner.smk"),
        yaml = os.path.join(config[KEY_OUTDIR],"config.yaml"),
        fasta = rules.seq_brownie.output.fasta,
        csv = rules.downsample_catchments.output.csv
    params:
        muts = os.path.join(config[KEY_TEMPDIR],"scorpio","scorpio.mutations.csv")
    output:
        csv = os.path.join(config[KEY_OUTDIR],"master_metadata.csv")
    log: os.path.join(config[KEY_TEMPDIR],"scorpio","scorpio.log")
    run:
        if config[KEY_MUTATIONS]:
            print(green("Running scorpio typing."))
            # spawn off a side snakemake with catchment wildcard that types the mutations of interest per catchment
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "{config[log_string]} "
                    "--directory {config[tempdir]:q} "
                    "--configfile {input.yaml:q} "
                    "--config fasta={input.fasta:q} csv={input.csv:q} "
                    "--cores {workflow.cores} &> {log:q}")
                    
        else:
            shell("cp {input.csv:q} {output.csv:q}")

rule tree_building:
    input:
        yaml = rules.merge_catchments.output.yaml,
        csv = rules.scorpio_type.output.csv,
        snakefile = os.path.join(workflow.current_basedir,"build_catchment_trees.smk")
    output:
        txt = os.path.join(config[KEY_TEMPDIR],"catchments","prompt.txt")
    run:
        if '3' in config[KEY_REPORT_CONTENT]:
            print(green("Running tree building pipeline."))
            # spawn off a side snakemake with catchment wildcard that builds an iqtree for each one
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "{config[log_string]} "
                    "--directory {config[tempdir]:q} "
                    "--configfile {input.yaml:q} "
                    f"--config csv='{input.csv}' "
                    "--cores {workflow.cores} && touch {output.txt:q}")
        else:
            shell("touch {output.txt:q}")

rule snipit:
    input:
        yaml = rules.merge_catchments.output.yaml,
        csv = rules.scorpio_type.output.csv,
        fasta = rules.seq_brownie.output.fasta,
        snakefile = os.path.join(workflow.current_basedir,"snipit_runner.smk")
    output:
        txt = os.path.join(config[KEY_TEMPDIR],"snipit","prompt.txt")
    run:
        if '4' in config[KEY_REPORT_CONTENT]:
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


rule global_snipit:
    input:
        fasta = rules.align_to_reference.output.fasta,
	yaml = rules.merge_catchments.output.yaml,
	prompt = rules.snipit.output.txt,
        snakefile = os.path.join(workflow.current_basedir,"global_snipit.smk")
    output:
        txt = os.path.join(config[KEY_TEMPDIR],"global_snipit","prompt.txt")
    run:
        if config[KEY_GLOBAL_SNIPIT]:
            print(green("Running global snipit"))
            # spawn off global snipit
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "{config[log_string]} "
                    "--directory {config[tempdir]:q} "
		    "--configfile {input.yaml:q} "
                    "--config fasta={input.fasta:q} "
                    "--cores {workflow.cores} && touch {output.txt:q}")
        else:
            shell("touch {output.txt:q}")


rule render_report:
    input:
        csv = rules.scorpio_type.output.csv,
        yaml = rules.merge_catchments.output.yaml,
        snipit = rules.snipit.output.txt,
        trees = rules.tree_building.output.txt,
	global_snipit = rules.global_snipit.output.txt
    output:
        html = os.path.join(config[KEY_OUTDIR],config[KEY_OUTPUT_REPORTS][0])
    run:
        with open(input.yaml, 'r') as f:
            config_loaded = yaml.safe_load(f)

        for report_to_generate in config[KEY_OUTPUT_REPORTS]:
            report_path = os.path.join(config[KEY_OUTDIR],report_to_generate)
            report.make_report(input.csv,report_path,config_loaded)
        
