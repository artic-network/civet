from civet.utils.log_colours import green,cyan,red
from civet.analysis_functions import catchment_parsing

import collections
from Bio import SeqIO

"""
To do: 

- input fasta, map to ref trim and pad- take note of insertions
account for ones that dont map

- merge and hash aligned input fasta and extracted matched fasta

- pull out local catchments from hashed sequences (configure distance and shape)

- find a sensible way to merge local catchments

-qs - are we building a tree?
    if yes - downsample local catchments (by what trait, enrich for things, protect things)
           - add an outgroup
           - build the tree with outgroup and downsample

- summarise non-downsampled catchments in the data

"""
rule all:
    input:
        os.path.join(config["outdir"],"catchments.merged.csv")

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
        hash_map = os.path.join(config["tempdir"],"hash_map.csv")
    run:
        records = 0
        
        seq_map = {}
        hash_map = collections.defaultdict(list)
        if config["matched_fasta"]:
            records = catchment_parsing.add_to_hash(config["matched_fasta"],seq_map,hash_map,records)

        if config["query_fasta"]:
            records = catchment_parsing.add_to_hash(input.query_fasta,seq_map,hash_map,records)
        
        with open(output.fasta,"w") as fseqs:
            for key in seq_map:
                fseqs.write(f">{key}\n{seq_map[key]}\n")

        with open(output.hash_map,"w") as fmap:
            fmap.write("key,sequences\n")
            for key in hash_map:
                fmap.write(f"{key},{';'.join(hash_map[key])}\n")

        print(green("Query sequences collapsed from ") + f"{records}" +green(" to ") + f"{len(seq_map)}" + green(" unique sequences."))

rule find_catchment:
    input:
        fasta = rules.seq_brownie.output.fasta
    output:
        catchments = os.path.join(config["outdir"],"catchments.csv")
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
        catchments = rules.find_catchment.output.catchments
    output:
        merged_catchments = os.path.join(config["outdir"],"catchments.merged.csv"),
        catchment_key = os.path.join(config["tempdir"],"catchment_key.csv")
    run:
        merged_catchments = catchment_parsing.get_merged_catchments(input.catchments,output.catchment_key,output.merged_catchments,config)
    
        print(green("Merged into ")+f'{len(merged_catchments)}' + green(" catchments."))
        

