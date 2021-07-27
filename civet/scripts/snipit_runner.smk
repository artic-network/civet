import os
from Bio import SeqIO
import csv
import collections

catchments = config["figure_catchments"]

rule all:
    input:
        expand(os.path.join(config["tempdir"],"snipit","{catchment}.snipit.svg"), catchment=catchments),
        os.path.join(config["tempdir"],"snipit","prompt.txt")

rule make_snipit_alignments:
    input:
        fasta = config["fasta"],
        csv = config["csv"]
    output:
        expand(os.path.join(config["tempdir"],"snipit","{catchment}.aln.fasta"), catchment=catchments)
    run:
        catchment_dict = collections.defaultdict(list)
        with open(input.csv,"r") as f:
            reader = csv.DictReader(f)
            
            if config["input_display_column"] in reader.fieldnames:
                column = config["input_display_column"]
            else:
                column = config["input_id_column"]

            for row in reader:
                if "query_boolean" in reader.fieldnames:
                    query = row["query_boolean"]
                
                if query=="True":
                    if row["source"] == "input_fasta":
                        if row["qc_status"] == "Pass":
                            catchment_dict[row["catchment"]].append((row[column], row["hash"]))
                    else:
                        catchment_dict[row["catchment"]].append((row[column], row["hash"]))

        sequences = {}
        for record in SeqIO.parse(input.fasta,"fasta"):
            sequences[record.id] = record.seq

        for record in SeqIO.parse(config["outgroup_fasta"],"fasta"):
            reference = record
        
        for catchment in catchment_dict:
            with open(os.path.join(config["tempdir"],"snipit",f"{catchment}.aln.fasta"),"w") as fw:
                fw.write(f">{reference.id}\n{reference.seq}\n")

                for query in catchment_dict[catchment]:
                    fw.write(f">{query[0]}\n{sequences[query[1]]}\n")

rule run_snipit:
    input:
        aln = os.path.join(config["tempdir"],"snipit","{catchment}.aln.fasta")
    params:
        out_stem =os.path.join(config["tempdir"],"snipit","{catchment}.snipit")
    output:
        os.path.join(config["tempdir"],"snipit","{catchment}.snipit.svg")
    shell:
        """
        snipit {input.aln:q} -r "outgroup" -o {params.out_stem} -f svg
        """

rule gather_graphs:
    input:
        expand(os.path.join(config["tempdir"],"snipit","{catchment}.snipit.svg"), catchment=catchments)
    output:
        os.path.join(config["tempdir"],"snipit","prompt.txt")
    shell:
        "touch {output:q}"
