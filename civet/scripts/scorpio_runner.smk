import os
import csv
import collections
from civet.utils.config import *

catchments = config['figure_catchments']
mutations = " ".join(config[KEY_MUTATIONS])

rule all:
    input:
        expand(os.path.join(config[KEY_TEMPDIR],"scorpio","scorpio.{catchment}.csv"), catchment=catchments),
        csv = os.path.join(config[KEY_OUTDIR],"master_metadata.csv")

rule scorpio_type_query:
    input:
        fasta = config["fasta"]
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"scorpio","scorpio.query.csv")
    log: os.path.join(config[KEY_TEMPDIR],"scorpio","scorpio.query.log")
    shell:
        f"scorpio haplotype "
        "-i {input.fasta:q} "
        f"--mutations {mutations} "
        "--append-genotypes "
        "-n mutations "
        "-o {output.csv:q} &> {log:q}"

rule scorpio_type_catchment:
    input:
        fasta = os.path.join(config[KEY_DATA_OUTDIR],"catchments","{catchment}.fasta")
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"scorpio","scorpio.{catchment}.csv")
    log: os.path.join(config[KEY_TEMPDIR],"scorpio","scorpio.{catchment}.log")
    shell:
        f"scorpio haplotype "
        "-i {input.fasta:q} "
        f"--mutations {mutations} "
        "--append-genotypes "
        "-n mutations "
        "-o {output.csv:q} &> {log:q}"

rule gather_mutation_calls:
    input:
        rules.scorpio_type_query.output.csv,
        expand(os.path.join(config[KEY_TEMPDIR],"scorpio","scorpio.{catchment}.csv"), catchment=catchments)
    output:
        csv = os.path.join(config[KEY_OUTDIR],"master_metadata.csv")
    run:
        mut_dict = collections.defaultdict(dict)
        for input_file in input:
            with open(input_file,"r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    for mutation in config[KEY_MUTATIONS]:
                        mut_dict[row[KEY_QUERY]][mutation]= row[mutation]
        with open(output.csv,"w") as fw:
            with open(config["csv"],"r") as f:
                reader = csv.DictReader(f)
                header = reader.fieldnames
                for mutation in config[KEY_MUTATIONS]:
                    if mutation not in header:
                        header.append(mutation)
                writer = csv.DictWriter(fw, fieldnames=header,lineterminator='\n')
                writer.writeheader()
                for row in reader:
                    new_row=row
                    if new_row[KEY_QUERY_BOOLEAN] == "True":
                        for mutation in config[KEY_MUTATIONS]:
                            new_row[mutation] = mut_dict[new_row[KEY_HASH]][mutation]
                    elif new_row["in_tree"] == "True":
                        for mutation in config[KEY_MUTATIONS]:
                            sequence_name = new_row[config[KEY_SEQUENCE_ID_COLUMN]]
                            new_row[mutation] = mut_dict[sequence_name][mutation]
                    else:
                        for mutation in config[KEY_MUTATIONS]:
                            new_row[mutation] = "Not typed"
                    writer.writerow(new_row)