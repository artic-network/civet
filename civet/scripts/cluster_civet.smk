from Bio import Phylo
from Bio import SeqIO
import csv
import collections

prefix = config["output_prefix"]

rule all:
    input:
        os.path.join(config["outdir"],f"{prefix}.csv")

rule get_sequence_names:
    input:
        query = config["query"],
        metadata = config["background_metadata"]
    output:
        names = os.path.join(config["outdir"], "cluster_civet","input.sequence_names.csv")
    run:
        ids = {}
        with open(input.query,"r") as f:
            reader = csv.DictReader(f)
            ids["header"] = reader.fieldnames
            for row in reader:
                ids[row[config["input_column"]]]=row

        with open(output.names, "w") as fw:
            
            header = ids["header"]
            header.append("sequence_name")
            if "new" not in header:
                header.append("new")
            writer = csv.DictWriter(fw, fieldnames=header,lineterminator='\n')
            writer.writeheader()
            
            with open(input.metadata,"r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row[config["data_column"]] in ids:
                        new_row = ids[row[config["data_column"]]]
                        new_row["sequence_name"] = row["sequence_name"]
                        new_row["new"] = False
                        writer.writerow(new_row)


rule find_common_ancestor:
    input:
        tree = config["background_tree"],
        query = rules.get_sequence_names.output.names
    params:
        outdir = os.path.join(config["tempdir"], "cluster_civet")
    output:
        tree = os.path.join(config["tempdir"], "cluster_civet",f"{prefix}_subtree_1.newick"),
        taxa = os.path.join(config["tempdir"], "cluster_civet",f"{prefix}_subtree_1.csv")
    shell:
        """
        jclusterfunk context \
        -i "{input.tree}" \
        -o "{params.outdir}" \
        --mrca \
        -f newick \
        -p {config[output_prefix]}_ \
        -m "{input.query}" \
        --output-taxa \
        --id-column sequence_name
        """

rule get_new_query:
    input:
        tree_taxa = rules.find_common_ancestor.output.taxa,
        query = rules.get_sequence_names.output.names,
        background_metadata= config["background_metadata"]
    output:
        new_metadata = os.path.join(config["outdir"],f"{prefix}.csv")
    run:
        old_cluster = []
        with open(output.new_metadata, "w") as fw:
            with open(input.query, newline="") as f:
                reader = csv.DictReader(f)
                header_names = reader.fieldnames

                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()
            
                for row in reader:
                    old_cluster.append(row["sequence_name"])
                    writer.writerow(row)
                
            new_sequences = []

            with open(input.tree_taxa, newline="") as f_tips:
                for l in f_tips:
                    l = l.rstrip('\n')
                    if l not in old_cluster:
                        new_sequences.append(l)
                
            with open(input.background_metadata, newline="") as f:
                reader = csv.DictReader(f)
                background_header = reader.fieldnames
                
                for row in reader:
                    if row["sequence_name"] in new_sequences:
                        
                        new_row = {}

                        for col in header_names:
                            if col == config["input_column"]:
                                new_row[col] = row[config["data_column"]]
                            elif col in background_header:
                                new_row[col] = row[col]
                            else:
                                new_row[col] = ""

                        new_row["new"] = "True"
                        writer.writerow(new_row)
            
