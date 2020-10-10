from Bio import Phylo
from Bio import SeqIO
import csv
import collections

prefix = config["output_prefix"]

rule all:
    input:
        os.path.join(config["outdir"],f"{prefix}.csv")

rule find_common_ancestor:
    input:
        tree = config["background_tree"],
        query = config["query"]
    params:
        outdir = os.path.join(config["tempdir"], "cluster_civet")
    output:
        tree = os.path.join(config["tempdir"], "cluster_civet",f"{prefix}_subtree_1.newick"),
        taxa = os.path.join(config["tempdir"], "cluster_civet",f"{prefix}_subtree_1.csv"),
    shell:
        """
        /Users/s1680070/repositories/jclusterfunk/release/jclusterfunk_v0.0.4/jclusterfunk context \
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
        query = config["query"],
        background_metadata= config["background_metadata"]
    output:
        new_metadata = os.path.join(config["outdir"],f"{prefix}.csv")
    run:
        old_cluster = []
        with open(output.new_metadata, "w") as fw:
            with open(input.query, newline="") as f:
                reader = csv.DictReader(f)
                header_names = reader.fieldnames

            
                if "new" not in header_names:
                    header_names.append("new")
                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()
            
                for row in reader:
                    new_row = row
                    new_row["new"] = "False"
                    old_cluster.append(new_row["sequence_name"])
                    writer.writerow(new_row)
                
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
                    if row["country"] == 
                    if row["sequence_name"] in new_sequences:
                        
                        new_row = {}
                        for col in header_names:
                            if col in background_header:
                                new_row[col] = row[col]
                            else:
                                new_row[col] = ""

                        new_row["new"] = "True"
                        writer.writerow(new_row)
            
