from Bio import Phylo
from Bio import SeqIO
import csv
import collections

prefix = config["output_prefix"]
today = config["today"]

cluster_file =  f"{prefix}_{today}.csv"

rule all:
    input:
        os.path.join(config["outdir"],cluster_file)

rule find_common_ancestor:
    input:
        tree = config["background_tree"],
        query = config["query"]
    params:
        outdir = os.path.join(config["tempdir"], "cluster_civet")
    output:
        tree = os.path.join(config["tempdir"], "cluster_civet",f"{prefix}_common_ancestor.tree")
    shell:
        """
        ./release/jclusterfunk_v0.0.4/jclusterfunk context \
        -i "{input.tree}" \
        -o "{params.outdir}" \
        --mrca \
        -f newick \
        -p tree_ \
        -m "{input.query}" \
        --id-column sequence_name
        """

rule extract_taxa:
    input:
        tree = rules.find_common_ancestor.output.tree
    output:
        tree_taxa = os.path.join(config["tempdir"], "cluster_civet",f"{prefix}_common_ancestor.csv")
    shell:
        "clusterfunk get_taxa -i {input.collapsed_tree:q} --in-format newick -o {output.tree_taxa:q} --out-format newick"

rule get_new_query:
    input:
        tree_taxa = rules.extract_tips.output.tree_taxa,
        query = config["query"]
    output:
        new_metadata = os.path.join(config["outdir"],cluster_file)
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
                header2 = reader.fieldnames
                
                for row in reader:
                    if row["sequence_name"] in new_sequences:
                        
                        new_row = row
                        new_row["new"] = "True"
                        writer.writerow(new_row)
            
