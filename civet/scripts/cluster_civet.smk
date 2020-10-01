from Bio import Phylo
from Bio import SeqIO
import csv
import collections

rule all:
    input:
        os.path.join(config["outdir"],"cluster_update.csv")

rule find_common_ancestor:
    input:
        tree = config["background_tree"],
        query = config["query"]
    params:
        outdir = os.path.join(config["tempdir"], "cluster_civet")
    output:
        tree = os.path.join(config["tempdir"], "cluster_civet","common_ancestor.tree")
    shell:
        """
        jclusterfunk context \
        -i "{input.tree}" \
        -o "{params.outdir}" \
        --max-parent 1 \
        --max-child 1000 \
        -f newick \
        -p tree_ \
        --ignore-missing \
        -m "{input.query}" \
        --id-column {config[input_column]} \
        && touch "{output.csv}" 
        """

rule extract_tips:
    input:
        tree = rules.find_common_ancestor.output.tree
    output:
        tips = os.path.join(config["tempdir"], "cluster_civet","common_ancestor.csv")
    shell:
        """
        jclusterfunk metadata
        """

rule get_new_query:
    input:
        tips = rules.extract_tips.output.tips,
        query = config["query"]
    output:
        new_metadata = os.path.join(config["outdir"], "cluster_civet","updated_cluster_query.csv")
    run:
        old_cluster = []

        with open(input.query, newline="") as f:
            reader = csv.DictReader(f)
            header_names = reader.fieldnames

            with open(output.new_metadata, "w") as fw:
                if "new_sequence" not in header_names:
                    header_names.append("new_sequence")
                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()
            
                for row in reader:
                    new_row = row
                    new_row["new_sequence"] = "False"
                    old_cluster.append(new_row["sequence_name"])
                    writer.writerow(new_row)
                
                with open(input.tips, newline="") as f_tips:
                    reader_tips = csv.DictReader(f_tips)
                    c = 0
                    for row in reader_tips:
                        if row["sequence_name"] not in old_cluster:
                            new_row = row
                            new_row["new_sequence"] = "True"
                            writer.writerow(new_row)
                            c+=1

                    print(qcfunk.green("New sequences detected in cluster:") + f"{c}")
                    config["update_cluster"] = True
                    