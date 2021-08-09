
from Bio import Phylo
from Bio import SeqIO
import csv
import collections

rule all:
    input:
        os.path.join(config["tempdir"], "local_trees_seqs", "all_focal_names.txt"),
	os.path.join(config["outdir"],"report","figures","genome_graph_global_focal.png")
	
# if the reference name is not supplied or not in ref_name, it will take the first entry in the focal alignment
rule get_all_sequence_names:
    input:
        alignment = config["focal_alignment"]
    output:
        seq_names = os.path.join(config["tempdir"], "local_trees_seqs", "all_focal_names.txt")
    params: 
        ref_name = config["reference_name"]
    run:
        first = "yes"
	ref_name = ["Ref", "Reference"]

        record_names = []
        fasta_sequences = SeqIO.parse(open(input.alignment), 'fasta')
        with open(output.seq_names, "w") as handle:
            handle.write(f"name,label\n")
            for record in fasta_sequences:
                record_names.append(record.name)
                if params.ref_name and record.name in params.ref_name:
                    first = "no" 
                    handle.write(f"{record.name},Reference\n")
                elif record.name in ref_name and first != "yes" and not params.ref_name:
                    handle.write(f"{record.name},Reference\n")
                elif first == "yes":
                    handle.write(f"{record_names[0]},Reference\n")
                    first = "no"
                else:
                    handle.write(f"{record.name},{record.name}\n")

rule snipit_all:
    input:
        aln = config["focal_alignment"],
	names = rules.get_all_sequence_names.output.seq_names
    params:
        out_path = os.path.join(config["outdir"],"report","figures","genome_graph_global_focal")
    output:
        os.path.join(config["outdir"],"report","figures","genome_graph_global_focal.png")
    shell:
        """
        snipit {input.aln:q} -o {params.out_path} -l {input.names}
        """

