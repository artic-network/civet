
from Bio import Phylo
from Bio import SeqIO
import csv
import collections

rule all:
    input:
        os.path.join(config["outdir"], "report", "global_snipit_labels.txt"),
	os.path.join(config["outdir"],"report","figures","genome_graph_global_focal.png")
	
# if the reference name is not supplied or not in ref_name, it will take the first entry in the focal alignment
rule get_all_sequence_names:
    input:
        alignment = config["focal_alignment"]
    output:
        seq_names = os.path.join(config["outdir"], "report", "global_snipit_labels.txt"),
	
    params: 
        ref_name = config["reference_name"]
    run:
        ref_done = "no"
       
        record_names = []
        fasta_sequences = SeqIO.parse(open(input.alignment), 'fasta')
        for record in fasta_sequences:
            record_names.append(record.name)
        with open(output.seq_names, "w") as handle:
            handle.write(f"name,label\n")
            if params.ref_name in record_names: 
                handle.write(f"{params.ref_name},Reference\n")
                ref_name = params.ref_name
                ref_done = "yes"
            for i in record_names: 
                 if i not in params.ref_name and ref_done == "no":
                     handle.write(f"{i},Reference\n")
                     ref_done = "yes"
                 elif i not in params.ref_name:
                     handle.write(f"{i},{i}\n")
                 else:
                     pass

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
	REF_NAME=$(grep "Reference" {input.names} | cut -d, -f1)
        snipit {input.aln:q} -o {params.out_path} -l {input.names} -r $REF_NAME
        """

