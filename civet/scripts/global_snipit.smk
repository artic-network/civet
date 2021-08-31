from civet.utils.log_colours import green,cyan,red
from civet.analysis_functions import catchment_parsing
from civet.utils import misc
from civet.report_functions import report
from Bio import Phylo
from Bio import SeqIO
import csv
import collections
from civet.utils.config import *

rule all:
    input:
        os.path.join(config["tempdir"], "snipit", "global_snipit_labels.txt"),
	os.path.join(config["tempdir"],"snipit", "global_snipit.svg")

rule gather_focal_seqs: 
    input: 
        outgroup_fasta = config[KEY_REFERENCE_SEQUENCE],
        focal_align = config["fasta"]
    output: 
        full_aln = os.path.join(config["tempdir"],"all_focal_aln.fasta")
    run: 
        with open(output.full_aln, "w") as handle: 
            for record in SeqIO.parse(input.outgroup_fasta, "fasta"):
                handle.write(f">Reference\n{record.seq}\n")
            for record in SeqIO.parse(input.focal_align, "fasta"):
                handle.write(f">{record.description}\n{record.seq}\n")


if not config["focal_alignment"]:
    ALIGN = config["fasta"]
else: 
    ALIGN = config["focal_alignment"]

rule get_all_sequence_names:
    input:
        alignment = ALIGN
    output:
        seq_names = os.path.join(config["tempdir"], "snipit", "global_snipit_labels.txt")
	
    params: 
        ref_name = config[KEY_REFERENCE_SEQUENCE]
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
        aln = ALIGN,
	names = rules.get_all_sequence_names.output.seq_names
    params:
        out_path = os.path.join(config["tempdir"],"snipit", "global_snipit")
    output:
        os.path.join(config["tempdir"],"snipit", "global_snipit.svg")
    shell:
        """
        snipit {input.aln:q} -o {params.out_path} -l {input.names} -r $(grep "Reference" {input.names} | head -n 1 | cut -d, -f1) -f svg
        """

