
from Bio import Phylo
from Bio import SeqIO
import csv
import collections

rule all:
    input:
        os.path.join(config["tempdir"], "report", "global_snipit_labels.txt"),
	os.path.join(config["outdir"],"report","figures","genome_graph_global_focal.png"),
        lambda wildcards: os.path.join(config["tempdir"],"all_focal.reference_mapped.sam") if not config["focal_alignment"] else [],
        lambda wildcards: os.path.join(config["tempdir"],"all_focal.reference_mapped.fasta") if not config["focal_alignment"] else [],
        lambda wildcards: os.path.join(config["tempdir"],"all_focal_aln.fasta") if not config["focal_alignment"] else []

rule minimap_all_focal: 
    input:
        focal = config["fasta"],
        ref = config["reference_fasta"]
    output: 
        sam = os.path.join(config["tempdir"],"all_focal.reference_mapped.sam")
    shell: 
        """
        minimap2 -a -x asm5 {input.ref:q} {input.focal:q} -o {output.sam:q}
        """

rule focal_sam_to_fasta: 
    input: 
        sam = rules.minimap_all_focal.output.sam,
        reference = config["reference_fasta"]
    output: 
        fasta = os.path.join(config["tempdir"],"all_focal.reference_mapped.fasta")
    shell: 
        """
        datafunk sam_2_fasta \
        -s {input.sam:q} \
        -r {input.reference:q} \
        -o {output.fasta:q} \
        --pad \
        --log-inserts
        """

rule gather_focal_seqs: 
    input: 
        outgroup_fasta = config["outgroup_fasta"],
        focal_align = rules.focal_sam_to_fasta.output.fasta
    output: 
        full_aln = os.path.join(config["tempdir"],"all_focal_aln.fasta")
    run: 
        with open(output.full_aln, "w") as handle: 
            for record in SeqIO.parse(input.outgroup_fasta, "fasta"):
                handle.write(f">Reference\n{record.seq}\n")
            for record in SeqIO.parse(input.focal_align, "fasta"):
                handle.write(f">{record.description}\n{record.seq}\n")


if not config["focal_alignment"]:
    ALIGN = rules.gather_focal_seqs.output.full_aln
else: 
    ALIGN = config["focal_alignment"]

rule get_all_sequence_names:
    input:
        alignment = ALIGN
    output:
        seq_names = os.path.join(config["tempdir"], "report", "global_snipit_labels.txt")
	
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
        aln = ALIGN,
	names = rules.get_all_sequence_names.output.seq_names
    params:
        out_path = os.path.join(config["outdir"],"report","figures","genome_graph_global_focal")
    output:
        os.path.join(config["outdir"],"report","figures","genome_graph_global_focal.png")
    shell:
        """
        snipit {input.aln:q} -o {params.out_path} -l {input.names} -r $(grep "Reference" {input.names} | head -n 1 | cut -d, -f1) --sort-by-mutation-number --high-to-low
        """

