#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
from civet.utils.config import *
"""
config = {
    KEY_DATA_OUTDIR:"civet_data",
    KEY_PRIMARY_FIELD_DELIMTER:"|",
    KEY_SECONDARY_FIELD_DELIMTER:"/",
    KEY_SECONDARY_FIELD_LOCATION:0,
    KEY_PRIMARY_METADATA_FIELDS:"sequence_name,gisaid_id,sample_date",
    KEY_SECONDARY_METADATA_FIELDS:"virus,country,sequence_id,year"
    # KEY_UNALIGNED_SEQUENCES:
    }
"""

def input_fasta_qc(input_fasta,output_fasta,output_notes,config):
    """
    Checking input fasta file for:
    - Minimum sequence length
    - Maximum ambiguities
    - Whether the records are duplicated
    """
    minlen = config[KEY_MIN_LENGTH]
    maxambig = config[KEY_MAX_AMBIGUITY]
    passed_qc = 0
    with open(output_notes,"w") as fw2:
        fw2.write("sequence_header,N_count,proportion_N,seq_length,QC_status\n")
        with open(output_fasta,"w") as fw:
            for record in SeqIO.parse(input_fasta,"fasta"):
                failed_qc = False
                length = len(record)
                num_N = str(record.seq).upper().count("N")
                prop_N = round((num_N)/length, 2)

                if length < minlen:
                    failed_qc = True

                if prop_N > maxambig: 
                    failed_qc = True
                
                if failed_qc:
                    fw2.write(f"{record.description},{num_N},{prop_N},{length},failed_qc\n")
                else:
                    fw2.write(f"{record.description},{num_N},{prop_N},{length},passed_qc\n")
                    fw.write(f">{record.description}\n{record.seq}\n")
                    passed_qc +=1

    if passed_qc == 0:
        sys.stderr.write(cyan(f"No supplied sequences pass the QC steps.\n") + f"""Please ensure sequences meet the following:
\t- Minimum sequence length (>{minlen} bases)
\t- Maximum ambiguities (<{maxambig} proportion N)
\t- Whether the records match a query supplied in an ID string or input csv
\t- Whether the records are duplicated in the file\n""" + cyan("You can change the default QC settings with `-n/--max-ambiguity` and `-l/--min-length`."))
        sys.exit(-1)
    
if config[KEY_BACKGROUND_DATA_ALIGN_ONLY]:
    rule all:
        input:
            os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_sequences.aln.fasta"),
            os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_snps.csv")
else:
    rule all:
        input:
            os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_sequences.aln.fasta"),
            os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_snps.csv"),
            os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_metadata.csv")
                    

rule qc_input_fasta:
    input:
        fasta = config[KEY_UNALIGNED_SEQUENCES]
    output:
        fasta = os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"qc_step","background_sequences.passed_qc.fasta"),
        notes = os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"qc_step","qc_notes.txt")
    run:
        input_fasta_qc(input.fasta,output.fasta,output.notes,config)

rule align_to_reference:
    input:
        fasta = rules.qc_input_fasta.output.fasta,
        reference = config[KEY_REFERENCE_SEQUENCE]
    params:
        trim_start = config[KEY_TRIM_START],
        trim_end = config[KEY_TRIM_END],
        sam = os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"mapped.sam")
    output:
        fasta = os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_sequences.aln.fasta")
    log:
        os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR], "logs","minimap2_sam.log")
    shell:
        """
        minimap2 -a -x asm20 --sam-hit-only --secondary=no --score-N=0  -t  {workflow.cores} {input.reference:q} '{input.fasta}' -o {params.sam:q} &> {log:q} 
        gofasta sam toMultiAlign \
            -s {params.sam:q} \
            -t {workflow.cores} \
            --reference {input.reference:q} \
            --trimstart {params.trim_start} \
            --trimend {params.trim_end} \
            --trim \
            --pad > '{output.fasta}'
        """

rule gofasta_SNPs:
    input:
        fasta = rules.align_to_reference.output.fasta,
        reference = config[KEY_REFERENCE_SEQUENCE]
    output:
        csv = os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_snps.csv")
    shell:
        """
        gofasta updown list \
            -q {input.fasta:q} \
            -o {output.csv:q} \
            --reference {input.reference:q} 
        """

rule generate_metadata:
    input:
        fasta = config[KEY_UNALIGNED_SEQUENCES]
    output:
        csv = os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_metadata.csv")
    run:
        with open(output.csv,"w") as fw:
            if config[KEY_SECONDARY_FIELDS]:
                header_string = "sequence_header,"+ config[KEY_PRIMARY_METADATA_FIELDS] + ',' + config[KEY_SECONDARY_METADATA_FIELDS]
            else:
                header_string = "sequence_header,"+config[KEY_PRIMARY_METADATA_FIELDS]

            header = header_string.split(",")

            writer = csv.DictWriter(fw, fieldnames=header,lineterminator='\n')
            writer.writeheader()

            for record in SeqIO.parse(input.fasta,"fasta"):
                primary_fields = record.description.split(config[KEY_PRIMARY_FIELD_DELIMTER])

                row = {"sequence_header":record.description}

                for col_name,field in zip(config[KEY_PRIMARY_METADATA_FIELDS].split(","),primary_fields):
                    row[col_name]=field
                
                if config[KEY_SECONDARY_FIELDS]:
                    location = config[KEY_SECONDARY_FIELD_LOCATION]
                    secondary_info = primary_fields[location]
                    secondary_fields = secondary_info.split(config[KEY_SECONDARY_FIELD_DELIMTER])
                    for col_name,field in zip(config[KEY_SECONDARY_METADATA_FIELDS].split(','),secondary_fields):
                        
                        row[col_name]=field

                writer.writerow(row)
