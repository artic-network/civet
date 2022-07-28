#!/usr/bin/env python

import csv
from Bio import SeqIO
import os

import sys
from civet.utils.log_colours import green,cyan,red
import civet.analysis_functions.background_curation as bc
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
        fasta = os.path.join(config[KEY_BACKGROUND_DATA_TEMPDIR],"qc_step","background_sequences.passed_qc.fasta"),
        notes = os.path.join(config[KEY_BACKGROUND_DATA_TEMPDIR],"qc_step","qc_notes.txt")
    run:
        bc.input_fasta_qc(input.fasta,output.fasta,output.notes,config)

rule align_to_reference:
    input:
        fasta = rules.qc_input_fasta.output.fasta,
        reference = config[KEY_REFERENCE_SEQUENCE]
    params:
        trim_start = config[KEY_TRIM_START],
        trim_end = config[KEY_TRIM_END],
        sam = os.path.join(config[KEY_BACKGROUND_DATA_TEMPDIR],"mapped.sam")
    output:
        fasta = os.path.join(config[KEY_BACKGROUND_DATA_OUTDIR],"background_sequences.aln.fasta")
    log:
        os.path.join(config[KEY_BACKGROUND_DATA_TEMPDIR], "logs","minimap2_sam.log")
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
