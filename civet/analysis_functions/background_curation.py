#!/usr/bin/env python

import csv
from Bio import SeqIO
import os

import sys
from civet.utils.log_colours import green,cyan,red
from civet.utils.config import *

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
  