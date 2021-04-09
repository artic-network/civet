#!/usr/bin/env python3
from Bio import SeqIO
import log_colours as colour
import sys


def input_fasta_qc(fasta,query_ids,minlen,maxambig):
    """
    Checking input fasta file for:
    - Minimum sequence length
    - Maximum ambiguities
    - Whether the records match a query id
    - Whether the records are duplicated
    """
    failed_qc = {}
    passed_qc = []

    for record in SeqIO.parse(fasta, "fasta"):
        if len(record) < minlen:
            failed_qc[record.id] = f"Sequence too short: Sequence length {len(record)}"
        else:
            num_N = str(record.seq).upper().count("N")
            prop_N = round((num_N)/len(record.seq), 2)
            if prop_N > maxambig: 
                failed_qc[record.id] = f"N content: {prop_N}"
            else:
                if record.id not in query_ids:
                    failed_qc[record.id] = "Sequence name doesn't match any query under investigation."
                else:
                    if record.id in passed_qc:
                        failed_qc[record.id] = "Sequence duplicated in input fasta file."
                    else:
                        passed_qc.append(record)
    
    print(colour.green(f"{len(passed_qc)} query sequences have passed QC."))
    print(colour.cyan(f"{len(failed_qc)} query sequences have failed QC."))

    if len(passed_qc) == 0:
        sys.stderr.write(cyan(f"No supplied sequences pass the QC steps.\n") + f"""Please ensure sequences meet the following:
\t- Minimum sequence length (>{minlen} bases)
\t- Maximum ambiguities (<{} proportion N)
\t- Whether the records match a query id
\t- Whether the records are duplicated in the file""")
        sys.exit(-1)
    return passed_qc, failed_qc

query_ids = []