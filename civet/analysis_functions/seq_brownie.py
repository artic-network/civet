
from civet.utils.log_colours import green,cyan,red
from civet.analysis_functions import catchment_parsing
from civet.utils import misc
from civet.report_functions import report
import collections
import sys
import yaml
from Bio import SeqIO
import csv

from civet.utils.config import *


def seq_brownie(query_fasta,output_fasta,output_csv,config):
    print(config[KEY_MATCHED_FASTA])
    print("trying to run code")
    records = 0
    
    seq_map = {}
    hash_map = collections.defaultdict(list)
    hash_map_for_metadata = {}

    if config[KEY_MATCHED_FASTA]:
        records = catchment_parsing.add_to_hash(config[KEY_MATCHED_FASTA],seq_map,hash_map,records)

    if config[KEY_QUERY_FASTA]:
        records = catchment_parsing.add_to_hash(query_fasta,seq_map,hash_map,records)
    
    with open(output_fasta,"w") as fseqs:
        for key in seq_map:
            fseqs.write(f">{key}\n{seq_map[key]}\n")

    for hash_str in hash_map:
        for record_id in hash_map[hash_str]:
            hash_map_for_metadata[record_id] = hash_str
            
    if config[KEY_QUERY_FASTA]:
        misc.add_col_to_metadata(KEY_HASH, hash_map_for_metadata, config[KEY_QUERY_METADATA], output_csv, config["input_id_column"], config)
    elif config[KEY_MATCHED_FASTA]:
        misc.add_col_to_metadata(KEY_HASH, hash_map_for_metadata, config[KEY_QUERY_METADATA], output_csv, config["sequence_id_column"], config)

    config[KEY_QUERY_METADATA] = output_csv
    print(green("Query sequences collapsed from ") + f"{records}" +green(" to ") + f"{len(seq_map)}" + green(" unique sequences."))
            