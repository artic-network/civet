import collections
import csv
from Bio import SeqIO

import os
import sys
from civet.utils.log_colours import green,cyan,red
from civet.utils import misc

from civet.utils.config import *

"""
    dc_group.add_argument("-bd","--generate-civet-background-data",dest="generate_civet_background_data",action="store",help="Create background metadata, sequence and SNP file from GISAID download")
    dc_group.add_argument("--background-data-checks",dest="debug",action="store_true",help="Run checks on custom background data files, not run by default")
    dc_group.add_argument("--unaligned-sequences",dest="unaligned_sequences",action="store",help="Sequence header fields to create metadata file from. Default: `sequence_name,gisaid_id,sample_date`")
    dc_group.add_argument("--data-outdir",dest="data_outdir",action="store",help="Directory to output the civet background data. Default: `civet_data`")
    dc_group.add_argument("--primary-field-delimiter",dest="primary_field_delimiter",action="store",help="Primary sequence header field delimiter to create metadata file from. Default: `|`")
    dc_group.add_argument("--primary-metadata-fields",dest="primary_metadata_fields",action="store",help="Primary sequence header fields to create metadata file from. Default: `sequence_name,gisaid_id,sample_date`")
    dc_group.add_argument("--secondary-field-delimiter",dest="secondary_field_delimiter",action="store",help="Secondary sequence header field delimiter to create metadata file from. Default: `/`")
    dc_group.add_argument("--secondary-field-location",dest="secondary_field_location",action="store",help="Secondary sequence header location within primary field list. Default: `0` (i.e. the first field)")
    dc_group.add_argument("--secondary-metadata-fields",dest="secondary_metadata_fields",action="store",help="Secondary sequence header fields to create metadata file from. Default: `virus,country,sequence_id,year`")
"""
def check_background_data_outdir(config):
    datadir = config[KEY_BACKGROUND_DATA_OUTDIR]
    if not os.path.exists(datadir):
        try:
            os.mkdir(datadir)
        except:
            sys.stderr.write(cyan(f"Cannot make background data outdir.\n"))
            sys.exit(-1)

def check_input_sequences(config):
    background_fasta = config["generate_civet_background_data"]
    config[KEY_UNALIGNED_SEQUENCES] = background_fasta
    ending = background_fasta.split(".")[-1]

    if ending not in ["fa","fasta","fas"]:
        sys.stderr.write(cyan(f"Please input unaligned sequences in fasta format, with file extension reflecting that.\n"))
        sys.exit(-1)
    
    with open(background_fasta) as f: 
        header = f.readline()

        primary_fields = header.split(config[KEY_PRIMARY_FIELD_DELIMTER])

        if len(primary_fields) != len(config[KEY_PRIMARY_METADATA_FIELDS].split(",")):
            sys.stderr.write(cyan(f"Number of primary fields in sequence header does not match number of primary field names. Check delimiter and fields specified.\n")+f"Primary field names: {config[KEY_PRIMARY_METADATA_FIELDS]}\nPrimary field delimiter: {config[KEY_PRIMARY_FIELD_DELIMTER]}\n")
            sys.exit(-1)

        if config[KEY_SECONDARY_FIELDS]:
            try:
                secondary_string = primary_fields[config[KEY_SECONDARY_FIELD_LOCATION]]
            except:
                sys.stderr.write(cyan(f"Secondary field location not found in primary fields.\n"))
                sys.exit(-1)

            secondary_fields = secondary_string.split(config[KEY_SECONDARY_FIELD_DELIMTER])
            if len(secondary_fields) != len(config[KEY_SECONDARY_METADATA_FIELDS].split(",")):
                sys.stderr.write(cyan(f"Number of secondary fields in sequence header does not match number of secondary field names. Check delimiter and fields specified.\n")+f"Secondary field names: {config[KEY_SECONDARY_METADATA_FIELDS]}\nSecondary field delimiter: {config[KEY_SECONDARY_FIELD_DELIMTER]}\n")
                sys.exit(-1)


def parse_generate_background_args(generate_civet_background_data,background_data_outdir,primary_field_delimiter,primary_metadata_fields,secondary_fields,secondary_field_delimiter,secondary_field_location,secondary_metadata_fields,config):
    misc.add_arg_to_config("generate_civet_background_data",generate_civet_background_data,config)
    misc.add_path_to_config(KEY_BACKGROUND_DATA_OUTDIR,background_data_outdir,config)
    misc.add_arg_to_config(KEY_PRIMARY_FIELD_DELIMTER,primary_field_delimiter,config)
    misc.add_arg_to_config(KEY_PRIMARY_METADATA_FIELDS,primary_metadata_fields,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELDS,secondary_fields,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELD_DELIMTER,secondary_field_delimiter,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELD_LOCATION,secondary_field_location,config)
    misc.add_arg_to_config(KEY_SECONDARY_METADATA_FIELDS,secondary_metadata_fields,config)

    check_input_sequences(config)

    check_background_data_outdir(config)



