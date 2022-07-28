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
    expanded_path = os.path.expanduser(config["cwd"])
    datadir = os.path.join(expanded_path, config[KEY_BACKGROUND_DATA_OUTDIR])
    config[KEY_BACKGROUND_DATA_OUTDIR] = datadir
    if not os.path.exists(datadir):
        try:
            os.mkdir(datadir)
        except:
            sys.stderr.write(cyan(f"Cannot make background data outdir.\n"))
            sys.exit(-1)

def parse_from_header(background_data_sequences,config):

    with open(background_data_sequences) as f: 
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


def check_bd_args(generate_background_data,background_data_sequences,background_data_metadata,config):
    
    misc.add_arg_to_config(KEY_GENERATE_BACKGROUND_DATA,generate_background_data,config)
    
    if not generate_background_data in VALUES_GENERATE_BACKGROUND_DATA:
        sys.stderr.write(cyan(f"{generate_background_data} not a valid option. Please specify either parse_seq_headers or align_curate\n"))
        sys.exit(-1)

    if generate_background_data == VALUES_GENERATE_BACKGROUND_DATA[0]:
        
        if not background_data_sequences:
            sys.stderr.write(cyan(f"Please supply a sequence file with --bd-seqs\n"))
            sys.exit(-1)
    elif generate_background_data == VALUES_GENERATE_BACKGROUND_DATA[1]:
        if not background_data_sequences:
            sys.stderr.write(cyan(f"Please supply a sequence and metadata file with --bd-seqs and --bd-metadata\n"))
            sys.exit(-1)
        if not background_data_metadata:
            sys.stderr.write(cyan(f"Please supply a sequence and metadata file with --bd-seqs and --bd-metadata\n"))
            sys.exit(-1)

        expanded_path = os.path.expanduser(config["cwd"])
        metadata = os.path.join(expanded_path, background_data_metadata)
    
        config[KEY_BACKGROUND_DATA_METADATA] = metadata

    
    expanded_path = os.path.expanduser(config["cwd"])
    unaligned_sequences = os.path.join(expanded_path, background_data_sequences)
    
    config[KEY_UNALIGNED_SEQUENCES] = unaligned_sequences
    ending = background_data_sequences.split(".")[-1]

    if ending not in ["fa","fasta","fas"]:
        sys.stderr.write(cyan(f"Please input unaligned sequences in fasta format, with file extension reflecting that.\n"))
        sys.exit(-1)


def check_seqs_metadata_match(background_sequences,background_metadata,sequence_id_column,outdir,config):
    record_dict = SeqIO.index(background_sequences, "fasta")
    
    match_dict = {}
    with open(background_metadata,"r") as f:
        reader = misc.read_csv_or_tsv(background_metadata, f)
        if not sequence_id_column in reader.fieldnames:
            sys.stderr.write(cyan(f"Please specify a sequence id column to match in the metadata with --sicol/--sequence-id-column\n"))
            sys.exit(-1)
        for row in reader:
            match = ""
            if row[sequence_id_column] not in record_dict:
                alt_matches = [row[sequence_id_column].replace(" ",""), row[sequence_id_column].replace(" ","_")]
                for alt_match in alt_matches:
                    if alt_match in record_dict:
                        match = alt_match
                        match_dict[row[sequence_id_column]] = alt_match

    if match_dict:
        config[KEY_SEQUENCE_ID_COLUMN] = f"modified_{sequence_id_column}"
        metadata_out = os.path.join(outdir, "modified_metadata.csv")
        with open(metadata_out,"w") as fw:
            with open(background_metadata,"r") as f:
                reader = misc.read_csv_or_tsv(background_metadata, f)
                header = reader.fieldnames
                header.append(f"modified_{sequence_id_column}")
                writer = csv.DictWriter(fw, fieldnames=header, lineterminator="\n")
                writer.writeheader()
                for row in reader:
                    if row[sequence_id_column] in match_dict:
                        new_row = row
                        new_row[f"modified_{sequence_id_column}"] = match_dict[row[sequence_id_column]]
                        writer.writerow(new_row)
                    else:
                        new_row = row
                        new_row[f"modified_{sequence_id_column}"] = row[sequence_id_column]
                        writer.writerow(new_row)
        return metadata_out
    else:
        return background_metadata


def sort_background_outdir(background_data_outdir,config):
    misc.add_path_to_config(KEY_BACKGROUND_DATA_OUTDIR,background_data_outdir,config)
    check_background_data_outdir(config)

def parse_metadata_from_seq_headers(primary_field_delimiter,primary_metadata_fields,secondary_fields,secondary_field_delimiter,secondary_field_location,secondary_metadata_fields,config):

    misc.add_arg_to_config(KEY_PRIMARY_FIELD_DELIMTER,primary_field_delimiter,config)
    misc.add_arg_to_config(KEY_PRIMARY_METADATA_FIELDS,primary_metadata_fields,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELDS,secondary_fields,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELD_DELIMTER,secondary_field_delimiter,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELD_LOCATION,secondary_field_location,config)
    misc.add_arg_to_config(KEY_SECONDARY_METADATA_FIELDS,secondary_metadata_fields,config)

    parse_from_header(config[KEY_BACKGROUND_SEQUENCES],config)



