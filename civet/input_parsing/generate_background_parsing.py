import collections
import csv
from Bio import SeqIO

import os
import sys
from civet.utils.log_colours import green,cyan,red
from civet.utils import misc

from civet.utils.config import *

fasta_written = {}

def check_background_data_outdir(config):
    datadir = config[KEY_BACKGROUND_DATA_OUTDIR]
    if not os.path.exists(datadir):
        try:
            os.mkdir(datadir)
        except:
            sys.stderr.write(cyan(f"Cannot make background data outdir.\n"))
            sys.exit(-1)

def get_full_path(cwd, infile):
    path_to_file = os.path.abspath(cwd)
    full_path = os.path.join(path_to_file,infile)
    return full_path

def check_file_ending(fasta):

    ending = fasta.split(".")[-1]

    if ending not in ["fa","fasta","fas"]:
        sys.stderr.write(cyan(f"Please input sequence file {fasta} in fasta format.\n"))
        sys.exit(-1)

def parse_background_generation_options(generate_civet_background_data,background_data_outdir,config):
    misc.add_arg_to_config(KEY_BACKGROUND_DATA_OUTDIR,background_data_outdir,config)
    check_background_data_outdir(config)

    misc.add_arg_to_config(KEY_GENERATE_CIVET_BACKGROUND_DATA,generate_civet_background_data,config)
    input_seq_option_dict = collections.defaultdict(list)
    input_metadata_option_dict = collections.defaultdict(list)
    for input_seq_option in config[KEY_GENERATE_CIVET_BACKGROUND_DATA]:
        try:
            option,input_file = input_seq_option.split("=")
        except:
            sys.stderr.write(cyan(f"{KEY_GENERATE_CIVET_BACKGROUND_DATA} option specified with option=seq_file.fasta pairs of sequence files. See docs for full instructions.\n"))
            sys.exit(-1)
        
        if option not in INPUT_SEQ_OPTIONS:
            if option not in INPUT_METADATA_OPTIONS:
                opt_str = ', '.join(INPUT_SEQ_OPTIONS)
                opt_str2 = ','.join(INPUT_METADATA_OPTIONS)
                sys.stderr.write(cyan(f"{KEY_GENERATE_CIVET_BACKGROUND_DATA} option given `{option}` not a valid option. Please specify one of {opt_str}, {opt_str2}.\n"))
                sys.exit(-1)
        
        if option in INPUT_SEQ_OPTIONS:
            input_seq_file = get_full_path(config[KEY_CWD], input_file)
            if not os.path.exists(input_seq_file):
                sys.stderr.write(cyan(f"{input_seq_file} not found, please check path of input sequence file.\n"))
                sys.exit(-1)
            check_file_ending(input_seq_file)
            
            input_seq_option_dict[option].append(input_seq_file)

        elif option in INPUT_METADATA_OPTIONS:
            input_metadata_file = get_full_path(config[KEY_CWD], input_file)
            input_metadata_option_dict[option] = input_metadata_file
            if not os.path.exists(input_metadata_file):
                sys.stderr.write(cyan(f"{input_metadata_file} not found, please check path of input metadata file.\n"))
                sys.exit(-1)
    
    if len(input_seq_option_dict) == 0:
        sys.stderr.write(cyan(f"Please supply at least one sequence file for background data generation.\n"))
        sys.exit(-1)
    else:
        print(green("Input sequence files provided for background data generation:"))
        for option in input_seq_option_dict:
            print(green(f"{option}"))
            input_seq_option_dict[option] = list(set(input_seq_option_dict[option]))
            files = input_seq_option_dict[option]
            for f in files:
                print(f"- {f}")

    if len(input_metadata_option_dict) != 0:
        print(green("Input metadata files provided for background data generation:"))
        for option in input_metadata_option_dict:
            print(green(f"{option}"))
            input_metadata_option_dict[option] = list(set(input_metadata_option_dict[option]))
            files = input_metadata_option_dict[option]
            for f in files:
                print(f"- {f}")

    return input_seq_option_dict, input_metadata_option_dict

def try_parse_gisaid_header(fasta):
    try:
        with open(fasta) as f:
            header = f.readline()
            seq_name,seq_id,date = header.split("|")
            vir,country,rec_id,year = seq_name.split("/")
    except:
        sys.stderr.write(cyan(f"{fasta} file header not in gisaid download format, please check option specified or explicitly specify fields and delimiters.\n"))
        sys.exit(-1)

def try_parse_auspice_header(fasta):
    try:
        with open(fasta) as f:
            header = f.readline()
            vir,country,rec_id,year = header.split("/")
    except:
        sys.stderr.write(cyan(f"{fasta} file header not in auspice download format, please check option specified or explicitly specify fields and delimiters.\n"))
        sys.exit(-1)

def parse_auspice_header(header,info):
    try:
        info["header"] = header
        vir,country,seq_id,year = header.split("/")
        country = country.replace(" ","_").replace("_","").title()
        # stripping spaces and also underscores because some cases of mismatch
        sequence_name = f"{country}/{seq_id}/{year}"
        info["sequence_name"] = sequence_name
        info["country"] = country.upper()
        return info
    except:
        return False

def parse_gisaid_header(header,info):
    info["header"] = header
    seq_name,gisaid_id,date = header.split("|")
    vir,country,seq_id,year = seq_name.split("/")
    country = country.replace(" ","_").replace("_","").title()
    sequence_name = f"{country}/{seq_id}/{year}"
    info["sequence_name"] = sequence_name
    info["gisaid_id"] = gisaid_id
    info["date"] = date
    info["country"] = country.upper()
    return info

def try_parse_fasta_header(config):
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

def parse_field_delimiter_args(primary_field_delimiter,primary_metadata_fields,secondary_fields,secondary_field_delimiter,secondary_field_location,secondary_metadata_fields,config):
    misc.add_arg_to_config(KEY_PRIMARY_FIELD_DELIMTER,primary_field_delimiter,config)
    misc.add_arg_to_config(KEY_PRIMARY_METADATA_FIELDS,primary_metadata_fields,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELDS,secondary_fields,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELD_DELIMTER,secondary_field_delimiter,config)
    misc.add_arg_to_config(KEY_SECONDARY_FIELD_LOCATION,secondary_field_location,config)
    misc.add_arg_to_config(KEY_SECONDARY_METADATA_FIELDS,secondary_metadata_fields,config)

def run_sequence_qc(seq,minlen,maxambig):
    pass_qc = True
    length = len(seq)
    num_N = str(seq).upper().count("N")
    prop_N = round((num_N)/length, 2)

    if length < minlen:
        pass_qc = False

    if prop_N > maxambig: 
        pass_qc = False
    
    return pass_qc,prop_N,length

def qc_and_write_record(seq,info,minlen,maxambig,writer,fasta_out_handle,fail_qc_count):
    sequence_name = info["sequence_name"]

    pass_qc,n_content,length = run_sequence_qc(seq,minlen,maxambig)
    if pass_qc:
        info["length"] = length
        info["n_content"] = n_content
        writer.writerow(info)
        if sequence_name not in fasta_written:
            fasta_out_handle.write(f">{sequence_name}\n{seq}\n")
            fasta_written[sequence_name] = 1
    else:
        fail_qc_count +=1
    
    return fail_qc_count

def parse_gisaid_file(fasta_file,writer,header,fasta_out_handle,fail_qc_count,minlen,maxambig):

    try_parse_gisaid_header(fasta_file)
    for record in SeqIO.parse(fasta_file,"fasta"):
        info = {}
        for i in header:
            info[i] = ""
        
        info = parse_gisaid_header(record.description,info)
        if info:
            # fails if there's extra fields in the seq header, which is the case if it's non-human or env
            fail_qc_count = qc_and_write_record(record.seq,info,minlen,maxambig,writer,fasta_out_handle,fail_qc_count)
        else:
            fail_qc_count += 1 
    return fail_qc_count


def parse_auspice_file(fasta_file,writer,header,fasta_out_handle,fail_qc_count,minlen,maxambig):
    try_parse_auspice_header(fasta_file)
    for record in SeqIO.parse(fasta_file,"fasta"):
        info = {}
        for i in header:
            info[i] = ""
        
        info = parse_auspice_header(record.description,info)
        if info:
            # fails if there's extra fields in the seq header, which is the case if it's non-human or env
            fail_qc_count = qc_and_write_record(record.seq,info,minlen,maxambig,writer,fasta_out_handle,fail_qc_count)
        else:
            fail_qc_count +=1 
    return fail_qc_count

def set_up_header_for_metadata_handle(header,input_metadata_files):
    header = set(header)
    if "auspice_source_tsv" in input_metadata_files:
        with open(input_metadata_files["auspice_source_tsv"][0],"r") as f:
            reader = csv.DictReader(f,delimiter="\t")
            for field in reader.fieldnames:
                header.add(field)
    
    if "csv" in input_metadata_files:
        for csv_file in input_metadata_files["csv"]:
            with open(csv_file,"r") as f:
                reader = csv.DictReader(f)
                for field in reader.fieldnames:
                    header.add(field)
    return list(header)

def set_up_header_for_handle(input_seq_files,config):
    header_fields = set()
    header_fields.add("header")
    if "gisaid" in input_seq_files:
        for i in GISAID_HEADER_METADATA_FIELDS:
            header_fields.add(i)
    
    if "auspice_source_fasta":
        for i in AUSPICE_SOURCE_HEADER_METADATA_FIELDS:
            header_fields.add(i)
    
    if "fasta" in input_seq_files:
        primary_fields = config[KEY_PRIMARY_METADATA_FIELDS]
        for i in primary_fields:
            header_fields.add(i)

        if config[KEY_SECONDARY_FIELDS]:
            secondary_fields = config[KEY_SECONDARY_METADATA_FIELDS]
            for i in secondary_fields:
                header_fields.add(i)
    
    header_fields = list(header_fields)
    header_fields.append("n_content")
    header_fields.append("length")

    return header_fields

def parse_all_background_fasta(input_seq_files,fasta_outfile,seq_metadata_outfile,config):
    
    minlen = config[KEY_MIN_LENGTH]
    maxambig = config[KEY_MAX_AMBIGUITY]

    fasta_out_handle = open(fasta_outfile,"w")

    with open(seq_metadata_outfile, "w") as fw:
        header = set_up_header_for_handle(input_seq_files,config)

        writer = csv.DictWriter(fw, fieldnames=header, lineterminator="\n",delimiter=",")
        writer.writeheader()

        fail_qc_count = 0
        if "auspice_source_fasta" in input_seq_files:
            auspice_source_fasta_files = input_seq_files["auspice_source_fasta"]
            for auspice_file in auspice_source_fasta_files:
                fail_qc_count += parse_auspice_file(auspice_file,writer,header,fasta_out_handle,fail_qc_count,minlen,maxambig)

        if "fasta" in input_seq_files:
            pass

        if "gisaid" in input_seq_files:
            gisaid_files = input_seq_files["gisaid"]
            for gisaid_file in gisaid_files:
                fail_qc_count += parse_gisaid_file(gisaid_file,writer,header,fasta_out_handle,fail_qc_count,minlen,maxambig)

    print(f"{fail_qc_count} sequences failed QC and have been omitted.")
    fasta_out_handle.close()
    fasta_written = {}

def parse_all_background_metadata(seq_metadata_outfile,input_metadata_files,config):
    record_occurance = {}

    with open(seq_metadata_outfile,"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames

        with open(background_metadata_outfile, "w") as fw:
            header = set_up_header_for_handle(header,input_metadata_files)
            writer = csv.DictWriter(fw, fieldnames=header, lineterminator="\n")

            #first pass store occurances
            # second pass store rows if occurance >1, but stream write out if not
            # and then then write out merge












