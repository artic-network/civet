#!/usr/bin/env python3
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import sys
import os
import csv
import pkg_resources
from collections import Counter
from Bio import SeqIO

def check_dir_for_file(file_list, extension,argument,key,config,required):
    match = []
    for fn in file_list:
        if fn.endswith(extension):
            match.append(fn)
    if len(match) > 1:
        sys.stderr.write(cyan(f"Error: More than 1 {extension} file found in {config['datadir']},\nplease specify which {extension} file using {argument}.\n"))
        sys.exit(-1)
    elif len(match) == 0:
        if required==True:
            sys.stderr.write(cyan(f"Error: No {extension} file found in {config['datadir']},\nplease specify which {extension} file using {argument}.\n"))
            sys.exit(-1)
    elif len(match) == 1:
        config[key] = match[0]


def background_data_load(config):
    file_list = []
    for r,d,f in os.walk(config["datadir"]):
        for fn in f:
            file_list.append(os.path.join(r,fn))
    if not config["background_metadata"]:
        check_dir_for_file(file_list, "metadata.csv","-bm/--background-metadata","background_metadata",config,True)
    if not config["background_sequences"]:
        check_dir_for_file(file_list, "fasta","-bseq/--background-sequences","background_sequences",config,True)

    if not config["background_snps"]:
        check_dir_for_file(file_list, "mutations.csv","-bsnp/--background-snps","background_snps",config,False)

    if not config["background_tree"]:
        check_dir_for_file(file_list, "newick","-bt/--background-tree","background_tree",config,False)

def check_datadir(config):
    datadir = config["datadir"]
    if not os.path.exists(datadir):
        if not config["background_metadata"] or not config["background_sequences"]:
            print(cyan(f"Error: data directory not found: {datadir}.\n")+"Please check the datadir path (-d/ --datadir) provided exists or supply files with -bm/--background-metadata and -bseq/--background-sequences.\n")
            sys.exit(-1)

def check_csv_file(argument,description,csv_file,background_column,fasta_column):
    
    ending = csv_file.split(".")[-1]

    if ending in ["xls","xlsx"]:
        sys.stderr.write(cyan(f"Error: it looks like you've provided an excel file as input.\n{argument} accepts a csv file.\n"))
        sys.exit(-1)
    elif ending == "csv":
        pass
    else:
        sys.stderr.write(cyan(f"Error: please provide a csv file for {argument}.\n"))
        sys.exit(-1)
    
    if os.path.isfile(csv_file):
        try:
            c = 0
            with open(csv_file,"r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    c +=1
            print(green(f"{c} rows in {description} file."))
        except:
            sys.stderr.write(cyan(f"Unable to read {description} file, please check your {description} is in the correct format with a header.\n"))
            sys.exit(-1)
    else:
        sys.stderr.write(cyan(f"Cannot find {description} file: ")+f"{csv_file}\n")
        sys.exit(-1)
    print(green(f"{description.replace('b','B')} file:") + f" {csv_file}")

    with open(csv_file,"r") as f:
        reader = csv.DictReader(f)
        if background_column in reader.fieldnames and fasta_column in reader.fieldnames:
            pass
        else:
            encode = False
            for i in reader.fieldnames:
                if "ufeff" in i:
                    sys.stderr.write(cyan(f"Error: it appears your {description} file may have been edited in Excel and now contains hidden characters.\n") + "Please remove said characters in a text editor and try again.")
                    sys.exit(-1)
            else:
                if background_column not in reader.fieldnames:
                    sys.stderr.write(cyan(f"Error: {background_column} column not found in {description} file. Please specifiy which column to match with `-bicol/--background-id-column.`\n"))
                    sys.exit(-1)
                elif fasta_column not in reader.fieldnames:
                    sys.stderr.write(cyan(f"Error: {fasta_column} column not found in {description} file. Please specifiy which column to match with `-sicol/--sequence-id-column.`\n"))
                    sys.exit(-1)
    return c

def check_background_fasta(background_fasta):
    ending = background_fasta.split(".")[-1]

    if ending not in ["fa","fasta","fas"]:
        sys.stderr.write(cyan(f"Please input background sequences in fasta format, with file extension reflecting that.\n"))
        sys.exit(-1)
    else:
        print(green(f"Background fasta file:") + f" {background_fasta}")
    lens = Counter()
    c = 0
    for record in SeqIO.parse(background_fasta, "fasta"):
        lens[len(record)] += 1
        c +=1

    if c == 0:
        sys.stderr.write(cyan(f"Error: no records found in fasta file, please check your background fasta is in the correct format.\n"))
        sys.exit(-1)
    else:
        print(green(f"{c} records in background fasta file."))
    
    if len(lens)>1:
        len_string = ""
        for i in lens:
            if lens[i] == 1:
                len_string += f"\t- {lens[i]} record is {i} bases\n"
            else:
                len_string += f"\t- {lens[i]} records are {i} bases\n"
        sys.stderr.write(cyan(f"Error: not all records are the same length, please check your background fasta is aligned correctly.\n") + len_string)
        sys.exit(-1)
    return c

def check_background_snps(config):
    records = 0
    if config["background_snps"]:
        fields = ['query',"SNPs","ambiguities","SNPcount","ambcount"]
        missing = []
        with open(config["background_snps"], "r") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            for field in fields:
                if field not in header:
                    missing.append(field)
            
                if len(missing) > 1:
                    sys.stderr.write(cyan(f"Error: some required fields missing from background SNPs file:\n") + "\n - ".join(missing) + "\n")
                    sys.exit(-1)
                elif len(missing) == 1:
                    sys.stderr.write(cyan(f"Error: background SNPs file missing field: ") + missing[0] + "\n")
                    sys.exit(-1)

                for row in reader:
                    records+=1
    
    return records
        


def data_group_parsing(debug,datadir,background_csv,background_SNPs,background_fasta,background_tree,background_column,fasta_column,config):
    """
    parses the data group arguments 
    --datadir (Default $DATADIR)
    --background-csv (Default False)
    --background-fasta (Default False)
    --background-id-column (Default central_sample_id)
    """
    if datadir:
        
        datadir = os.path.abspath(datadir)

    # if command line arg, overwrite config value
    misc.add_arg_to_config("datadir",datadir,config)
    misc.add_file_to_config("background_metadata",background_csv,config)
    misc.add_file_to_config("background_snps",background_SNPs,config)
    misc.add_file_to_config("background_sequences",background_fasta,config)
    misc.add_file_to_config("background_tree",background_tree,config)
    misc.add_arg_to_config("background_id_column",background_column,config)
    misc.add_arg_to_config("sequence_id_column",fasta_column,config)

    # needs either datadir specified or both the files specified
    if not config["datadir"]:
        if not config["background_metadata"] and not config["background_sequences"]:
            sys.stderr.write(cyan(f"Error: insufficient background data supplied, please provide a background csv and fasta file.\n"))
            sys.exit(-1)
        
    if config["datadir"]:
        # does the path exist? if not and if neither file specified exit
        check_datadir(config)
    
    # finds the files in the datadir, exits if unclear which files are which
    background_data_load(config)
    
    # check if it's the right file type, 
    # see if you can read it, 
    # check if background_column in csv header, 
    # count the number of records

    if not config["sequence_id_column"]:
        config["sequence_id_column"] = config["background_id_column"]
    
    if config["background_snps"]:
        config["background_search_file"] = config["background_snps"]
        check_background_snps(config)
    else:
        config["background_search_file"] = config["background_sequences"]

    if debug:
        csv_record_count = check_csv_file("-bm/--background-metadata","background csv",config["background_metadata"],config["background_id_column"],config["sequence_id_column"])
        fasta_record_count = check_background_fasta(config["background_sequences"])
        
        if csv_record_count != fasta_record_count:
            sys.stderr.write(cyan(f"Error: different number of background csv and background fasta records.\n")+"Please provide a fasta record for each row in the background metadata file.\n")
            sys.exit(-1)

        if config["background_snps"]:
            SNP_record_count = check_background_snps(config)
            if csv_record_count != SNP_record_count:
                sys.stderr.write(cyan(f"Error: different number of background csv and background SNP records.\n")+"Please provide a SNP record for each row in the background metadata file.\n")
                sys.exit(-1)
