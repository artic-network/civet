#!/usr/bin/env python3

import os
import argparse
import csv
import datetime 
import sys
from Bio import SeqIO
from datetime import datetime 
from tempfile import gettempdir
import tempfile
import pkg_resources

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'

def get_colours():
    colours = ['viridis', 'plasma', 'inferno', 'magma', 'cividis','Greys', 
            'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
            'twilight', 'twilight_shifted', 'hsv',
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c',
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']
    return colours

def colour(text, text_colour):
    bold_text = 'bold' in text_colour
    text_colour = text_colour.replace('bold', '')
    underline_text = 'underline' in text_colour
    text_colour = text_colour.replace('underline', '')
    text_colour = text_colour.replace('_', '')
    text_colour = text_colour.replace(' ', '')
    text_colour = text_colour.lower()
    if 'red' in text_colour:
        coloured_text = RED
    elif 'green' in text_colour:
        coloured_text = GREEN
    elif 'yellow' in text_colour:
        coloured_text = YELLOW
    elif 'dim' in text_colour:
        coloured_text = DIM
    elif 'cyan' in text_colour:
        coloured_text = 'cyan'
    else:
        coloured_text = ''
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text

def red(text):
    return RED + text + END_FORMATTING

def cyan(text):
    return CYAN + text + END_FORMATTING

def green(text):
    return GREEN + text + END_FORMATTING

def yellow(text):
    return YELLOW + text + END_FORMATTING

def bold_underline(text):
    return BOLD + UNDERLINE + text + END_FORMATTING

def get_snakefile(thisdir):
    snakefile = os.path.join(thisdir, 'scripts','Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}\n Check installation'))
        sys.exit(-1)
    return snakefile

# def get_seqs_for_aln(seqs_arg,cwd):
#     if not seqs_arg:
#         sys.stderr.write(cyan(f"""Error: please input fasta file for alignment\n"""))
#         sys.exit(-1)
#     else:
#         seqs = os.path.join(cwd, seqs_arg)

#     if not os.path.exists(seqs):
#         sys.stderr.write(cyan(f"""Error: cannot find sequence file at {seqs}\n"""))
#         sys.exit(-1)
#     return seqs

def get_outgroup_sequence(outgroup_arg, cwd, config):
    if outgroup_arg:
        reference_fasta = os.path.join(cwd, outgroup_arg)
        if not os.path.isfile(reference_fasta):
            sys.stderr.write(cyan(f"""Error: cannot find specified outgroup file at {outgroup_arg}\n"""))
            sys.exit(-1)
        else:
            config["reference_fasta"] = reference_fasta
    else:
        reference_fasta = pkg_resources.resource_filename('llama', 'data/reference.fasta')
        config["reference_fasta"] = reference_fasta

def get_query_fasta(fasta_arg,cwd,config):
    if fasta_arg:
        fasta = os.path.join(cwd, fasta_arg)
        if not os.path.exists(fasta):
            sys.stderr.write(cyan(f'Error: cannot find fasta query at {fasta}\n'))
            sys.exit(-1)
        else:
            print(green(f"Input fasta file:") + f" {fasta}")
    else:
        fasta = ""
    config["fasta"] = fasta 

def get_outdir(outdir_arg,cwd,config):
    outdir = ''
    if outdir_arg:
        rel_outdir = outdir_arg #for report weaving
        outdir = os.path.join(cwd, outdir_arg)
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
    else:
        timestamp = str(datetime.now().isoformat(timespec='milliseconds')).replace(":","").replace(".","").replace("T","-")
        outdir = os.path.join(cwd, timestamp)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        rel_outdir = os.path.join(".",timestamp)
    print(green(f"Output dir:") + f" {outdir}")
    config["outdir"] = outdir 
    config["rel_outdir"] = rel_outdir 
    return outdir 
        
def get_temp_dir(tempdir_arg,no_temp_arg, cwd,config):
    tempdir = ''
    if no_temp_arg:
        print(green(f"--no-temp:") + f" All intermediate files will be written to {outdir}")
        tempdir = outdir
    elif tempdir_arg:
        to_be_dir = os.path.join(cwd, tempdir_arg)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name
    else:
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
        tempdir = temporary_directory.name
    
    config["tempdir"] = tempdir 
    return tempdir

# def check_data_dir(datadir,no_seqs,cwd,config):
#     data_dir = os.path.join(cwd, datadir)
    
#     seqs = os.path.join(data_dir,"alignment.fasta")
    
#     metadata = os.path.join(data_dir,"metadata.csv")

#     tree = os.path.join(data_dir,"global.tree")
#     if no_seqs:
#         if not os.path.isfile(metadata) or not os.path.isfile(tree):
#             sys.stderr.write(cyan(f"""Error: cannot find correct data files at {data_dir}\nThe directory should contain the following files:\n\
#     - global.tree\n\
#     - metadata.csv\n"""))
#             sys.exit(-1)
#         else:
#             config["metadata"] = metadata
#             config["tree"] = tree
#             seqs = ""
#             print(green("Found input data files:") + f"\n - {metadata}\n - {tree}")
#     else:
#         if not os.path.isfile(seqs) or not os.path.isfile(metadata) or not os.path.isfile(tree):
#             sys.stderr.write(cyan(f"""Error: cannot find correct data files at {data_dir}\nThe directory should contain the following files:\n\
#     - alignment.fasta\n\
#     - global.tree\n\
#     - metadata.csv\n"""))
#             sys.exit(-1)
#         else:
#             config["seqs"] = seqs
#             config["metadata"] = metadata
#             config["tree"] = tree

#             print(green("Found input data files:") + f" - {seqs}\n - {metadata}\n - {tree}")
#     return metadata,seqs,tree

# def parse_from_metadata_arg(metadata, from_metadata, data_column, config):
#     queries = []
#     query_dict = {}
#     column_names =""
#     config["input_column"] = data_column
#     with open(metadata, newline="") as f:
#         reader = csv.DictReader(f)
#         column_names = reader.fieldnames
        
#         for factor in from_metadata:
#             column_name,to_search = factor.split("=")
#             if column_name in column_names:
#                 query_dict[column_name] = to_search
#             else:
#                 cols = "\n- ".join(column_names)
#                 cols = cols + "\n"
#                 sys.stderr.write(cyan(f"""Error: `--from-metadata` argument contains a column {column_name} that is not found in the metadata file supplied.
# Columns that were found:\n{cols}"""))
#                 sys.exit(-1)
    
#     rows_to_search = []
    
#     for column_name in query_dict:
#         to_search = query_dict[column_name]
#         if ':' in to_search:
#             if to_search.startswith("2020-") or to_search.startswith("2019-") or to_search.startswith("2021-"):
#                 print(f"Date range detected: {to_search}")
#                 date_range = to_search.split(":")
#                 start_date = datetime.strptime(date_range[0], "%Y-%m-%d").date()
#                 end_date = datetime.strptime(date_range[1], "%Y-%m-%d").date()
#                 if rows_to_search == []:
#                     with open(metadata, newline="") as f:
#                         reader = csv.DictReader(f)
#                         c =0
#                         for row in reader:
#                             c +=1
#                             row_date = row[column_name]
#                             try:
#                                 check_date = datetime.strptime(row_date, "%Y-%m-%d").date()
#                             except:
#                                 sys.stderr.write(cyan(f"Error: Metadata field `{row_date}` [at column: {column_name}, row: {c}] contains unaccepted date format\Please use format `2020-05-19`\n"))
#                                 sys.exit(-1)
#                             if start_date <= check_date <= end_date:
#                                 rows_to_search.append((row,c))
#                 else:
#                     last_rows_to_search = rows_to_search
#                     new_rows_to_search = []
#                     for row,c in last_rows_to_search:
#                         row_date = row[column_name]
#                         try:
#                             check_date = datetime.strptime(row_date, "%Y-%m-%d").date()
#                         except:
#                             sys.stderr.write(cyan(f"Error: Metadata field `{row_date}` [at column: {column_name}, row: {c}] contains unaccepted date format\Please use format `YYYY-MM-DD`\n"))
#                             sys.exit(-1)
#                         if start_date <= check_date <= end_date:
#                             new_rows_to_search.append((row,c))
                    
#                     rows_to_search = new_rows_to_search
#         else:
#             if rows_to_search == []:
#                 with open(metadata, newline="") as f:
#                     reader = csv.DictReader(f)
#                     c =0
#                     for row in reader:
#                         c +=1
#                         row_info = row[column_name]
                        
#                         if row_info == to_search:
#                             rows_to_search.append((row,c))
#             else:
#                 last_rows_to_search = rows_to_search
#                 new_rows_to_search = []
#                 for row,c in last_rows_to_search:
#                     row_info = row[column_name]
                    
#                     if row_info == to_search:
#                         new_rows_to_search.append((row,c))
#                 rows_to_search = new_rows_to_search
#     query = os.path.join(config["outdir"], "from_metadata_query.csv")
#     with open(query,"w") as fw:
#         writer = csv.DictWriter(fw, fieldnames=column_names,lineterminator='\n')
#         writer.writeheader()
#         count = 0
#         query_ids = []
#         for row,c in rows_to_search:
#             writer.writerow(row)
#             count +=1
#             query_ids.append(row[data_column])
#         if count == 0:
#             sys.stderr.write(cyan(f"Error: No sequences meet the criteria defined with `--from-metadata`.\nExiting\n"))
#             sys.exit(-1)
#         print(green(f"Number of sequences matching defined query:") + f" {count}")
#         if len(query_ids) < 100:
#             for i in query_ids:
#                 print(f" - {i}")
#     config["query"] = query
#     return query

def parse_input_query(query_arg,ids_arg,cwd,config):
    query = os.path.join(cwd,query_arg)
            
    if not os.path.exists(query):
        if ids_arg:
            id_list = query_arg.split(",")
            query = os.path.join(config["tempdir"], "query.csv")
            with open(query,"w") as fw:
                in_col = config["input_column"]
                fw.write(f"{in_col}\n")
                for i in id_list:
                    fw.write(i+'\n')
        else:
            sys.stderr.write(cyan(f"Error: cannot find query file at {query}\nCheck if the file exists, or if you're inputting an id string (e.g. EPI12345,EPI23456), please use in conjunction with the `--id-string` flag\n."))
            sys.exit(-1)

    ending = query.split(".")[-1]
    if ending in ["yaml","yml"]:
        print(green(f"Input config file:") + f" {query}")
    elif ending == "csv":
        print(green(f"Input file:") + f" {query}")
    config["query"] = query
    return query
        

def check_label_and_colour_fields(query_file, query_arg, colour_fields, label_fields,display_arg, input_column, config):
    acceptable_colours = get_colours()
    queries = []
    colour_field_list = []
    labels = []

    graphics_list = []
    with open(query_file, newline="") as f:
        reader = csv.DictReader(f)
        column_names = reader.fieldnames
        if input_column not in column_names:
            sys.stderr.write(cyan(f"Error: Query file missing header field {input_column}\n"))
            sys.exit(-1)

        if query_arg:
            print(green("Input querys to process:"))
            queries = []
            for row in reader:
                queries.append(row[input_column])
                print(f" - {row[input_column]}")
            print(green(f"Number of queries:") + f" {len(queries)}")

        if not colour_fields:
            colour_field_list.append("adm1")
        else:
            desired_fields = colour_fields.split(",")
            for field in desired_fields:
                if field in column_names:
                    colour_field_list.append(field)
                else:
                    sys.stderr.write(cyan(f"Error: {field} field not found in metadata file\n"))
                    sys.exit(-1)
        
        colour_field_str = ",".join(colour_field_list)
        print(f"Going to colour by: {colour_field_str}")
        config["fields"] = colour_field_str
        
        if not label_fields:
            labels.append("NONE")
        else:
            label_fields = label_fields.split(",")
            for label_f in label_fields:
                if label_f in column_names:
                    labels.append(label_f)
                else:
                    sys.stderr.write(cyan(f"Error: {label_f} field not found in metadata file\n"))
                    sys.exit(-1)

        labels_str = ",".join(labels)
        print(f"Going to label by: {labels_str}")
        config["label_fields"] = labels_str

        if display_arg:
            sections = display_arg.split(",")
            for item in sections:
                splits = item.split("=")
                graphic_trait = splits[0]
                if graphic_trait in column_names:
                    if len(splits) == 1:
                        graphics_list.append(graphic_trait + ":default")
                    else:
                        colour_scheme = splits[1]
                        if colour_scheme in acceptable_colours:
                            graphics_list.append(graphic_trait + ":" + colour_scheme)
                        else:
                            sys.stderr.write(cyan(f"Error: {colour_scheme} not a matplotlib compatible colour scheme\n"))
                            sys.stderr.write(cyan(f"Please use one of {acceptable_colours}"))
                            sys.exit(-1)
                else:
                    sys.stderr.write(cyan(f"Error: {graphic_trait} field not found in metadata file\n"))
                    sys.exit(-1)
        else:
            graphics_list.append("adm1:default")

        config["graphic_dict"] = ",".join(graphics_list)

def map_sequences_config(map_sequences,mapping_trait,x_col,y_col,input_crs,query,config):

        if map_sequences:
            config["map_sequences"] = True
            if not x_col or not y_col:
                sys.stderr.write('Error: coordinates not supplied for mapping sequences. Please provide --x-col and --y-col')
                sys.exit(-1)
            elif not input_crs:
                sys.stderr.write('Error: input coordinate system not provided for mapping. Please provide --input-crs eg EPSG:3395')
                sys.exit(-1)
            else:
                config["x_col"] = x_col
                config["y_col"] = y_col
                config["input_crs"] = input_crs

            with open(query, newline="") as f:
                reader = csv.DictReader(f)
                column_names = reader.fieldnames
                relevant_cols = [x_col, y_col, mapping_trait]
                for map_arg in relevant_cols:

                    if map_arg not in column_names:
                        sys.stderr.write(f"Error: {map_arg} field not found in metadata file")
                        sys.exit(-1)

            if mapping_trait:
                config["mapping_trait"] = mapping_trait
            else:
                config["mapping_trait"] = False
                
        else:
            config["map_sequences"] = False
            config["x_col"] = False
            config["y_col"] = False
            config["input_crs"] = False
            config["mapping_trait"] = False

def local_lineages_config(local_lineages, query, config):
        if local_lineages:
            config['local_lineages'] = "True"
            with open(query, newline="") as f:
                reader = csv.DictReader(f)
                header = reader.fieldnames
                if not "adm2" in header:
                    sys.stderr.write(f"Error: --local-lineages argument called, but input csv file doesn't have an adm2 column. Please provide that to have local lineage analysis.\n")
                    sys.exit(-1)
        else:
            config['local_lineages'] = "False"


def check_summary_fields(full_metadata, summary_field, config):

    with open(full_metadata, newline="") as f:
        reader = csv.DictReader(f)
        column_names = reader.fieldnames

        if not summary_field:
                summary = "lineage"
        else:
            if summary_field in column_names:
                summary = summary_field
            else:
                sys.stderr.write(cyan(f"Error: {summary_field} field not found in metadata file\n"))
                sys.exit(-1)
        
        print(f"Going to summarise collapsed nodes by: {summary}")
        config["node_summary"] = summary


def input_file_qc(fasta,minlen,maxambig,config):
    post_qc_query = ""
    qc_fail = ""
    if fasta:
        do_not_run = []
        run = []
        for record in SeqIO.parse(fasta, "fasta"):
            if len(record) <minlen:
                record.description = record.description + f" fail=seq_len:{len(record)}"
                do_not_run.append(record)
                print(cyan(f"    - {record.id}\tsequence too short: Sequence length {len(record)}"))
            else:
                num_N = str(record.seq).upper().count("N")
                prop_N = round((num_N)/len(record.seq), 2)
                if prop_N > maxambig: 
                    record.description = record.description + f" fail=N_content:{prop_N}"
                    do_not_run.append(record)
                    print(cyan(f"    - {record.id}\thas an N content of {prop_N}"))
                else:
                    run.append(record)

        post_qc_query = os.path.join(config["outdir"], 'query.post_qc.fasta')
        with open(post_qc_query,"w") as fw:
            SeqIO.write(run, fw, "fasta")
        qc_fail = os.path.join(config["outdir"],'query.failed_qc.csv')
        with open(qc_fail,"w") as fw:
            fw.write("name,reason_for_failure\n")
            for record in do_not_run:
                desc = record.description.split(" ")
                for i in desc:
                    if i.startswith("fail="):
                        fw.write(f"{record.id},{i}\n")

    config["post_qc_query"] = post_qc_query
    config["qc_fail"] = qc_fail

def get_package_data(cog_report,thisdir,config):
    reference_fasta = pkg_resources.resource_filename('civet', 'data/reference.fasta')
    outgroup_fasta = pkg_resources.resource_filename('civet', 'data/outgroup.fasta')
    polytomy_figure = pkg_resources.resource_filename('civet', 'data/polytomies.png')
    footer_fig = pkg_resources.resource_filename('civet', 'data/footer.png')
    clean_locs = pkg_resources.resource_filename('civet', 'data/mapping_files/adm2_cleaning.csv')
    map_input_1 = pkg_resources.resource_filename('civet', 'data/mapping_files/gadm36_GBR_2.json')
    map_input_2 = pkg_resources.resource_filename('civet', 'data/mapping_files/channel_islands.json')  
    map_input_3 = pkg_resources.resource_filename('civet', 'data/mapping_files/NI_counties.geojson')  
    map_input_4 = pkg_resources.resource_filename('civet', 'data/mapping_files/Mainland_HBs_gapclosed_mapshaped_d3.json')
    map_input_5 = pkg_resources.resource_filename('civet', 'data/mapping_files/urban_areas_UK.geojson')
    spatial_translations_1 = pkg_resources.resource_filename('civet', 'data/mapping_files/HB_Translation.pkl')
    spatial_translations_2 = pkg_resources.resource_filename('civet', 'data/mapping_files/adm2_regions_to_coords.csv')
    config["reference_fasta"] = reference_fasta
    config["outgroup_fasta"] = outgroup_fasta
    config["polytomy_figure"] = polytomy_figure
    config["footer"] = footer_fig
    
    config["clean_locs"] = clean_locs
    config["uk_map"] = map_input_1
    config["channels_map"] = map_input_2
    config["ni_map"] = map_input_3
    config["uk_map_d3"] = map_input_4
    config["urban_centres"] = map_input_5
    config["HB_translations"] = spatial_translations_1
    config["PC_translations"] = spatial_translations_2

    if cog_report:
        report_template = os.path.join(thisdir, 'scripts','COG_template.pmd')
    else:
        report_template = os.path.join(thisdir, 'scripts','civet_template.pmd')
    
    if not os.path.exists(report_template):
        sys.stderr.write(cyan(f'Error: cannot find report_template at {report_template}\n'))
        sys.exit(-1)
    config["report_template"] = report_template

    