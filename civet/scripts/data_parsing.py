#!/usr/bin/env python3
from collections import defaultdict
from collections import Counter
import pandas as pd
import csv
from tabulate import tabulate
import baltic as bt
import os
import datetime as dt
import math
import matplotlib.pyplot as plt

class taxon():

    def __init__(self, name, global_lin, uk_lin, phylotype, label_fields, tree_fields):

        self.name = name

        self.sample_date = "NA"

        self.date_dict = {}
        
        if global_lin == "":
            self.global_lin = "NA"
        else:
            self.global_lin = global_lin
        
        if uk_lin == "":
            self.uk_lin = "NA"
        else:
            self.uk_lin = uk_lin
       
        if phylotype == "":
            self.phylotype = "NA"
        else:
            self.phylotype = phylotype
       
        self.in_cog = False
        self.attribute_dict = {}
        self.attribute_dict["adm1"] = "NA"
        
        for i in label_fields:
            self.attribute_dict[i] = "NA"
        for i in tree_fields:
            self.attribute_dict[i] = "NA"

        self.tree = "NA"

        self.closest_distance = "NA"
        self.snps = "NA"

class lineage():
    
    def __init__(self, name, taxa):
        
        self.name = name
        self.taxa = taxa
        self.dates = []
        self.global_lins = set()
        
        for tax in taxa:
            if tax.sample_date != "NA":
                self.dates.append(tax.date_dt)
            self.global_lins.add(tax.global_lin)
                
        if self.dates == []:
            self.first_date = "NA"
        else:
            self.first_date = min(self.dates)

def analyse_inputs(inputs):

    desired_fields, label_fields, graphic_dict, node_summary_option, map_sequences, mapping_trait, map_inputs = inputs 

    print("Showing " + ",".join(desired_fields) + " on the tree.")
    print(",".join(list(graphic_dict.keys())) + " fields are displayed graphically using " + ",".join(list(graphic_dict.values())) + " colour schemes respectively.")
    
    if label_fields != "NONE":
        print("Labelling by " + ",".join(label_fields) + " on the tree.")
    
    print("Summarising nodes by " + node_summary_option)

    if map_sequences != "False":
        map_args = map_inputs.split(",")
        if len(map_args) == 2:
            if mapping_trait != "False":
                print("Mapping sequences using columns " + map_args[0] + " " + map_args[1] + " for x values and y values respectively, and colouring by " + mapping_trait)
            else:
                print("Mapping sequences using columns " + map_args[0] + " " + map_args[1] + " for x values and y values respectively.")
        else:
            if mapping_trait != "False":
                print("Mapping sequences using columns " + map_args[0] + " for outer postcodes, and colouring by " + mapping_trait)
            else:
                print("Mapping sequences using columns " + map_args[0] + " for outer postocdes.")

    
    

def prepping_adm2_adm1_data(full_metadata):

    official_adm2_adm1 = {'BARNSLEY': 'England', 'BATH AND NORTH EAST SOMERSET': 'England', 'BEDFORDSHIRE': 'England', 'BIRMINGHAM': 'England', 'BLACKBURN WITH DARWEN': 'England', 'BLACKPOOL': 'England', 'BOLTON': 'England', 'BOURNEMOUTH': 'England', 'BRACKNELL FOREST': 'England', 'BRADFORD': 'England', 'BRIGHTON AND HOVE': 'England', 'BRISTOL': 'England', 'BUCKINGHAMSHIRE': 'England', 'BURY': 'England', 'CALDERDALE': 'England', 'CAMBRIDGESHIRE': 'England', 'CENTRAL BEDFORDSHIRE': 'England', 'CHESHIRE EAST': 'England', 'CHESHIRE WEST AND CHESTER': 'England', 'CORNWALL': 'England', 'COVENTRY': 'England', 'CUMBRIA': 'England', 'DARLINGTON': 'England', 'DERBY': 'England', 'DERBYSHIRE': 'England', 'DEVON': 'England', 'DONCASTER': 'England', 'DORSET': 'England', 'DUDLEY': 'England', 'DURHAM': 'England', 'EAST RIDING OF YORKSHIRE': 'England', 'EAST SUSSEX': 'England', 'ESSEX': 'England', 'GATESHEAD': 'England', 'GLOUCESTERSHIRE': 'England', 'GREATER LONDON': 'England', 'HALTON': 'England', 'HAMPSHIRE': 'England', 'HARTLEPOOL': 'England', 'HEREFORDSHIRE': 'England', 'HERTFORDSHIRE': 'England', 'ISLE OF WIGHT': 'England', 'ISLES OF SCILLY': 'England', 'KENT': 'England', 'KINGSTON UPON HULL': 'England', 'KIRKLEES': 'England', 'KNOWSLEY': 'England', 'LANCASHIRE': 'England', 'LEEDS': 'England', 'LEICESTER': 'England', 'LEICESTERSHIRE': 'England', 'LINCOLNSHIRE': 'England', 'LUTON': 'England', 'MANCHESTER': 'England', 'MEDWAY': 'England', 'MIDDLESBROUGH': 'England', 'MILTON KEYNES': 'England', 'NEWCASTLE UPON TYNE': 'England', 'NORFOLK': 'England', 'NORTH LINCOLNSHIRE': 'England', 'NORTH SOMERSET': 'England', 'NORTH TYNESIDE': 'England', 'NORTH YORKSHIRE': 'England', 'NORTHAMPTONSHIRE': 'England', 'NORTHUMBERLAND': 'England', 'NOTTINGHAM': 'England', 'NOTTINGHAMSHIRE': 'England', 'OLDHAM': 'England', 'OXFORDSHIRE': 'England', 'PETERBOROUGH': 'England', 'PLYMOUTH': 'England', 'POOLE': 'England', 'PORTSMOUTH': 'England', 'READING': 'England', 'REDCAR AND CLEVELAND': 'England', 'ROCHDALE': 'England', 'ROTHERHAM': 'England', 'RUTLAND': 'England', 'SAINT HELENS': 'England', 'SALFORD': 'England', 'SANDWELL': 'England', 'SEFTON': 'England', 'SHEFFIELD': 'England', 'SHROPSHIRE': 'England', 'SLOUGH': 'England', 'SOLIHULL': 'England', 'SOMERSET': 'England', 'SOUTH GLOUCESTERSHIRE': 'England', 'SOUTH TYNESIDE': 'England', 'SOUTHAMPTON': 'England', 'SOUTHEND-ON-SEA': 'England', 'STAFFORDSHIRE': 'England', 'STOCKPORT': 'England', 'STOCKTON-ON-TEES': 'England', 'STOKE-ON-TRENT': 'England', 'SUFFOLK': 'England', 'SUNDERLAND': 'England', 'SURREY': 'England', 'SWINDON': 'England', 'TAMESIDE': 'England', 'TELFORD AND WREKIN': 'England', 'THURROCK': 'England', 'TORBAY': 'England', 'TRAFFORD': 'England', 'WAKEFIELD': 'England', 'WALSALL': 'England', 'WARRINGTON': 'England', 'WARWICKSHIRE': 'England', 'WEST BERKSHIRE': 'England', 'WEST SUSSEX': 'England', 'WIGAN': 'England', 'WILTSHIRE': 'England', 'WINDSOR AND MAIDENHEAD': 'England', 'WIRRAL': 'England', 'WOKINGHAM': 'England', 'WOLVERHAMPTON': 'England', 'WORCESTERSHIRE': 'England', 'YORK': 'England', 'ANTRIM AND NEWTOWNABBEY': 'Northern Ireland', 'ARMAGH, BANBRIDGE AND CRAIGAVON': 'Northern Ireland', 'BELFAST': 'Northern Ireland', 'CAUSEWAY COAST AND GLENS': 'Northern Ireland', 'DERRY AND STRABANE': 'Northern Ireland', 'FERMANAGH AND OMAGH': 'Northern Ireland', 'LISBURN AND CASTLEREAGH': 'Northern Ireland', 'MID AND EAST ANTRIM': 'Northern Ireland', 'MID ULSTER': 'Northern Ireland', 'NEWRY, MOURNE AND DOWN': 'Northern Ireland', 'NORTH DOWN AND ARDS': 'Northern Ireland', 'ABERDEEN': 'Scotland', 'ABERDEENSHIRE': 'Scotland', 'ANGUS': 'Scotland', 'ARGYLL AND BUTE': 'Scotland', 'CLACKMANNANSHIRE': 'Scotland', 'DUMFRIES AND GALLOWAY': 'Scotland', 'DUNDEE': 'Scotland', 'EAST AYRSHIRE': 'Scotland', 'EAST DUNBARTONSHIRE': 'Scotland', 'EAST LOTHIAN': 'Scotland', 'EAST RENFREWSHIRE': 'Scotland', 'EDINBURGH': 'Scotland', 'EILEAN SIAR': 'Scotland', 'FALKIRK': 'Scotland', 'FIFE': 'Scotland', 'GLASGOW': 'Scotland', 'HIGHLAND': 'Scotland', 'INVERCLYDE': 'Scotland', 'MIDLOTHIAN': 'Scotland', 'MORAY': 'Scotland', 'NORTH AYRSHIRE': 'Scotland', 'NORTH LANARKSHIRE': 'Scotland', 'ORKNEY ISLANDS': 'Scotland', 'PERTHSHIRE AND KINROSS': 'Scotland', 'RENFREWSHIRE': 'Scotland', 'SCOTTISH BORDERS': 'Scotland', 'SHETLAND ISLANDS': 'Scotland', 'SOUTH AYRSHIRE': 'Scotland', 'SOUTH LANARKSHIRE': 'Scotland', 'STIRLING': 'Scotland', 'WEST DUNBARTONSHIRE': 'Scotland', 'WEST LOTHIAN': 'Scotland', 'ANGLESEY': 'Wales', 'BLAENAU GWENT': 'Wales', 'BRIDGEND': 'Wales', 'CAERPHILLY': 'Wales', 'CARDIFF': 'Wales', 'CARMARTHENSHIRE': 'Wales', 'CEREDIGION': 'Wales', 'CONWY': 'Wales', 'DENBIGHSHIRE': 'Wales', 'FLINTSHIRE': 'Wales', 'GWYNEDD': 'Wales', 'MERTHYR TYDFIL': 'Wales', 'MONMOUTHSHIRE': 'Wales', 'NEATH PORT TALBOT': 'Wales', 'NEWPORT': 'Wales', 'PEMBROKESHIRE': 'Wales', 'POWYS': 'Wales', 'RHONDDA, CYNON, TAFF': 'Wales', 'SWANSEA': 'Wales', 'TORFAEN': 'Wales', 'VALE OF GLAMORGAN': 'Wales', 'WREXHAM': 'Wales'}

    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIR": "Northern_Ireland"}

    illegal_values = ["", "NOT FOUND", "NONE","OTHER", "WALES", "UNKNOWN", "UNKNOWN SOURCE"]

    test_set = set()

    adm2_adm1 = official_adm2_adm1.copy()

    with open(full_metadata) as f:
        r = csv.DictReader(f)
        in_data = [x for x in r]
        for seq in in_data:
            if seq["country"] == "UK":
                adm1 = contract_dict[seq["adm1"].split("-")[1]]
                adm2 = seq["adm2"]
            
                if adm2.upper() not in illegal_values and adm2 not in official_adm2_adm1:
                    adm2_adm1[adm2] = adm1

    return adm2_adm1
    

def parse_filtered_metadata(metadata_file, tip_to_tree, label_fields, tree_fields):
    
    query_dict = {}
    query_id_dict = {}
    present_lins = set()

    tree_to_tip = defaultdict(list)

    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIR": "Northern_Ireland"}

    with open(metadata_file, "r") as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            glob_lin = sequence['lineage']
            uk_lineage = sequence['uk_lineage']
            
            adm1_prep = sequence["adm1"].split("-")[1]
            adm1 = contract_dict[adm1_prep]

            query_id = sequence['query_id']
            query_name = sequence['query']
            closest_name = sequence["closest"]
            closest_distance = sequence["closest_distance"]
            snps = sequence['snps']

            phylotype = sequence["phylotype"]
            sample_date = sequence["sample_date"]

            
            new_taxon = taxon(query_name, glob_lin, uk_lineage, phylotype, label_fields, tree_fields)

            new_taxon.query_id = query_id

            present_lins.add(uk_lineage)

            if query_name == closest_name: #if it's in COG, get it's sample date
                new_taxon.in_cog = True
                new_taxon.sample_date = sample_date
                new_taxon.attribute_dict["adm1"] = adm1
                new_taxon.closest = "NA"

            else:
                new_taxon.closest = closest_name
                new_taxon.closest_distance = closest_distance
                new_taxon.snps = snps
                for k,v in contract_dict.items():
                    if k in query_name or v in query_name: #if any part of any country name is in the query name it will pick it up assign it
                        new_taxon.attribute_dict["adm1"] = v
                
            new_taxon.attribute_dict["country"] = "UK"

            relevant_tree = tip_to_tree[query_name]
            new_taxon.tree = relevant_tree

            tree_to_tip[relevant_tree].append(new_taxon)
           
            query_dict[query_name] = new_taxon
            query_id_dict[query_id] = new_taxon
            
    return query_dict, query_id_dict, present_lins, tree_to_tip
    

def convert_date(date_string):
    bits = date_string.split("-")
    date_dt = dt.date(int(bits[0]),int(bits[1]), int(bits[2]))
    
    return date_dt


def parse_input_csv(input_csv, query_id_dict, desired_fields, label_fields, date_fields, adm2_adm1_dict, cog_report):
    new_query_dict = {}
    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIR": "Northern_Ireland"}
    cleaning = {"SCOTLAND":"Scotland", "WALES":"Wales", "ENGLAND":"England", "NORTHERN_IRELAND": "Northern_Ireland", "NORTHERN IRELAND": "Northern_Ireland"}

    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())

    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            
            name = sequence["name"]

            if name in query_id_dict.keys():
                taxon = query_id_dict[name]

                for field in date_fields:
                    if field in reader.fieldnames:
                        if sequence[field] != "":
                            date_dt = convert_date(sequence[field])
                            taxon.date_dict[field] = date_dt 

                #keep this separate to above, because sample date is specifically needed
                if "sample_date" in col_names: #if it's not in COG but date is provided (if it's in COG, it will already have been assigned a sample date.)
                    if sequence["sample_date"] != "":
                        taxon.sample_date = sequence["sample_date"]

                for col in col_names: #Add other metadata fields provided
                    if col in label_fields:
                        if sequence[col] == "":
                                taxon.attribute_dict[col] = "NA"     
                        else:
                            taxon.attribute_dict[col] = sequence[col]
                    else:
                        if col != "name" and col in desired_fields and col != "adm1":
                            if sequence[col] != "":
                                taxon.attribute_dict[col] = sequence[col]
            
                        if col == "adm1":
                            if "UK" in sequence[col]:
                                adm1_prep = sequence[col].split("-")[1]
                                adm1 = contract_dict[adm1_prep]
                            else:
                                if sequence[col].upper() in cleaning.keys():
                                    adm1 = cleaning[sequence[col].upper()]
                                else:
                                    adm1 = sequence[col]

                            taxon.attribute_dict["adm1"] = adm1

                        if col == "adm2" and "adm1" not in col_names: #or sequence["adm1"] == ""):
                            if sequence[col] in adm2_adm1_dict.keys():
                                adm1 = adm2_adm1_dict[sequence[col]]
                                taxon.attribute_dict["adm1"] = adm1

                if cog_report:
                    taxon.attribute_dict["adm2"]= sequence["adm2"] 
                    if "collection_date" in reader.fieldnames:
                        date_string = sequence["collection_date"]
                    else:
                        date_string = sequence["sample_date"]

                    taxon.sample_date = date_string
                    
                    bits = date_string.split("-")
                    taxon.date_dt = dt.date(int(bits[0]),int(bits[1]), int(bits[2])) 


                new_query_dict[taxon.name] = taxon
                
    return new_query_dict 


def parse_tree_tips(tree_dir):

    tips = []
    tip_to_tree = {}
    tree_list = []

    for fn in os.listdir(tree_dir):
        if fn.endswith("tree"):
            tree_name = fn.split(".")[0]
            tree = bt.loadNewick(tree_dir + "/" + fn, absoluteTime=False)
            tree_list.append(tree_name)
            for k in tree.Objects:
                if k.branchType == 'leaf' and "inserted" not in k.name:
                    tips.append(k.name)
                    tip_to_tree[k.name] = tree_name

        elif fn.endswith(".txt"):
            with open(tree_dir + "/" + fn) as f:
                for l in f:
                    tip_string = l.strip("\n").split("\t")[1]
                    tip_list = tip_string.split(",")
                    tips.extend(tip_list)

    return tips, tip_to_tree, tree_list

def parse_full_metadata(query_dict, label_fields, tree_fields,full_metadata, present_lins, present_in_tree, node_summary_option, date_fields):

    full_tax_dict = query_dict.copy()

    with open(full_metadata, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())

    with open(full_metadata, 'r') as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            uk_lin = sequence["uk_lineage"]
            seq_name = sequence["sequence_name"]

            date = sequence["sample_date"]
            adm2 = sequence["adm2"]
            country = sequence["country"]

            glob_lin = sequence["lineage"]
            phylotype = sequence["phylotype"]

            node_summary_trait = sequence[node_summary_option]


            if (uk_lin in present_lins or seq_name in present_in_tree) and seq_name not in query_dict.keys():
                new_taxon = taxon(seq_name, glob_lin, uk_lin, phylotype, label_fields, tree_fields)
                if date == "":
                    date = "NA"
                
                new_taxon.sample_date = date

                new_taxon.node_summary = node_summary_trait

                
                new_taxon.attribute_dict["adm2"] = adm2
                new_taxon.attribute_dict["country"] = country

                full_tax_dict[seq_name] = new_taxon

            #There may be sequences not in COG tree but that are in the full metadata, so we want to pull out the additional information if it's not in the input csv
            #Remember, then name has to match fully, so if it's not the x/y/z name this section won't work
            if seq_name in query_dict.keys(): 
                tax_object = query_dict[seq_name]
                if tax_object.sample_date == "NA" and date != "":
                    tax_object.sample_date = date
                    tax_object.all_dates.append(convert_date(date))
                if "adm2" not in tax_object.attribute_dict.keys() and adm2 != "":
                    tax_object.attribute_dict["adm2"] = adm2

                for field in date_fields:
                    if field in reader.fieldnames:
                        if sequence[field] != "" and field not in tax_object.date_dict.keys():
                            date_dt = convert_date(sequence[field])
                            tax_object.date_dict[field] = date_dt 
                    
                for field in label_fields:
                    if field in col_names:
                        if tax_object.attribute_dict[field] == "NA" and sequence[field] != "":
                            tax_object.attribute_dict[field] = sequence[field]

                full_tax_dict[seq_name] = tax_object
                    
    return full_tax_dict
    

def make_initial_table(query_dict, desired_fields, label_fields, cog_report):

    df_dict_incog = defaultdict(list)
    df_dict_seqprovided = defaultdict(list)

    incog = 0
    seqprovided = 0
    incogs = False
    seqprovideds = False

    for query in query_dict.values():

        if query.in_cog:
            df_dict = df_dict_incog
            incog += 1
        else:
            df_dict = df_dict_seqprovided
            seqprovided += 1
        
        df_dict["Query ID"].append(query.query_id.replace("|","\|"))
        
        if query.in_cog: 
            df_dict["Sequence name in Tree"].append(query.name)
        # else:
        #     df_dict["Sequence name in Tree"].append("NA")
        

        df_dict["Sample date"].append(query.sample_date)

        if not query.in_cog and not cog_report:
            df_dict["Closest sequence in Tree"].append(query.closest)
            df_dict["Distance to closest sequence"].append(query.closest_distance)
            df_dict["SNPs"].append(query.snps)

        df_dict["UK lineage"].append(query.uk_lin)
        df_dict["Global lineage"].append(query.global_lin)
        df_dict["Phylotype"].append(query.phylotype)

        
        if query.tree != "NA":
            tree_number = query.tree.split("_")[-1]
            pretty_tree = "Tree " + str(tree_number)
            df_dict["Tree"].append(pretty_tree)
        else:
            df_dict["Tree"].append("NA") #this should never happen, it's more error catching

        if desired_fields != []:
            for i in desired_fields:
                df_dict[i].append(query.attribute_dict[i])
        
        if label_fields != []:
            for i in label_fields: 
                if i not in desired_fields and i != "sample_date" and i != "name":
                    df_dict[i].append(query.attribute_dict[i])

        if cog_report:
            df_dict['adm2'].append(query.attribute_dict["adm2"])

    if incog != 0:
        df_incog = pd.DataFrame(df_dict_incog)
        df_incog.set_index("Query ID", inplace=True)
        incogs = True
    
    if seqprovided != 0:
        df_seqprovided = pd.DataFrame(df_dict_seqprovided)
        df_seqprovided.set_index("Query ID", inplace=True)
        seqprovideds = True

    if seqprovideds and incogs:
        return df_incog, df_seqprovided
    elif seqprovideds and not incogs:
        return df_seqprovided
    elif incogs and not seqprovideds:
        return df_incog

def investigate_QC_fails(QC_file):

    fail_dict = {}

    with open(QC_file) as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            name = sequence["name"]
            reason = sequence["reason_for_failure"]

            if "seq_len" in reason:
                length = reason.split(":")[1]
                final_reason = "Sequence too short: only " + length + " bases."
            elif "N_content" in reason:
                n_content = reason.split(":")[1]
                final_reason = "Sequence has too many Ns: " + str(float(round(float(n_content)*100))) + "\% of bases"

            fail_dict[name] = final_reason


    return fail_dict

def print_missing_seqs(missing_seqs_file):
    
    failed_names = []

    with open(missing_seqs_file) as f:
        for l in f:
            name = l.strip("\n").split(",")[0]
            
            failed_names.append(name)

    return failed_names

def find_new_introductions(query_dict, min_date): #will only be called for the COG sitrep, and the query dict will already be filtered to the most recent sequences

    lin_to_tax = defaultdict(list)
    intro_to_regions = defaultdict(dict)
    df_dict = defaultdict(list)
    new_intros = []
    lin_dict = {}

    df = "" #so that if there are no new ones then it returns an empty string

    for taxon in query_dict.values():
        lin = taxon.uk_lin
        lin_to_tax[lin].append(taxon)

    for k,v in lin_to_tax.items():
        new_lin = lineage(k,v)
        lin_dict[k] = new_lin    

    for key, value in lin_dict.items():
        if value.first_date > min_date:
            new_intros.append(value)


    for intro in new_intros:
        adm2s = []
        trees = set()
        for i in intro.taxa:
            if i.attribute_dict["adm2"] != "":
                adm2s.append(i.attribute_dict["adm2"])
                
            trees.add(i.tree)

        adm2_counts = Counter(adm2s)


        place_string = ""
        for place, count in adm2_counts.items():
            place_string += place + " (" + str(count) + ") "  
        

        df_dict["Name"].append(intro.name)
        df_dict["Size"].append(len(intro.taxa))
        df_dict["Locations"].append(place_string)
        df_dict["Global lineage"].append(intro.global_lins)
        df_dict["Trees"].append(trees)
        

    df = pd.DataFrame(df_dict)

    new_lins = []
    for i in df["Global lineage"]:
        new_lin = str(i).strip("{").strip("}").replace("'","")
        new_lins.append(new_lin)
    df["Global lineage"] = new_lins

    new_trees = []
    for i in df["Trees"]:
        new_tree = str(i).strip("{").strip("}").replace("'","")
        new_trees.append(new_tree)
    df["Trees"] = new_trees

    df.set_index("Name", inplace=True)

    return new_intros, df
                







