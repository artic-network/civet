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

from class_definitions import taxon,lineage

def convert_date(date_string):
    bits = date_string.split("-")
    date_dt = dt.date(int(bits[0]),int(bits[1]), int(bits[2]))
    
    return date_dt

def parse_tree_tips(tree_dir):

    tips = []
    tip_to_tree = {}

    for fn in os.listdir(tree_dir):
        if fn.endswith("tree"):
            tree_name = fn.split(".")[0]
            tree = bt.loadNewick(tree_dir + "/" + fn, absoluteTime=False)
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

    return tips, tip_to_tree

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



def parse_input_csv(input_csv, query_id_dict, desired_fields, label_fields, date_fields, adm2_adm1_dict): 
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
                elif "collection_date" in col_names:
                    if sequence["collection_date"] != "": #for the COG report 
                        taxon.sample_date = sequence["collection_date"]

                for col in col_names: #Add other metadata fields provided
                    if col in label_fields:
                        if sequence[col] != "":
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

                        if col == "adm2":
                            tax.attribute_dict["adm2"] = adm2
                            if "adm1" not in col_names:
                                if sequence[col] in adm2_adm1_dict.keys():
                                    adm1 = adm2_adm1_dict[sequence[col]]
                                    taxon.attribute_dict["adm1"] = adm1

                new_query_dict[taxon.name] = taxon

      
    return new_query_dict 


def parse_full_metadata(query_dict, label_fields, tree_fields, full_metadata, present_lins, present_in_tree, node_summary_option, date_fields):

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
                # elif "adm2" not in tax_object.attribute_dict.keys() and adm2 == "":
                #     tax_object.attribute_dict['adm2'] = "NA"

                for field in date_fields:
                    if field in reader.fieldnames:
                        if sequence[field] != "" and field not in tax_object.date_dict.keys():
                            date_dt = convert_date(sequence[field])
                            tax_object.date_dict[field] = date_dt 
                    
                for field in label_fields:
                    if field in col_names:
                        if tax_object.attribute_dict[field] == "NA" and sequence[field] != "NA": #this means it's not in the input file
                                tax_object.attribute_dict[field] = sequence[field]


                full_tax_dict[seq_name] = tax_object
                    
    return full_tax_dict
    

def parse_all_metadata(treedir, filtered_cog_metadata, full_metadata_file, input_csv, label_fields, tree_fields, date_fields, node_summary_option, adm2_to_adm1):

    present_in_tree, tip_to_tree = parse_tree_tips(treedir)
    
    #parse the metadata with just those queries found in cog
    query_dict, query_id_dict, present_lins, tree_to_tip = parse_filtered_metadata(filtered_cog_metadata, tip_to_tree, label_fields, tree_fields) 

    if input_csv != '':
         #Any query information they have provided
        query_dict = parse_input_csv(input_csv, query_id_dict, tree_fields, label_fields, date_fields, adm2_to_adm1)
    
    #parse the full background metadata
    full_tax_dict = parse_full_metadata(query_dict, label_fields, tree_fields, full_metadata_file, present_lins, present_in_tree, node_summary_option, date_fields)

    return full_tax_dict, query_dict, tree_to_tip

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

        if not query.in_cog and not cog_report: #poss just needs to be if not query.in_cog now
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
                







