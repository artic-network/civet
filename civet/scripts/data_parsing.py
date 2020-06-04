from collections import defaultdict
import pandas as pd
import csv


class taxon():

    def __init__(self, name, global_lin, uk_lin, phylotype):

        self.name = name
        self.global_lin = global_lin
        self.uk_lin = uk_lin
        self.phylotype = phylotype
       
        self.in_cog = False

        self.attribute_dict = {}

def parse_big_metadata(metadata_file):
    
    query_dict = {}
    full_tax_dict = {}

    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIRE": "Northern_Ireland"}

    with open(metadata_file, "r") as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            if sequence['country'] == "UK":

                glob_lin = sequence['lineage']
                uk_lineage = sequence['uk_lineage']
                
                adm1_prep = sequence["adm1"].split("-")[1]
                adm1 = contract_dict[adm1_prep]

                seq_name = sequence['sequence_name']
                closest_name = sequence["closest"]
                query_name = sequence['query']

                phylotype = sequence["phylotype"]

                if query_name != "":

                    new_taxon = taxon(query_name, glob_lin, uk_lineage, phylotype)
                    
                    if query_name == closest_name:
                        new_taxon.in_cog = True

                    query_dict[query_name] = new_taxon
                    full_tax_dict[query_name] = new_taxon
                
                else:

                    new_taxon = taxon(seq_name, glob_lin, uk_lineage)
                    new_taxon.sample_date = sequence["sample_date"]
                    full_tax_dict[seq_name] = new_taxon


    return query_dict, full_tax_dict

def parse_their_metadata(input_metadata, query_dict, desired_fields):

    new_query_dict = {}
    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIRE": "Northern_Ireland"}

    with open(input_metadata, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())
        in_data = [r for r in reader]
        for sequence in in_data:

            name = sequence["sequence_name"]
            sample_date = sequence["sample_date"]

            taxon = query_dict[name]

            taxon.sample_date = sample_date

            for col in col_names:
                if desired_fields != []:
                    if col != "sequence_name" and col != "sample_date" and col in desired_fields:
                        taxon.attribute_dict[col] = sequence[col]
                if col == "adm1":
                    if "UK" in sequence[col]:
                        adm1_prep = sequence[col].split("-")[1]
                        adm1 = contract_dict[adm1_prep]
                    else:
                        adm1 = sequence[col]

                    taxon.adm1 = adm1

            new_query_dict[name] = taxon
                
    return new_query_dict


def make_initial_table(query_dict):

    df_dict = defaultdict(list)

    for query in query_dict.values():

        df_dict["Sequence name"].append(query.name)
        df_dict["Part of COG"].append(query.in_cog)
        df_dict["UK lineage"].append(query.uk_lin)
        df_dict["Global lineage"].append(query.global_lin)
        df_dict["Phylotype"].append(query.phylotype)

        for key, value in query.attribute_dict.items():
            df_dict[key].append(value)

    df = pd.DataFrame(df_dict)

    df.set_index("Sequence name", inplace=True)

    return df

    



            















                

        


