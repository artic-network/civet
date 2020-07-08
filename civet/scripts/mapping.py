import geopandas
import pandas as pd
import math
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt


def prep_data(tax_dict, clean_locs_file):

    adm2s = []

    for tax in tax_dict.values():
        if tax.attribute_dict["adm2"] != "":
        
            adm2s.append(tax.attribute_dict["adm2"])

    merged_locations = defaultdict(list)
    metadata_multi_loc = {}
    straight_map = {}

    with open(clean_locs_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split(",")
            toks [:] = [x for x in toks if x]
            metadata_loc = toks[0]
            real_locs = toks[1:]   
            
            if len(real_locs) == 1:
                straight_map[metadata_loc] = real_locs[0].upper()
            else:
                merged_locations[metadata_loc] = real_locs
                for i in real_locs:
                    metadata_multi_loc[i] = metadata_loc.upper()

    return adm2s, metadata_multi_loc, straight_map

def prep_mapping_data(mapping_input, metadata_multi_loc):

    uk_file = mapping_input[0]
    channel_file = mapping_input[1]
    ni_file = mapping_input[2]

    UK = geopandas.read_file(uk_file)
    NI = geopandas.read_file(ni_file)
    channels = geopandas.read_file(channel_file)

    ni_name = []
    for i in range(len(NI["CountyName"])):
        ni_name.append("Northern Ireland C")

    NI["NAME_2"] = NI["CountyName"]
    NI["NAME_1"] = ni_name  

    all_uk = UK.append(channels).append(NI)


    for_merging = []

    for location in all_uk["NAME_2"]:
        if location in metadata_multi_loc.keys():
            new_loc = metadata_multi_loc[location]
        else:
            new_loc = location.upper()
            
        for_merging.append(new_loc)

    all_uk["Multi_loc"] = for_merging

    merged_locs = all_uk.dissolve(by="Multi_loc")

    result = pd.merge(merged_locs, all_uk, how="outer")

    return all_uk, result


def make_centroids(result,adm2s, straight_map):
    
    centroid_dict = {}

    for name, bigger, geometry in zip(result["NAME_2"], result["Multi_loc"], result["geometry"]):
        if type(bigger) != str:
            centroid_dict[name.upper()] = geometry.centroid
        else:
            centroid_dict[bigger.upper()] = geometry.centroid


    centroids = []

    centroid_counts = Counter(adm2s)

    centroid_df = defaultdict(list)

    for adm2, count in centroid_counts.items():
        if adm2 in straight_map.keys():
            new = straight_map[adm2]
            centroids = centroid_dict[new]
        else:
            centroid = centroid_dict[adm2]
            
        centroid_df["Adm2"].append(adm2)
        centroid_df["geometry"].append(centroid)
        centroid_df["seq_count"].append(count)
        
        
    centroid_geo = geopandas.GeoDataFrame(centroid_df)

    return centroid_geo

def make_map(centroid_geo, all_uk):

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(20, 15)

    all_uk = all_uk.to_crs("EPSG:3395")
    centroid_geo.crs = "EPSG:4326"
    centroids_final = centroid_geo.to_crs(all_uk.crs)

    base = all_uk.plot(ax=ax, color="steelblue")

    centroids_final.plot(ax=base, color='goldenrod', markersize=centroids_final["seq_count"]*10)

    ax.axis("off")


def run_map_functions(tax_dict, clean_locs_file, mapping_json_files):

    adm2s, metadata_multi_loc, straight_map = prep_data(tax_dict, clean_locs_file)

    all_uk, result = prep_mapping_data(mapping_json_files, metadata_multi_loc)

    centroid_geo = make_centroids(result, adm2s, straight_map)

    make_map(centroid_geo, all_uk)