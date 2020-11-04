import geopandas
import pandas as pd
import math
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
from shapely.geometry.point import Point
import csv
import adjustText as aT
import os


def generate_all_uk_dataframe(mapping_input):

    uk_file = mapping_input[0]
    channel_file = mapping_input[1]
    ni_file = mapping_input[2]

    UK = geopandas.read_file(uk_file)
    NI = geopandas.read_file(ni_file)
    channels = geopandas.read_file(channel_file)

    ###GET NI COUNTIES

    ni_name = []
    for i in range(len(NI["CountyName"])):
        ni_name.append("Northern Ireland C")

    NI["NAME_2"] = NI["CountyName"]
    NI["NAME_1"] = ni_name  

    all_uk = UK.append(channels).append(NI)

    return all_uk


def find_ambiguities(adm2s):

    ambiguous = []
    ambiguous_dict = defaultdict(set)
    clusters = []

    for adm2 in adm2s:
        if "|" in adm2:
            ambiguous.append(set(adm2.split("|")))

    for group in ambiguous:
        for group2 in ambiguous:
            if group & group2:
                group |= group2

        clusters.append(group)

    for cluster in clusters:
        for place in cluster:
            ambiguous_dict[place] = "|".join(sorted(cluster))


    return ambiguous_dict


def prep_mapping_data(mapping_input, tax_dict):

    adm2s = []

    for tax in tax_dict.values():
        if tax.attribute_dict["adm2"] != "":
            adm2s.append(tax.attribute_dict["adm2"]) #should already be upper and with underscores

    if len(adm2s) == 0:
        return False

    ambiguous_dict = find_ambiguities(adm2s)

    all_uk = generate_all_uk_dataframe(mapping_input)

    ###CAPITALISE GEOJSON DATA

    uppers = []
    for i in all_uk["NAME_2"]:
        uppers.append(i.upper().replace(" ","_"))

    all_uk["NAME_2"] = uppers

    original = all_uk.copy()

    ##DEAL WITH MERGED LOCATIONS EG WEST MIDLANDS

    for_merging = []

    for location in all_uk["NAME_2"]:
        if location in ambiguous_dict:
            for_merging.append(ambiguous_dict[location])
        else:
            for_merging.append(location)

    all_uk["Multi_loc"] = for_merging

    merged_locs = all_uk.dissolve(by="Multi_loc")

    mergeds = []
    for multi_loc in merged_locs.index:
        mergeds.append(multi_loc)

    merged_locs["NAME_2"] = mergeds    

    ###ADD MERGED AND NON-MERGED DATABASES TOGETHER

    result = pd.merge(merged_locs, original, how="outer")

    return all_uk, result, adm2s, ambiguous_dict

def make_centroids_get_counts(result, adm2s, ambiguous_dict):

    not_mappable = ["WALES", "OTHER", "UNKNOWN", "UNKNOWN_SOURCE", "NOT_FOUND", "GIBRALTAR", "FALKLAND_ISLANDS", "CITY_CENTRE"]

    centroid_df = defaultdict(list)
    centroid_dict = {}

    for name, geometry in zip(result["NAME_2"], result["geometry"]):
        centroid_dict[name.upper()] = geometry.centroid

    centroid_counts = Counter(adm2s)
    to_remove = set()

    for k,v in centroid_counts.items():
        if k in ambiguous_dict:
            testing = ambiguous_dict[k]
            for location in centroid_counts.keys():
                if "|" in location:
                    if any([i for i in testing.split("|") if i in location.split("|")]): #in case it's in there in a different order to the value in the ambiguity dict
                        centroid_counts[location] += v
                        to_remove.add(k)
                        break

    for loc in to_remove:
        del centroid_counts[loc]


    for adm2, count in centroid_counts.items():
        centroid_df["Adm2"].append(adm2)
        centroid_df["geometry"].append(centroid_dict[adm2])
        centroid_df["seq_count"].append(count)

    centroid_geo = geopandas.GeoDataFrame(centroid_df)

    return centroid_geo, centroid_counts

def prep_data_old(tax_dict, clean_locs_file):

    adm2s = []

    for tax in tax_dict.values():
        if tax.attribute_dict["adm2"] != "":
            adm2s.append(tax.attribute_dict["adm2"].upper())

    metadata_multi_loc = {}
    straight_map = {}

    with open(clean_locs_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split(",")
            toks [:] = [x for x in toks if x]
            metadata_loc = toks[0]
            real_locs = toks[1:]   
            
            if metadata_loc == 'RHONDDA CYNON TAF':
                straight_map[metadata_loc] = "RHONDDA, CYNON, TAFF" 
            else:
                if len(real_locs) == 1:
                    straight_map[metadata_loc] = real_locs[0].upper()
                else:
                    for i in real_locs:
                        metadata_multi_loc[i.upper()] = metadata_loc.upper()
    
    return adm2s, metadata_multi_loc, straight_map

def prep_mapping_data_old(mapping_input, metadata_multi_loc):

    all_uk = generate_all_uk_dataframe(mapping_input)
    
    ###CAPITALISE GEOJSON DATA

    uppers = []
    for i in all_uk["NAME_2"]:
        uppers.append(i.upper())
        
    all_uk["NAME_2"] = uppers

    original = all_uk.copy()

    ##DEAL WITH MERGED LOCATIONS EG WEST MIDLANDS

    for_merging = []

    for location in all_uk["NAME_2"]:
        if location in metadata_multi_loc.keys():
            new_loc = metadata_multi_loc[location]
        else:
            new_loc = location.upper()
            
        for_merging.append(new_loc)

    all_uk["Multi_loc"] = for_merging

    merged_locs = all_uk.dissolve(by="Multi_loc")

    mergeds = []
    for multi_loc in merged_locs.index:
        mergeds.append(multi_loc)

    merged_locs["NAME_2"] = mergeds

    ###ADD MERGED AND NON-MERGED DATABASES TOGETHER

    result = pd.merge(merged_locs, original, how="outer")

    return all_uk, result


def make_centroids_old(result,adm2s, straight_map):

    not_mappable = ["WALES", "OTHER", "UNKNOWN", "UNKNOWN SOURCE", "NOT FOUND", "GIBRALTAR", "FALKLAND ISLANDS", "CITY CENTRE"]

    centroid_dict = {}

    for name, geometry in zip(result["NAME_2"], result["geometry"]):
        centroid_dict[name.upper()] = geometry.centroid

    centroids = []

    centroid_counts = Counter(adm2s)

    centroid_df = defaultdict(list)

    for adm2, count in centroid_counts.items():
        try:
            if adm2 in straight_map.keys():
                new = straight_map[adm2]
                centroid = centroid_dict[new]
            else:
                centroid = centroid_dict[adm2]
        except KeyError:
            if adm2 != "" and adm2 not in not_mappable:
                print(adm2 + " is not associated with an correct adm2 region so cannot be plotted yet.")
                
        try:
            centroid_df["Adm2"].append(adm2)
            centroid_df["geometry"].append(centroid)
            centroid_df["seq_count"].append(count)
        except UnboundLocalError:
            return False
        
        
    centroid_geo = geopandas.GeoDataFrame(centroid_df)

    return centroid_geo, centroid_counts

def make_map(centroid_geo, all_uk, figdir):

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(20, 15)

    all_uk = all_uk.to_crs("EPSG:3395")
    centroid_geo.crs = "EPSG:4326"
    centroids_final = centroid_geo.to_crs(all_uk.crs)

    base = all_uk.plot(ax=ax, color="steelblue")

    centroids_final.plot(ax=base, color="goldenrod", markersize=centroids_final["seq_count"]*10)

    ax.axis("off")

    plt.savefig(figdir + "/map_adm2.svg", format="svg")



def map_adm2(tax_dict, clean_locs_file, mapping_json_files, figdir, old_data): #So this takes adm2s and plots them onto the whole UK

    if old_data:
        adm2s, metadata_multi_loc, straight_map = prep_data_old(tax_dict, clean_locs_file)
        all_uk, result = prep_mapping_data_old(mapping_json_files, metadata_multi_loc)
        output = make_centroids_old(result, adm2s, straight_map)

    
        if type(output) == bool:
            print("None of the sequences provided have adequate adm2 data and so cannot be mapped")
            return
        else:
            centroid_geo, adm2_counter = output
    else:
        output = prep_mapping_data(mapping_json_files, tax_dict)

        if type(output) == bool:
            print("None of the sequences provided have adequate adm2 data and so cannot be mapped")
            return
        else:
            all_uk, result, adm2s, ambiguous_dict = output

        centroid_geo, adm2_counter = make_centroids_get_counts(result, adm2s, ambiguous_dict)

    make_map(centroid_geo, all_uk, figdir)

    adm2_percentages = {}

    total = len(adm2_counter)

    adm2_to_label = {}
    for taxa in tax_dict.values():
        if taxa.attribute_dict["adm2"].upper() != taxa.attribute_dict["location_label"].upper():
            adm2_to_label[taxa.attribute_dict["adm2"]] =  taxa.attribute_dict["location_label"]

    for adm2, count in adm2_counter.items():
        adm2_percentages[adm2] = round(((count/total)*100),2)

    return adm2_counter, adm2_percentages, adm2_to_label

def get_coords_from_file(input_csv, input_crs, colour_map_trait, x_col, y_col):

    ##READ IN TRAITS##

    name_to_coords = {}
    name_to_trait = {}

    with open(input_csv) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        
        for seq in data:
            name = seq["name"]
            x = seq[x_col]
            y = seq[y_col]
            
            if colour_map_trait:
                trait = seq[colour_map_trait]
            
            if x != "" and y != "":
                #If we have the actual coordinates
                name_to_coords[name] = (float(x),float(y))

                if colour_map_trait:
                    name_to_trait[name] = trait

    return name_to_coords, name_to_trait

def generate_coords_from_outer_postcode(pc_file, input_csv, postcode_col, colour_map_trait):

    pc_to_coords = {}
    name_to_coords = {}
    name_to_trait = {}

    with open(pc_file) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        
        for line in data:

            pc = line["outcode"]
            x = float(line["longitude"])
            y = float(line["latitude"])

            pc_to_coords[pc] = ((x,y))
    
    
    with open(input_csv) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        
        for seq in data:
            name = seq["name"]
            outer_postcode = seq[postcode_col]
            
            if colour_map_trait:
                trait = seq[colour_map_trait]
            
            if outer_postcode != "":
                if outer_postcode in pc_to_coords.keys():
                    name_to_coords[name] = pc_to_coords[outer_postcode]

                    if colour_map_trait:
                        name_to_trait[name] = trait

                else:
                    pass


    return name_to_coords, name_to_trait


def plot_coordinates(mapping_json_files, urban_centres, name_to_coords, name_to_trait, input_crs, colour_map_trait, figdir):

    ##MAKE DATAFRAME##

    all_uk = generate_all_uk_dataframe(mapping_json_files)
    all_uk = all_uk.to_crs("EPSG:3395")

    urban = geopandas.read_file(urban_centres)

    df_dict = defaultdict(list)

    for name, point in name_to_coords.items():
        df_dict["geometry"].append(Point(point))
        if colour_map_trait:
            df_dict[colour_map_trait].append(name_to_trait[name])
        
    crs = {'init':input_crs}

    df = geopandas.GeoDataFrame(df_dict, crs=crs)    
    df_final = df.to_crs(all_uk.crs)

    ##IDENTIFY WHICH ADM2 ARE PRESENT##

    adm2_present = []

    for i in df_final["geometry"]:
        for l,j in zip(all_uk["NAME_2"], all_uk["geometry"]):
            if j.contains(i):
                adm2_present.append(l)

    if len(adm2_present) == 0:
        return

    adm2_counter = Counter(adm2_present)

    total=len(adm2_present)

    adm2_percentages = {}

    for adm2, count in adm2_counter.items():
        adm2_percentages[adm2] = round(((count/total)*100),2)

    ##PREP THE DIFFERENT LAYERS##

    filtered = all_uk[all_uk.NAME_2.isin(list(adm2_counter.keys()))]

    filtered_shape = filtered.dissolve(by="NAME_0")
    filtered_urban = urban[(urban["geometry"].bounds["minx"] > float(filtered_shape["geometry"].bounds.minx)) & (urban["geometry"].bounds["maxx"] < float(filtered_shape["geometry"].bounds.maxx)) & (urban["geometry"].bounds["miny"] > float(filtered_shape["geometry"].bounds.miny)) & (urban["geometry"].bounds["maxy"] < float(filtered_shape["geometry"].bounds.maxy))]

    expanded_filter = all_uk[(all_uk["geometry"].bounds["minx"] > float(filtered_shape["geometry"].bounds.minx)) & (all_uk["geometry"].bounds["maxx"] < float(filtered_shape["geometry"].bounds.maxx)) & (all_uk["geometry"].bounds["miny"] > float(filtered_shape["geometry"].bounds.miny)) & (all_uk["geometry"].bounds["maxy"] < float(filtered_shape["geometry"].bounds.maxy))]

    ##MAKE MAP##

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(10, 10)

    base = filtered.plot(ax=ax, color="whitesmoke", edgecolor="darkgrey")
    expanded_filter.plot(ax=base, color="whitesmoke", edgecolor="darkgrey")
    filtered_urban.plot(ax=base, color="lightgrey")

    if colour_map_trait:
        df_final.plot(ax=base, column=colour_map_trait, legend=True, markersize=10, legend_kwds={"fontsize":10, "bbox_to_anchor":(1.8,1), 'title':colour_map_trait, 'title_fontsize':10})
    else:
        df_final.plot(ax=base, markersize=10)

    ##ADD TEXT##

    filtered["rep"] = filtered["geometry"].representative_point()
    filtered_wth_centre = filtered.copy()
    filtered_wth_centre.set_geometry("rep", inplace=True)

    texts = []
    adm2s_already = []

    for x, y, label in zip(filtered_wth_centre.geometry.x, filtered_wth_centre.geometry.y, filtered_wth_centre["NAME_2"]):
        texts.append(plt.text(x, y, label, fontsize = 15))
        adm2s_already.append(label)

    expanded_filter["rep"] = expanded_filter["geometry"].representative_point()
    filtered_wth_centre2 = expanded_filter.copy()
    filtered_wth_centre2.set_geometry("rep", inplace=True)

    for x, y, label in zip(filtered_wth_centre2.geometry.x, filtered_wth_centre2.geometry.y, filtered_wth_centre2["NAME_2"]):
        if label not in adm2s_already:
            texts.append(plt.text(x, y, label, fontsize = 15))

    aT.adjust_text(texts, force_points=0.3, force_text=0.8, expand_points=(1,1), expand_text=(1,1))
                # arrowprops=dict(arrowstyle="-", color='grey', lw=0.5))


    ax.axis("off") 

    plt.savefig(figdir + "/map_coordinate_or_outer_postcode.svg", format='svg')

    return adm2_counter, adm2_percentages
        

def map_sequences_using_coordinates(input_csv, mapping_json_files, urban_centres, pc_file,colour_map_trait, map_inputs, input_crs, figdir):

    cols = map_inputs.split(",")
    if len(cols) == 2:
        x_col = cols[0].replace(" ","")
        y_col = cols[1].replace(" ","")
        name_to_coords, name_to_trait = get_coords_from_file(input_csv, input_crs, colour_map_trait, x_col, y_col)
    elif len(cols) == 1:
        postcode_col = cols[0]
        name_to_coords, name_to_trait = generate_coords_from_outer_postcode(pc_file, input_csv, postcode_col, colour_map_trait)
    
    mapping_output = plot_coordinates(mapping_json_files, urban_centres, name_to_coords, name_to_trait, input_crs, colour_map_trait, figdir)

    return mapping_output


def convert_str_to_list(string, rel_dir):

    lst = []
    
    prep = string.split(",") 
    for i in prep:
        new = i.replace("'","").replace(" ","")
        newer = new.strip("[").strip("]")
        if rel_dir:
            even_newer = newer.split("/figures")[1]
            final = ("./figures" + even_newer)
        else:
            final = newer
        
        lst.append(final)
    
    return lst

def local_lineages_section(lineage_maps, lineage_tables):

    print("These figures show the background diversity of lineages in the local area to aid with identifying uncommon lineages.")
    
    big_list = convert_str_to_list(lineage_tables,False)
    centralLoc = [t for t in big_list if "_central_" in t][0]
    tableList = [t for t in big_list if "_central_" not in t]
    centralName = centralLoc.split('/')[-1].split("_")[0]

    with open(centralLoc) as f:
        for l in f:
            centralName = l.strip("### ").strip("\n")
            break
    
    linmapList = convert_str_to_list(lineage_maps, True)        

    print(f'Based on the sample density for submitted sequences with adm2 metadata, **{centralName}** was determined to be the focal NHS Health-board.\n')
    print(f'The below figure visualises the relative proportion of assigned UK-Lineages for samples sampled in **{centralName}** for the defined time-frame.')
    print ("![]("+linmapList[0]+")\n")
    print(f'The below figure visualises the relative proportions of assigned UK-Lineages for samples collected in the whole region for the defined time-frame. Plot-size demonstrates relative numbers of sequences across given NHS healthboards.')
    print ("![]("+linmapList[2]+")\n")
    #print(f'The below figure visualises the relative proportion of assigned UK-Lineages for samples collected and sequenced within neighbouring healthboard regions for the defined time-frame.')
    #print ("![]("+linmapList[1]+")")
    #print('\n')
    print(f'Tabulated lineage data for the **central** health-board region:\n')

    with open(centralLoc[0], 'r') as file:
        count = 0
        for l in file:
            if count != 0:
                l = l.rstrip("\n")
            else:
                l=l
            print(l)
            count += 1
    print(f'Tabulated lineage data for the **neighbouring** health-board regions:\n')

    for each in tableList:
        with open(each, "r") as file: 
            count = 0               
            for l in file:
                if count != 0:
                    l = l.rstrip("\n")
                else:
                    l=l
                
                print(l)
                count += 1