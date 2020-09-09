import geopandas
import pandas as pd
import math
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
from shapely.geometry.point import Point
import csv
import adjustText as aT


def prep_data(tax_dict, clean_locs_file):

    adm2s = []

    for tax in tax_dict.values():
        if tax.attribute_dict["adm2"] != "":
        
            adm2s.append(tax.attribute_dict["adm2"])

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

def prep_mapping_data(mapping_input, metadata_multi_loc):

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


def make_centroids(result,adm2s, straight_map):

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
                centroids = centroid_dict[new]
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


def map_adm2(tax_dict, clean_locs_file, mapping_json_files): #So this takes adm2s and plots them onto the whole UK

    adm2s, metadata_multi_loc, straight_map = prep_data(tax_dict, clean_locs_file)

    all_uk, result = prep_mapping_data(mapping_json_files, metadata_multi_loc)

    centroid_geo = make_centroids(result, adm2s, straight_map)
    
    if type(centroid_geo) == bool:
        print("None of the sequences provided have adequate adm2 data and so cannot be mapped")
        return

    make_map(centroid_geo, all_uk)

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
            
            if colour_map_trait != "False":
                trait = seq[colour_map_trait]
            
            if x != "" and y != "":
                #If we have the actual coordinates
                name_to_coords[name] = (float(x),float(y))

                if colour_map_trait != "False":
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
            
            if colour_map_trait != "False":
                trait = seq[colour_map_trait]
            
            if outer_postcode != "":
                if outer_postcode in pc_to_coords.keys():
                    name_to_coords[name] = pc_to_coords[outer_postcode]

                    if colour_map_trait != "False":
                        name_to_trait[name] = trait

                else:
                    pass


    return name_to_coords, name_to_trait


def plot_coordinates(mapping_json_files, urban_centres, name_to_coords, name_to_trait, input_crs, colour_map_trait):

    ##MAKE DATAFRAME##

    all_uk = generate_all_uk_dataframe(mapping_json_files)
    all_uk = all_uk.to_crs("EPSG:3395")

    urban = geopandas.read_file(urban_centres)

    df_dict = defaultdict(list)

    for name, point in name_to_coords.items():
        df_dict["geometry"].append(Point(point))
        if colour_map_trait != "False":
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

    if colour_map_trait != "False":
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

    return adm2_counter, adm2_percentages
        

def map_sequences_using_coordinates(input_csv, mapping_json_files, urban_centres, pc_file,colour_map_trait, map_inputs, input_crs):

    cols = map_inputs.split(",")
    if len(cols) == 2:
        x_col = cols[0]
        y_col = cols[1]
        name_to_coords, name_to_trait = get_coords_from_file(input_csv, input_crs, colour_map_trait, x_col, y_col)
    elif len(cols) == 1:
        postcode_col = cols[0]
        name_to_coords, name_to_trait = generate_coords_from_outer_postcode(pc_file, input_csv, postcode_col, colour_map_trait)
    
    adm2_counter, adm2_percentages = plot_coordinates(mapping_json_files, urban_centres, name_to_coords, name_to_trait, input_crs, colour_map_trait)

    return adm2_counter, adm2_percentages