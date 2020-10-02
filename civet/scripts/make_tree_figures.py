#!/usr/bin/env python3
import os
from matplotlib import font_manager as fm, rcParams
import baltic as bt

import re
import copy

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import numpy as np
from scipy.special import binom
import math

import itertools
import requests

from io import StringIO as sio
from io import BytesIO as csio
from Bio import Phylo
from collections import defaultdict

import datetime as dt
from collections import Counter
from collections import defaultdict
import math

thisdir = os.path.abspath(os.path.dirname(__file__))

def find_tallest_tree(input_dir):
    tree_heights = []
    
    for r,d,f in os.walk(input_dir):
        for fn in f:
            if fn.endswith(".tree"):
                num_taxa = 0
                intro_name = ""
                with open(r + '/' + fn,"r") as f:
                    for l in f:
                        l = l.rstrip("\n")

                        if l.startswith(" Dimensions NTax="):

                            num_taxa = int(l.rstrip(";").split("=")[1])
                            intro_name = fn.rstrip(".tree")
                
                if num_taxa > 1:
                    tree_file = os.path.join(r, fn)
                    tree = bt.loadNewick(tree_file,absoluteTime=False)
                    tips = []
                    
                    for k in tree.Objects:
                        if k.branchType == 'leaf':
                            tips.append(k.name)
                  
                    tree_heights.append(tree.treeHeight)
    
    max_height = sorted(tree_heights, reverse=True)[0]
    return max_height

def display_name(tree, tree_name, tree_dir, full_taxon_dict, query_dict, custom_tip_fields):
    for k in tree.Objects:
        if k.branchType == 'leaf':
            name = k.name
            display_name = ""
            if "inserted" in name:
                collapsed_node_info, number_nodes = summarise_collapsed_node_for_label(tree_dir, name, tree_name, full_taxon_dict)
                k.traits["display"] = collapsed_node_info
                k.node_number = number_nodes
            else:
                if name in full_taxon_dict:
                    taxon_obj = full_taxon_dict[name]
                
                    date = taxon_obj.sample_date
                    k.traits["display"] = f"{name}|{date}"
                    
                    if "adm2" in taxon_obj.attribute_dict.keys():
                        if taxon_obj.attribute_dict["adm2"] != "NA" and taxon_obj.attribute_dict["adm2"] != "":
                            adm2 = taxon_obj.attribute_dict["adm2"]
                            k.traits["display"] = f"{name}|{adm2}|{date}"
                        else:
                            k.traits["display"] =  f"{name}|{date}"

                    count = 0
                    if len(custom_tip_fields) > 0: 
                        if name in query_dict.keys(): 
                            for label_element in custom_tip_fields:
                                if count == 0:
                                    display_name = taxon_obj.attribute_dict[label_element] 
                                else:   
                                    display_name = display_name + "|" + taxon_obj.attribute_dict[label_element]
                                count += 1
                            
                            k.traits["display"] = display_name

                    k.node_number = 1
                
                
                else:
                    if name.startswith("subtree"):
                        number = name.split("_")[-1]
                        display_name = f"Tree {number}"
                        k.traits["display"] = display_name
                        k.node_number = 1
                    else:
                        k.traits["display"] = name + "|" + "not in dict"
                        k.node_number = 1


def find_colour_dict(query_dict, trait, colour_scheme):

    attribute_options = set()

    if colour_scheme == "default":
        cmap = cm.get_cmap("Paired")
    else:
        cmap = cm.get_cmap(colour_scheme)

    if trait == "adm1":
        colour_dict = {"Wales":"darkseagreen",
                "England":"indianred",
                "Scotland":"steelblue",
                "Northern_Ireland":"skyblue",
                "NA": "dimgrey"}
        return colour_dict

    else:
        for query in query_dict.values():
            attribute_options.add(query.attribute_dict[trait])

    if colour_scheme == "default":     
        if len(attribute_options) == 2:
            options = sorted(list(attribute_options))
            colour_dict = {options[0]: "indianred",
                            options[1]:"steelblue"}
            if "NA" in options:
                colour_dict["NA"] = "dimgrey"
            return colour_dict

        elif len(attribute_options) == 3:
            options = sorted(list(attribute_options))
            colour_dict = {options[0]: "darkseagreen",
                            options[1]:"indianred",
                            options[2]:"steelblue"}
            if "NA" in options:
                colour_dict["NA"] = "dimgrey"
            return colour_dict

    else:
        #get the right number of colours, then loop through the set
        colour_dict = {}
        count = 0
        colors = cmap(np.linspace(0, 1, len(attribute_options)))
        for option in sorted(attribute_options):
            colour_dict[option] = colors[count]
            count += 1
        if "NA" in attribute_options:
            colour_dict["NA"] = "dimgrey"
        return colour_dict

    
def make_scaled_tree(My_Tree, tree_name, tree_dir, num_tips, colour_dict_dict, desired_fields, tallest_height, taxon_dict, query_dict, custom_tip_labels, graphic_dict):

    display_name(My_Tree, tree_name, tree_dir, taxon_dict, query_dict, custom_tip_labels) 
    My_Tree.uncollapseSubtree()

    if num_tips < 10:
        page_height = num_tips
    else:
        page_height = num_tips/2  

    offset = tallest_height - My_Tree.treeHeight
    space_offset = tallest_height/10
    absolute_x_axis_size = tallest_height+space_offset+space_offset + tallest_height #changed from /3 
    
    tipsize = 40
    c_func=lambda k: 'dimgrey' ## colour of branches
    l_func=lambda k: 'lightgrey' ## colour of dotted lines
    s_func = lambda k: tipsize*5 if k.name in query_dict.keys() else (0 if k.node_number > 1 else tipsize)
    z_func=lambda k: 100
    b_func=lambda k: 2.0 #branch width
    so_func=lambda k: tipsize*5 if k.name in query_dict.keys() else 0
    zo_func=lambda k: 99
    zb_func=lambda k: 98
    zt_func=lambda k: 97
    font_size_func = lambda k: 25 if k.name in query_dict.keys() else 15
    kwargs={'ha':'left','va':'center','size':12}

    #Colour by specified trait. If no trait is specified, they will be coloured by UK country
    #The first trait will colour the tips, and additional dots are added to the right of the tip
    
    if len(graphic_dict) == 1 and "adm1" in graphic_dict.keys(): #if they didn't specify any graphics
        trait = "adm1"
    else:
        key_iterator = iter(graphic_dict.keys())
        trait = next(key_iterator) #so always have the first trait as the first colour dot

    first_trait = trait
    colour_dict = colour_dict_dict[trait]
    cn_func = lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_dict.keys() else 'dimgrey'
    co_func=lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_dict.keys() else 'dimgrey' 
    outline_colour_func = lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_dict.keys() else 'dimgrey' 

    x_attr=lambda k: k.height + offset
    y_attr=lambda k: k.y

    y_values = []
    for k in My_Tree.Objects:
        y_values.append(y_attr(k))
    min_y_prep = min(y_values)
    max_y_prep = max(y_values)
    vertical_spacer = 0.5 
    full_page = page_height + vertical_spacer + vertical_spacer
    min_y,max_y = min_y_prep-vertical_spacer,max_y_prep+vertical_spacer

    x_values = []
    for k in My_Tree.Objects:
        x_values.append(x_attr(k))
    max_x = max(x_values)
    
    
    fig,ax = plt.subplots(figsize=(20,page_height),facecolor='w',frameon=False, dpi=200)
    
    My_Tree.plotTree(ax, colour_function=c_func, x_attr=x_attr, y_attr=y_attr, branchWidth=b_func)
    
    My_Tree.plotPoints(ax, x_attr=x_attr, colour_function=cn_func,y_attr=y_attr, size_function=s_func, outline_colour=outline_colour_func)
    My_Tree.plotPoints(ax, x_attr=x_attr, colour_function=co_func, y_attr=y_attr, size_function=so_func, outline_colour=outline_colour_func)

    blob_dict = {}

    for k in My_Tree.Objects:
        if "display" in k.traits:
            name=k.traits["display"]
            
            x=x_attr(k)
            y=y_attr(k)
        
            if k.node_number > 1:
                new_dot_size = tipsize*(1+math.log(k.node_number)) 
                ax.scatter(x, y, s=new_dot_size, marker="s", zorder=3, color="dimgrey")

            height = My_Tree.treeHeight+offset
            text_start = tallest_height+space_offset+space_offset

            if len(desired_fields) > 1:
                
                division = (text_start - tallest_height)/(len(desired_fields))
                tip_point = tallest_height+space_offset

                if k.name in query_dict.keys():
                    
                    count = 0
                    
                    for trait in desired_fields:
                        
                        if trait != first_trait:

                            x_value = tip_point + count
                            count += division

                            option = query_dict[k.name].attribute_dict[trait]
                            
                            if trait in graphic_dict.keys():
                                colour_dict = colour_dict_dict[trait]
                                trait_blob = ax.scatter(x_value, y, tipsize*5, color=colour_dict[option])  
                            else:
                                trait_text = ax.text(x_value, y, option, size=15, ha="left", va="center", fontweight="light")
                            
                            blob_dict[trait] = x_value

                    ax.text(text_start+division, y, name, size=font_size_func(k), ha="left", va="center", fontweight="light")
                    
                    if x != max_x:
                        ax.plot([x+space_offset,tallest_height],[y,y],ls='--',lw=1,color=l_func(k))

                else:

                    ax.text(text_start+division, y, name, size=font_size_func(k), ha="left", va="center", fontweight="light")
                    if x != max_x:
                        ax.plot([x+space_offset,tallest_height],[y,y],ls='--',lw=1,color=l_func(k))

                #This section adds a line in between each trait in the tree
                # for blob_x in blob_dict.values():
                #     line_x = blob_x - (division/2)
                #     ax.plot([line_x,line_x],[min_y,max_y],ls='--',lw=3,color=l_func(k))
            
            
            else:
                ax.text(text_start, y, name, size=font_size_func(k), ha="left", va="center", fontweight="ultralight")
                ax.plot([x+space_offset,tallest_height+space_offset],[y,y],ls='--',lw=1,color=l_func(k))

    #Adds labels to the top of the tree to indicate what each labelled trait is
    if len(desired_fields) > 1:

        blob_dict[first_trait] = tallest_height
        
        for trait, blob_x in blob_dict.items():
            y = max_y
            x = blob_x

            ax.text(x,y,trait, rotation=90, size=15,ha="center", va="bottom")
    
    if num_tips < 10:
        fig2,ax2 =  plt.subplots(figsize=(20,page_height/5),facecolor='w',frameon=False, dpi=200)
    else:
        fig2,ax2 =  plt.subplots(figsize=(20,page_height/10),facecolor='w',frameon=False, dpi=200)

    length = 0.00003

    ax2.plot([0,length], [0.5,0.5], ls='-', lw=1, color="dimgrey")
    ax2.text(0.000015,0.15,"1 SNP",size=20, ha="center", va="center")

    ax.spines['top'].set_visible(False) ## make axes invisible
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax2.spines['top'].set_visible(False) ## make axes invisible
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    ax.set_xlim(-space_offset,absolute_x_axis_size)
    ax.set_ylim(min_y,max_y)
    ax2.set_xlim(-space_offset,absolute_x_axis_size)
    ax2.set_ylim(0,1)

    fig.tight_layout()



def sort_trees_index(tree_dir):
    b_list = []
    d_list = []
    for r,d,f in os.walk(tree_dir):
        for thing in f:
            if thing.endswith("tree"):
                a = thing.split(".")[0]
                b = a.split("_")[-1]
                b_list.append(int(b))
        
    c = sorted(b_list, key=int)
        
    return c

def make_all_of_the_trees(input_dir, tree_name_stem, taxon_dict, query_dict, desired_fields, custom_tip_labels, graphic_dict, min_uk_taxa=3):

    tallest_height = find_tallest_tree(input_dir)

    too_tall_trees = []
    colour_dict_dict = defaultdict(dict)

    overall_df_dict = defaultdict(dict)

    overall_tree_count = 0
    
    tree_order = sort_trees_index(input_dir)

    for trait, colour_scheme in graphic_dict.items():
        colour_dict = find_colour_dict(query_dict, trait, colour_scheme)
        
        colour_dict_dict[trait] = colour_dict

    for fn in tree_order:
        lineage = fn
        treename = f"{tree_name_stem}_{fn}"
        treefile = f"{tree_name_stem}_{fn}.tree"
        nodefile = f"{tree_name_stem}_{fn}"
        num_taxa = 0
        intro_name = ""
        with open(input_dir + "/" + treefile,"r") as f:
            for l in f:
                l = l.rstrip("\n")
                if l.startswith(" Dimensions NTax="):
                    num_taxa = int(l.rstrip(";").split("=")[1])
                    intro_name = fn

        if num_taxa > 1: 
            tree = bt.loadNewick(input_dir + "/" + treefile, absoluteTime=False)

            #make root line
            old_node = tree.root
            new_node = bt.node()
            new_node.children.append(old_node)
            old_node.parent = new_node
            old_node.length=0.000015
            new_node.height = 0
            new_node.y = old_node.y
            tree.root = new_node

            tree.Objects.append(new_node)

            tips = []
            
            for k in tree.Objects:
                if k.branchType == 'leaf':
                    tips.append(k.name)
            if len(tips) < 1000:

                df_dict = summarise_node_table(input_dir, treename, taxon_dict)

                overall_df_dict[treename] = df_dict

                overall_tree_count += 1      
                
                make_scaled_tree(tree, treename, input_dir, len(tips), colour_dict_dict, desired_fields, tallest_height, taxon_dict, query_dict, custom_tip_labels, graphic_dict)     
            
            else:
                too_tall_trees.append(lineage)
                continue

    return too_tall_trees, overall_tree_count, colour_dict_dict, overall_df_dict, tree_order

def summarise_collapsed_node_for_label(tree_dir, focal_node, focal_tree, full_tax_dict): 
    
    focal_tree_file = focal_tree + ".txt"

    with open(tree_dir + "/" + focal_tree_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split("\t")
            node_name = toks[0]
            members = toks[1]
        
            if node_name == focal_node:
                summaries = []
                
                member_list = members.split(",")
                number_nodes = str(len(member_list)) + " nodes"

                for tax in member_list:
                    if tax in full_tax_dict.keys():
                        taxon_obj = full_tax_dict[tax]
                        
                        summaries.append(taxon_obj.node_summary)

                summary_counts = Counter(summaries)

                most_common_summary = []

                if len(summary_counts) > 5:
                    
                    remaining = len(summary_counts) - 5
                    
                    most_common_tups = summary_counts.most_common(5)
                    for i in most_common_tups:
                        most_common_summary.append(i[0])

                    pretty_summary_prep = str(most_common_summary).lstrip("[").rstrip("]").replace("'", "")
                    
                    if remaining == 1:
                        pretty_summary = pretty_summary_prep + " and " + str(remaining) + " other"
                    else:
                        pretty_summary = pretty_summary_prep + " and " + str(remaining) + " others"
                
                else:
                    pretty_summary = str(list(summary_counts.keys())).lstrip("[").rstrip("]").replace("'", "")


                node_number = node_name.lstrip("inserted_node")
                pretty_node_name = "Collapsed node " + node_number

                info = pretty_node_name + ": " + number_nodes + " in " + pretty_summary

    return info, len(member_list)

def summarise_node_table(tree_dir, focal_tree, full_tax_dict):

    focal_tree_file = focal_tree + ".txt"

    df_dict = defaultdict(list)

    with open(tree_dir + "/" + focal_tree_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split("\t")
            node_name = toks[0]
            members = toks[1]
        
            dates = []
            countries = []
            adm2_present = []
            uk_present = False

            node_number = node_name.lstrip("inserted_node")
            
            member_list = members.split(",")

            for tax in member_list:
                if tax in full_tax_dict.keys():
                    taxon_obj = full_tax_dict[tax]
                
                    if taxon_obj.sample_date != "NA":
                        date_string = taxon_obj.sample_date
                        date = dt.datetime.strptime(date_string, "%Y-%m-%d").date()
                        dates.append(date)
                    
                    countries.append(taxon_obj.attribute_dict["country"])

                    if taxon_obj.attribute_dict["country"] == "UK":
                        if "adm2" in taxon_obj.attribute_dict.keys():
                            if taxon_obj.attribute_dict["adm2"] != "" and taxon_obj.attribute_dict["adm2"] != "NA":
                                adm2_present.append(taxon_obj.attribute_dict["adm2"])
            
            if len(adm2_present) != 0:
                adm2_counts = Counter(adm2_present)

            country_counts = Counter(countries)

            most_commons = country_counts.most_common(5)

            country_str = ""

            elem_count = 0

            for country, count in most_commons:
                elem_count += 1
                if elem_count == len(most_commons):
                    elem = country + " (" + str(count) + ")"
                    country_str += elem
                else:
                    elem = country + " (" + str(count) + "), "
                    country_str += elem

            if len(adm2_present) != 0:
                adm2_string = ""
                elem_count = 0
                for adm2, c in adm2_counts.items():
                    elem_count += 1
                    if elem_count == len(adm2_counts):
                        elem = adm2 + " (" + str(c) + ")"
                        adm2_string += elem
                    else:
                        elem = adm2 + " (" + str(c) + "), "
                        adm2_string += elem

            else:
                adm2_string = "NA"
                

            min_date = str(min(dates))
            max_date = str(max(dates))

            if "UK" in countries:
                uk_present = True

            size = len(member_list)

            df_dict["Node number"].append(node_number)
            df_dict["UK present"].append(uk_present)
            df_dict["Number of sequences"].append(size)
            df_dict["Date range"].append(min_date + " to " + max_date)
            df_dict["Countries"].append(country_str)
            df_dict["Admin 2 regions"].append(adm2_string)

    return df_dict

def make_legend(colour_dict_dict):
    num_colours = []
    num_traits = 0
    for trait, colour_dict in colour_dict_dict.items():
        num_colours.append(len(colour_dict))
        num_traits +=1 
    y = 0
    x = 0
    max_colours = sorted(num_colours, reverse=True)[0]
    height = math.sqrt(num_traits)*0.75
    fig,ax = plt.subplots(figsize=(len(colour_dict)+1,height), dpi=700)
    
    for trait, colour_dict in colour_dict_dict.items():
        y +=2
        plt.text(-0.5,y,trait, fontsize=5, ha="right",va="center")
        last_option = ""
        for option in sorted(colour_dict):
            if option == "NA":
                last_option = "NA"
            else:
                x += 1
                col = c=np.array([colour_dict[option]])
                plt.scatter([x], [y], s =10, c=col) #((xloc, yloc), radius) relative to overall plot size
                # ax.add_artist(circle)
                plt.text(x,y-1,option, fontsize=5,ha="center",va="center")
                
        if last_option != "":
            x += 1
            col = c=np.array([colour_dict["NA"]])
            plt.scatter([x], [y], s = 10, c=col) #((xloc, yloc), radius) relative to overall plot size
            # ax.add_artist(circle)
            plt.text(x,y-1,"NA", fontsize=5,ha="center",va="center")
        x = 0
        
    plt.xlim(-1,max_colours+1)
    plt.ylim(0,y+1)

    ax.spines['top'].set_visible(False) ## make axes invisible
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    plt.show()

def describe_tree_background(full_tax_dict, tree_name_stem, tree_dir):

    tree_lst = sort_trees_index(tree_dir)

    hidden_countries = defaultdict(list)

    figure_count = 0

    for fn in tree_lst:
        focal_tree = f"{tree_name_stem}_{fn}"
        focal_tree_file = f"{tree_dir}/{focal_tree}.txt"
        pretty_focal = "Tree " + str(fn)

        collapsed_dict = defaultdict(list)

        with open(focal_tree_file) as f:
            next(f)
            for l in f:
                toks = l.strip("\n").split("\t")
                seqs = toks[1].split(",")

                node_number = toks[0].lstrip("inserted_node")
                new_name = "Collapsed node" + node_number
                collapsed_dict[new_name] = seqs
    
            ndes_country_counts = defaultdict(dict)
            nodes = []

            for nde, seqs in collapsed_dict.items():
                countries = []
                for i in seqs:
                    if i in full_tax_dict.keys():
                        obj = full_tax_dict[i]
                        countries.append(obj.attribute_dict["country"])

                    else:
                        pass

                country_counts = Counter(countries)

                                
                if len(country_counts) > 10:
                    keep_countries = dict(country_counts.most_common(10))
                    if "UK" in countries and "UK" not in keep_countries.keys():
                        keep_countries["UK"] = country_counts["UK"]

                    hidden_countries[focal_tree].append(nde)
                    
                else:
                    keep_countries = country_counts
                    
                if len(country_counts) > 1:
                    
                    ndes_country_counts[nde] = keep_countries
                    nodes.append(nde)
                            
            if len(ndes_country_counts) > 1:
                
                figure_count += 1

                plt.rc('ytick', labelsize=5)
                
                count = 0

                rows = math.ceil(len(ndes_country_counts)/5)
                
                # fig.tight_layout()


                if rows == 1:
                    fig, axs = plt.subplots(rows,5, figsize=(10,2), dpi=250)

                    fig.tight_layout()
                    count = 0      
                    for nde, country_counts in ndes_country_counts.items():

                        x = country_counts.keys()
                        y = country_counts.values()

                        # print("1 counts" + str(count))
                        axs[count].bar(x,y, color="goldenrod")
                        axs[count].set_title(nde, size=8)
                        axs[count].set_xticklabels(x,rotation=90, size=5)
                        #axs[count].set_yticklabels(size=5)
                        axs[count].spines['top'].set_visible(False) ## make axes invisible
                        axs[count].spines['right'].set_visible(False)
                        count += 1

                 
                    fig.suptitle(pretty_focal,y=1.1,x=0.05, size=10)
                
                else:
                    fig, axs = plt.subplots(rows,5,figsize=(10,10), dpi=250)
                    fig.subplots_adjust(hspace=1.0, wspace=0.7)
                    # fig.tight_layout()
                    
                    
                    for nrow in range(0,rows):
                        for i in range(0,5):
                            try:
                                relevant_nde = nodes[(nrow*5) + i]
                                
                                x = ndes_country_counts[relevant_nde].keys()
                                y = ndes_country_counts[relevant_nde].values()

                                axs[nrow][i].bar(x,y, color="goldenrod")
                                axs[nrow][i].set_title(relevant_nde, size=8)
                                axs[nrow][i].set_xticklabels(x,rotation=70, size=5)
                                # axs[nrow][i].set_yticklabels(y, size=5)
                                axs[count].spines['top'].set_visible(False) ## make axes invisible
                                axs[count].spines['right'].set_visible(False)
                            except IndexError:
                                continue

               
                    fig.suptitle(pretty_focal,y=0.95,x=0.1, size=10)

                if len(ndes_country_counts) != rows*5:
                    number_empty_ones = rows*5 - len(ndes_country_counts)
                    for_removal = [i for i in range((rows*5-number_empty_ones),rows*5)]

                    for j in for_removal:
                         fig.delaxes(axs.flatten()[j])

                    
            elif len(ndes_country_counts) == 1:
                
                figure_count += 1
                fig, ax = plt.subplots(figsize=(2,2), dpi=250)

                for nde, country_counts in ndes_country_counts.items():
                    
                    x = country_counts.keys()
                    y = country_counts.values()

                    plt.bar(x,y, color="goldenrod")
                    ax.spines['top'].set_visible(False) ## make axes invisible
                    ax.spines['right'].set_visible(False)
                    # plt.title(nde)
                    plt.xticks(size=5, rotation=90)
                    plt.yticks(size=5)
                    
                    plt.title(pretty_focal + ": " + nde, size=5)


              

    return figure_count

