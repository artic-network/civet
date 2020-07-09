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

def display_name(tree, tree_name, tree_dir, query_id_dict, full_taxon_dict):
    for k in tree.Objects:
        if k.branchType == 'leaf':
            name = k.name
            
            if "inserted" in name:
                collapsed_node_info = summarise_collapsed_node_for_label(tree_dir, name, tree_name, full_taxon_dict)
                k.traits["display"] = collapsed_node_info
            else:
                if name in full_taxon_dict:
                    taxon_obj = full_taxon_dict[name]
                
                    date = taxon_obj.sample_date
                    k.traits["display"] = f"{name}|{date}"
                    
                    if "adm2" in taxon_obj.attribute_dict.keys():
                        adm2 = taxon_obj.attribute_dict["adm2"]
                        k.traits["display"] = f"{name}|{adm2}|{date}"
                
                
                else:
                    k.traits["display"] = name + "|" + "not in dict"


def find_colour_dict(query_dict, trait):

    attribute_options = set()

    cmap = cm.get_cmap("viridis")

    if trait == "adm1":
        colour_dict = {"Wales":"darkseagreen",
                "England":"indianred",
                "Scotland":"steelblue",
                "Northern_Ireland":"skyblue",
                "NA": "goldenrod"}
        return colour_dict

    else:
        for query in query_dict.values():
            attribute_options.add(query.attribute_dict[trait])
            
    if len(attribute_options) == 2:
        colour_dict = {list(attribute_options)[0]: "goldenrod",
                        list(attribute_options)[1]:"midnightblue"}
        return colour_dict

    else:
        #get the right number of colours, then loop through the set
        colour_dict = {}
        count = 0
        colors = cmap(np.linspace(0, 1, len(attribute_options)))
        for option in attribute_options:
            colour_dict[option] = colors[count]
            count += 1
    
        return colour_dict
    
def make_scaled_tree_without_legend(My_Tree, tree_name, tree_dir, num_tips, colour_dict, trait, tallest_height,lineage, taxon_dict, query_id_dict, query_dict):

    display_name(My_Tree, tree_name, tree_dir, query_id_dict, taxon_dict) 
    My_Tree.uncollapseSubtree()

    # closest_names = []

    # for query in query_id_dict.values():
    #     closest_names.append(query.closest)

    if num_tips < 10:
        #page_height = num_tips/2
        page_height = num_tips
    else:
        #page_height = num_tips/4 
        page_height = num_tips/2  

    offset = tallest_height - My_Tree.treeHeight
    space_offset = tallest_height/10
    absolute_x_axis_size = tallest_height+space_offset+space_offset + tallest_height #changed from /3 
    
    tipsize = 40
    c_func=lambda k: 'dimgrey' ## colour of branches
    l_func=lambda k: 'lightgrey' ## colour of branches
    s_func = lambda k: tipsize*5 if k.name in query_id_dict.keys() or k.name in query_dict.keys() else tipsize
    z_func=lambda k: 100
    b_func=lambda k: 2.0 #branch width
    so_func=lambda k: tipsize*5 if k.name in query_id_dict.keys() or k.name in query_dict.keys() else 0
    zo_func=lambda k: 99
    zb_func=lambda k: 98
    zt_func=lambda k: 97
    font_size_func = lambda k: 25 if k.name in query_id_dict.keys() or k.name in query_dict.keys() else 15
    kwargs={'ha':'left','va':'center','size':12}

    #Colour by specified trait. If no trait is specified, they will be coloured by UK country
    cn_func = lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_id_dict.keys() or k.name in query_dict.keys() else 'dimgrey'
    co_func=lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_id_dict.keys() or k.name in query_dict.keys() else 'dimgrey' 
    outline_colour_func = lambda k: colour_dict[query_dict[k.name].attribute_dict[trait]] if k.name in query_id_dict.keys() or k.name in query_dict.keys() else 'dimgrey' 

    x_attr=lambda k: k.height + offset
    y_attr=lambda k: k.y

    y_values = []
    for k in My_Tree.Objects:
        y_values.append(y_attr(k))
    y_values = sorted(y_values)
    vertical_spacer = 0.5 
    full_page = page_height + vertical_spacer + vertical_spacer
    min_y,max_y = y_values[0]-vertical_spacer,y_values[-1]+vertical_spacer
    
    
    fig2,ax2 = plt.subplots(figsize=(20,page_height),facecolor='w',frameon=False, dpi=100)
    

    My_Tree.plotTree(ax2, colour_function=c_func, x_attr=x_attr, y_attr=y_attr, branchWidth=b_func)
    My_Tree.plotPoints(ax2, x_attr=x_attr, colour_function=cn_func,y_attr=y_attr, size_function=s_func, outline_colour=outline_colour_func)
    My_Tree.plotPoints(ax2, x_attr=x_attr, colour_function=co_func, y_attr=y_attr, size_function=so_func, outline_colour=outline_colour_func)


    for k in My_Tree.Objects:
        
        # if k.branchType == 'leaf': 
        #     if k.name in query_dict.keys():
        #         print(k.name)
        # if "country" in k.traits: #don't need these for now because we're only doing UK trees
        #     if k.traits["country"] in colour_dict:
        if "display" in k.traits:
            name=k.traits["display"]

            x=x_attr(k)
            y=y_attr(k)
        
            height = My_Tree.treeHeight+offset
            
            ax2.text(tallest_height+space_offset+space_offset, y, name, size=font_size_func(k), ha="left", va="center", fontweight="ultralight")
            ax2.plot([x+space_offset,tallest_height+space_offset],[y,y],ls='--',lw=0.5,color=l_func(k))

            # if k.name in query_dict.keys() or k.name in query_id_dict.keys():
            #     tree_to_query[tree_name].append(k.name)

    ax2.spines['top'].set_visible(False) ## make axes invisible
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    
    ax2.set_xlim(-space_offset,absolute_x_axis_size)
    ax2.set_ylim(min_y,max_y)

    plt.yticks([])
    plt.xticks([])
    #plt.title(lineage, size=25)

    fig2.tight_layout()



def sort_trees_index(tree_dir):
    b_list = []
    d_list = []
    for r,d,f in os.walk(tree_dir):
        for thing in f:
            if thing.endswith("tree"):
                a = thing.split(".")[0]
                b = a.split("_")[1]
                b_list.append(int(b))
        
    c = sorted(b_list, key=int)
        
    return c

def make_all_of_the_trees(input_dir, taxon_dict, query_id_dict, query_dict, desired_fields, min_uk_taxa=3):

    tallest_height = find_tallest_tree(input_dir)

    too_tall_trees = []
    colour_dict_dict = defaultdict(dict)

    overall_df_dict = defaultdict(dict)

    overall_tree_count = 0
    
    lst = sort_trees_index(input_dir)

    for fn in lst:
        lineage = fn
        treename = "tree_" + str(fn)
        treefile = "tree_" + str(fn) + ".tree"
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

            old_node = tree.root
            new_node = bt.node()
            new_node.children.append(old_node)
            old_node.parent = new_node
            old_node.length=2.0
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
                
                if desired_fields == []:
                    colour_by = ["adm1"]
                
                else:
                    colour_by = desired_fields
            
                for trait in colour_by:
                    colour_dict = find_colour_dict(query_dict, trait)
                    colour_dict_dict[trait] = colour_dict
                    make_scaled_tree_without_legend(tree, treename, input_dir, len(tips), colour_dict, trait, tallest_height, lineage, taxon_dict, query_id_dict, query_dict)     
            else:
                too_tall_trees.append(lineage)
                continue

    return too_tall_trees, overall_tree_count, colour_dict_dict, overall_df_dict

def summarise_collapsed_node_for_label(tree_dir, focal_node, focal_tree, full_tax_dict): 
    
    focal_tree_file = focal_tree + ".txt"

    with open(tree_dir + "/" + focal_tree_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split("\t")
            node_name = toks[0]
            members = toks[1]
        
            if node_name == focal_node:
                countries = []
                
                member_list = members.split(",")
                number_nodes = str(len(member_list)) + " nodes"

                for tax in member_list:
                    if tax in full_tax_dict.keys():
                        taxon_obj = full_tax_dict[tax]
                        
                        countries.append(taxon_obj.attribute_dict["country"])
                    
                    else: #should always be in the full metadata now
                        print("tax missing from full metadata")
                    #     country = tax.split("/")[0]
                    #     countries.append(country)
                    

                country_counts = Counter(countries)

                most_common_countries = []

                if len(country_counts) > 5:
                    
                    remaining = len(country_counts) - 5
                    
                    most_common_tups = country_counts.most_common(5)
                    for i in most_common_tups:
                        most_common_countries.append(i[0])

                    pretty_countries_prep = str(most_common_countries).lstrip("[").rstrip("]").replace("'", "")
                    
                    if remaining == 1:
                        pretty_countries = pretty_countries_prep + " and " + str(remaining) + " other"
                    else:
                        pretty_countries = pretty_countries_prep + " and " + str(remaining) + " others"
                
                else:
                    pretty_countries = str(list(country_counts.keys())).lstrip("[").rstrip("]").replace("'", "")


                node_number = node_name.lstrip("inserted_node")
                pretty_node_name = "Collapsed node " + node_number

                info = pretty_node_name + ": " + number_nodes + " in " + pretty_countries

    return info

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
                            if taxon_obj.attribute_dict["adm2"] != "":
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

def make_legend(colour_dict):
    
    fig,ax = plt.subplots(figsize=(len(colour_dict)+1,1))

    plt.gca().set_aspect('equal', adjustable='box')
    plt.text
    
    x = 0
    for option in colour_dict.keys():
        circle = plt.Circle((x, 0.5), 0.05, color=colour_dict[option]) #((xloc, yloc), radius) relative to overall plot size
        ax.add_artist(circle)
        plt.text(x-0.1,0.3,option, fontsize=5)
        x += 1
        
        
    length = len(colour_dict)

    plt.xlim(-1,length)
    plt.ylim(0,1)

    ax.spines['top'].set_visible(False) ## make axes invisible
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


    plt.yticks([])
    plt.xticks([])
    plt.show()

def describe_tree_background(full_tax_dict, tree_dir):

    tree_lst = sort_trees_index(tree_dir)

    hidden_countries = defaultdict(list)

    figure_count = 0

    for fn in tree_lst:
        focal_tree = "tree_" + str(fn)
        focal_tree_file = tree_dir + "/" + focal_tree + ".txt"
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
                    obj = full_tax_dict[i]
                    countries.append(obj.attribute_dict["country"])

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
                    fig, axs = plt.subplots(rows,5, figsize=(10,2)) 

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
                        
                        count += 1

                 
                    fig.suptitle(pretty_focal,y=1.1,x=0.05, size=10)
                
                else:
                    fig, axs = plt.subplots(rows,5,figsize=(10,10))
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
                plt.figure(figsize=(2,2))

                for nde, country_counts in ndes_country_counts.items():
                    
                    x = country_counts.keys()
                    y = country_counts.values()

                    plt.bar(x,y, color="goldenrod")
                    # plt.title(nde)
                    plt.xticks(size=5, rotation=90)
                    plt.yticks(size=5)



                    plt.title(pretty_focal + ": " + nde, size=5)


              

    return figure_count

