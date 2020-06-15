import os
from matplotlib import font_manager as fm, rcParams
import baltic as bt

import re
import copy

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

import numpy as np
from scipy.special import binom

import itertools
import requests

from io import StringIO as sio
from io import BytesIO as csio
from Bio import Phylo
from collections import defaultdict

import datetime as dt
from collections import Counter


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
                collapsed_node_info = summarise_collapsed_node(tree_dir, name, tree_name, full_taxon_dict)
                k.traits["display"] = collapsed_node_info
            else:
                # taxon_obj = query_id_dict[name] #I think this isn't needed because the ones in COG will just have the full name in the tree
                if name in full_taxon_dict:
                    taxon_obj = full_taxon_dict[name]
                    if "sample_date" in taxon_obj.attribute_dict.keys():
                        date = taxon_obj.attribute_dict["sample_date"]
                        k.traits["display"] = f"{name}|{date}"
                    else:
                        k.traits["display"] = name
                else:
                    #taxon_obj = ""
                    k.traits["display"] = name + "|" + "not in dict"

                
                
                # if "county" in taxon_obj.attribute_dict.keys(): 
                #     county = taxon_obj.attribute_dict["county"]
                # elif "adm2" in taxon_obj.attribute_dict.keys():
                #     county = taxon_obj.attribute_dict["adm2"]
                # else:
                #     county = country
                
            #     if "country" in k.traits: 
            #             k.traits["display"]= f"{name}|{date}"
            #         else:
            #             k.traits["display"]=""
            #     else:
            #         k.traits["display"]=""

            # else:
def find_colour_dict(query_dict, trait):

    attribute_options = set()

    if trait == "adm1":
        colour_dict = {"Wales":"darkseagreen",
                "England":"indianred",
                "Scotland":"steelblue",
                "Northern_Ireland":"skyblue"}
        return colour_dict

    else:
        for query in query_dict.values():
            attribute_options.add(query.attribute_dict[trait])
            
    if len(attribute_options) == 2:
        colour_dict = {list(attribute_options)[0]: "goldenrod",
                        list(attribute_options)[1]:"midnightblue"}
        return colour_dict
    
    
def make_scaled_tree_without_legend(My_Tree, tree_name, tree_dir, num_tips, colour_dict, trait, tallest_height,lineage, taxon_dict, query_id_dict, query_dict, tree_to_query):

    display_name(My_Tree, tree_name, tree_dir, query_id_dict, taxon_dict) #this is id dict for when the ids are in the tree.
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
    space_offset = tallest_height/100
    absolute_x_axis_size = tallest_height+space_offset+space_offset + tallest_height #changed from /3 
    
    tipsize = 40
    c_func=lambda k: 'dimgrey' ## colour of branches
    l_func=lambda k: 'lightgrey' ## colour of branches
    cn_func = lambda k: 'fuchsia' if k.name in query_id_dict.keys() or k.name in query_dict.keys() else 'dimgrey'
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

            if k.name in query_dict.keys() or k.name in query_id_dict.keys():
                tree_to_query[lineage].append(k.name)

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

    return tree_to_query

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

def make_all_of_the_trees(input_dir, taxon_dict, query_id_dict, query_dict, tree_to_query, desired_fields, min_uk_taxa=3):
    
    tallest_height = find_tallest_tree(input_dir)

    too_tall_trees = []

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
            tips = []
            
            for k in tree.Objects:
                if k.branchType == 'leaf':
                    tips.append(k.name)
            if len(tips) < 1000:
                overall_tree_count += 1      
                if desired_fields == []:
                    colour_by = ["adm1"]
                else:
                    colour_by = desired_fields
            
                for trait in colour_by:
                    colour_dict = find_colour_dict(query_dict, trait)
                    tree_to_query = make_scaled_tree_without_legend(tree, treename, input_dir, len(tips), colour_dict, trait, tallest_height, lineage, taxon_dict, query_id_dict, query_dict, tree_to_query)     
            else:
                too_tall_trees.append(lineage)
                continue

    return too_tall_trees, overall_tree_count, tree_to_query

def summarise_collapsed_node(tree_dir, focal_node, focal_tree, full_tax_dict):

    focal_tree_file = focal_tree + ".txt"

    with open(tree_dir + "/" + focal_tree_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split("\t")
            node_name = toks[0]
            members = toks[1]
        
            if node_name == focal_node:
                dates = []
                countries = []
                
                member_list = members.split(",")
                number_nodes = str(len(member_list)) + " nodes"

                for tax in member_list:
                    if tax in full_tax_dict.keys():
                        taxon_obj = full_tax_dict[tax]
                    
                        if taxon_obj.sample_date != "NA":
                            date_string = taxon_obj.sample_date
                            date = dt.datetime.strptime(date_string, "%Y-%m-%d").date()
                            dates.append(date)
                        
                        countries.append(taxon_obj.attribute_dict["country"])
                    
                    else:
                        country = tax.split("/")[0]
                        countries.append(country)

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


                #taken out for now until we can get dates for all the tips
                # min_date = str(min(dates))
                # max_date = str(max(dates))
                #info = number_nodes + ", ranging from " + min_date + " to " + max_date + " in " + pretty_countries

                info = number_nodes + " in " + pretty_countries



    return info

