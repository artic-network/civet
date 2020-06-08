import os
from matplotlib import font_manager as fm, rcParams
import baltic as bt

# from IPython.display import HTML
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
import matplotlib.font_manager as font_manager


thisdir = os.path.abspath(os.path.dirname(__file__))
font_file = os.path.join(thisdir, "..","data/HelveticaNeue.ttf")
print(font_file)
font_list = font_manager.fontManager.addfont(font_file)

mpl.rcParams['font.family'] = 'helveticaneue'
mpl.rcParams['font.weight']=300
mpl.rcParams['axes.labelweight']=300


def find_tallest_tree(input_dir):
    tree_heights = []
    
    for r,d,f in os.walk(input_dir):
        for fn in f:
            if fn.endswith("tree"):
                treename = fn
                num_taxa = 0
                intro_name = ""
                with open(r + '/' + fn,"r") as f:
                    for l in f:
                        l = l.rstrip("\n")

                        if l.startswith("    DIMENSIONS NTAX"):

                            num_taxa = int(l.rstrip(";").split("=")[1])
                            intro_name = fn.rstrip(".tree")
                if num_taxa > 1:
                    tree = bt.loadNewick(r + '/' + fn,absoluteTime=False)
                    tips = []
                    uk_tips = []
                    
                    for k in tree.Objects:
                        if k.branchType == 'leaf':
                            if k.traits["country_uk"]=="True": #this may change
                                uk_tips.append(k.name)
                            tips.append(k.name)
                    if len(uk_tips)>3:
                        tree_heights.append(tree.treeHeight)
    
    max_height = sorted(tree_heights, reverse=True)[0]

    return max_height

def display_name(tree, taxon_dict):
    for k in tree.Objects:
        if k.branchType == 'leaf':
            name = k.name
            
            if name in taxon_dict.keys():
                taxon_obj = taxon_dict[name]
                date = taxon_obj.sample_date
                country = taxon_obj.adm1
                
                # if "county" in taxon_obj.attribute_dict.keys(): 
                #     county = taxon_obj.attribute_dict["county"]
                # elif "adm2" in taxon_obj.attribute_dict.keys():
                #     county = taxon_obj.attribute_dict["adm2"]
                # else:
                #     county = country
                
                if "country" in k.traits: 
                    if k.traits["country"] in ["Wales","Scotland","England","Northern_Ireland"]:
                        #county = county.capitalize()
                        k.traits["display"]= f"{name}|{date}"
                    else:
                        k.traits["display"]=""
                else:
                    k.traits["display"]=""

            else:
                k.traits["display"]=name

def annotate_country(tree, taxon_dict):
    for k in tree.Objects:
        if k.branchType == 'leaf':
            name = k.name
            
            country = name.split("/")[0]

            # tax_object = taxon_dict[name]
            # country = tax_object.adm1
            
            if k.traits["country_uk"] == "True":
                k.traits["country"]=country
                if country == "WALES":
                    k.traits["country"]="Wales"
                else:
                    k.traits["country"]=country
        else:
            k.traits["country"]= ""

def make_scaled_tree_without_legend(My_Tree, num_tips, colour_dict, tallest_height,lineage, taxon_dict, query_dict):
    
    annotate_country(My_Tree, taxon_dict)
    display_name(My_Tree, query_dict)
    My_Tree.uncollapseSubtree()

    if num_tips < 10:
        #page_height = num_tips/2
        page_height = num_tips
    else:
        #page_height = num_tips/4 
        page_height = num_tips/2  

    offset = tallest_height - My_Tree.treeHeight
    space_offset = tallest_height/100
    absolute_x_axis_size = tallest_height+space_offset+space_offset + tallest_height #changed from /3 
    
    tipsize = 20
    c_func=lambda k: 'dimgrey' ## colour of branches
    l_func=lambda k: 'lightgrey' ## colour of branches
    # cn_func=lambda k: colour_dict[k.traits["country"]] if k.traits["country"] in colour_dict else 'dimgrey'
    cn_func = lambda k: 'fuchsia' if k.name in query_dict.keys() else 'dimgrey'
    #s_func=lambda k: tipsize*5 if k.traits["country"] in colour_dict else tipsize
    s_func = lambda k: tipsize*5 if k.name in query_dict.keys()  else tipsize
    z_func=lambda k: 100
    b_func=lambda k: 0.5 #branch width
    #co_func=lambda k: colour_dict[k.traits["country"]] if k.traits["country"] in colour_dict else 'dimgrey' ## for plotting a black outline of tip circles
    co_func=lambda k: 'fuchsia' if k.name in query_dict.keys() else 'dimgrey' 
    #so_func=lambda k: tipsize*5 if k.traits["country"] in colour_dict else 0 #plots the uk tips over the grey ones
    so_func=lambda k: tipsize*5 if k.name in query_dict.keys() else 0
    zo_func=lambda k: 99
    # outline_func = lambda k: None
    #outline_colour_func = lambda k: colour_dict[k.traits["country"]] if k.traits["country"] in colour_dict else 'dimgrey'
    outline_colour_func = lambda k: "fuchsia" if k.name in query_dict.keys() else 'dimgrey'
    zb_func=lambda k: 98
    zt_func=lambda k: 97
    font_size_func = lambda k: 25 if k.name in query_dict.keys() else 15
    kwargs={'ha':'left','va':'center','size':12}
    
    x_attr=lambda k: k.height + offset
    y_attr=lambda k: k.y

    y_values = []
    for k in My_Tree.Objects:
        y_values.append(y_attr(k))
    y_values = sorted(y_values)
    vertical_spacer = 0.5 
    full_page = page_height + vertical_spacer + vertical_spacer
    min_y,max_y = y_values[0]-vertical_spacer,y_values[-1]+vertical_spacer
    
    fig2,ax2 = plt.subplots(figsize=(20,page_height),facecolor='w',frameon=False)

    My_Tree.plotTree(ax2, colour_function=c_func, x_attr=x_attr, y_attr=y_attr, branchWidth=b_func)
    My_Tree.plotPoints(ax2, x_attr=x_attr, colour_function=cn_func,y_attr=y_attr, size_function=s_func, outline_colour=outline_colour_func)
    My_Tree.plotPoints(ax2, x_attr=x_attr, colour_function=co_func, y_attr=y_attr, size_function=so_func, outline_colour=outline_colour_func)

    for k in My_Tree.Objects:
        # if "country" in k.traits: #don't need these for now because we're only doing UK trees
        #     if k.traits["country"] in colour_dict:
        if "display" in k.traits:
            name=k.traits["display"] 

            x=x_attr(k)
            y=y_attr(k)
        
            height = My_Tree.treeHeight+offset
            
            ax2.text(tallest_height+space_offset+space_offset, y, name, size=font_size_func(k), ha="left", va="center", fontweight="bold")
            ax2.plot([x+space_offset,tallest_height+space_offset],[y,y],ls='--',lw=0.5,color=l_func(k))

    ax2.spines['top'].set_visible(False) ## make axes invisible
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    
    ax2.set_xlim(-space_offset,absolute_x_axis_size)
    ax2.set_ylim(min_y,max_y)

    plt.yticks([])
    plt.xticks([])
    plt.title(lineage, size=25)

    fig2.tight_layout()

def sort_trees_index(tree_dir):
    b_list = []
    d_list = []
    for r,d,f in os.walk(tree_dir):
        for thing in f:
            if "tree" in thing:
                a = thing.split(".")[0]
                b = a.split("_")[2]
                b2 = b.strip("UK")
                b_list.append(b2)
        
    c = sorted(b_list, key=int)
    for i in c:
        d = "UK" + i
        d_list.append(d)
        
    return d_list

    
 
def make_all_of_the_trees(input_dir, taxon_dict, query_dict, min_uk_taxa=3):
    
    tallest_height = find_tallest_tree(input_dir)
    
    # colour_dict = {"Wales":"darkseagreen",
    #            "England":"indianred",
    #            "Scotland":"steelblue",
    #            "Northern_Ireland":"skyblue"}

    lst = sort_trees_index(input_dir)

    for fn in lst:
        lineage = fn
        treename = "uk_lineage_" + lineage + ".tree"
        
        num_taxa = 0
        intro_name = ""
        with open(input_dir + "/" + treename,"r") as f:
            for l in f:
                l = l.rstrip("\n")
                if l.startswith("    DIMENSIONS NTAX"):
                    num_taxa = int(l.rstrip(";").split("=")[1])
                    intro_name = fn.rstrip(".tree")

        if num_taxa > 1: 
            tree = bt.loadNewick(input_dir + "/" + treename, absoluteTime=False)
            tips = []
            uk_tips = []
            for k in tree.Objects:
                if k.branchType == 'leaf':
                    if k.traits["country_uk"]=="True":
                        uk_tips.append(k.name)
                    tips.append(k.name)
            if len(uk_tips)>min_uk_taxa:
                #outfile = f'{outdir}/{treename}.pdf'
                #make_scaled_tree_without_legend(tree, outfile,len(tips),colour_dict, tallest_height)
                try:
                    make_scaled_tree_without_legend(tree,len(tips),colour_dict, tallest_height, lineage, taxon_dict, query_dict)     
                except ValueError:
                    pass
    return lst

def sort_fig(fig_dir):
    b_list = []
    d_list = []
    for fig_file in os.listdir(fig_dir):
        if "trees" in fig_file:
            a = fig_file.split(".")
            b = a[0].split("_")
            b_list.append(b[-1])
            stem = b[:len(b)-1]
            stem = "_".join(b[:len(b)-1])

    c = sorted(b_list, key=int)
    
    for i in c:
        d = stem + "_" + i + ".png"
        d_list.append(d)

    return d_list
                
