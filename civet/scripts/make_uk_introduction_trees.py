import os
from matplotlib import font_manager as fm, rcParams
import utils.baltic as bt

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

font_dirs = ['utils/helveticaneue', ]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
font_list = font_manager.createFontList(font_files)
font_manager.fontManager.ttflist.extend(font_list)

mpl.rcParams['font.family'] = 'helveticaneue'
mpl.rcParams['font.weight']=300
mpl.rcParams['axes.labelweight']=300
mpl.rcParams['font.size']=20


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
            taxon_obj = taxon_dict[name]
            date = taxon_obj.sample_date
            country = taxon_obj.adm1
            
            if "county" in taxon_obj.attribute_dict.keys(): 
                county = taxon_obj.attribute_dict["county"]
            elif "adm2" in taxon_obj.attribute_dict.keys():
                county = taxon_obj.attribute_dict["adm2"]

            if county == "":
                county = country
            
            if "country" in k.traits: 
                if k.traits["country"] in ["Wales","Scotland","England","Northern_Ireland"]:
                    county = county.capitalize()
                    k.traits["display"]= f"{name}|{county}|{date}"
                else:
                    k.traits["display"]=""
            else:
                k.traits["display"]=""

def annotate_country(tree, taxon_dict):
    for k in tree.Objects:
        if k.branchType == 'leaf':
            name = k.name
            tax_object = taxon_dict[name]
            country = tax_object.adm1
            
            if k.traits["country_uk"] == "True":
                k.traits["country"]=country
                if country == "WALES":
                    k.traits["country"]="Wales"
                else:
                    k.traits["country"]=country
        else:
            k.traits["country"]= ""

def make_scaled_tree_without_legend(My_Tree, num_tips, colour_dict, tallest_height,intro, taxon_dict):
    
    annotate_country(My_Tree)
    display_name(My_Tree, taxon_dict)
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
    cn_func=lambda k: colour_dict[k.traits["country"]] if k.traits["country"] in colour_dict else 'dimgrey'
    s_func=lambda k: tipsize*5 if k.traits["country"] in colour_dict else tipsize
    z_func=lambda k: 100
    b_func=lambda k: 0.5 #branch width
    co_func=lambda k: colour_dict[k.traits["country"]] if k.traits["country"] in colour_dict else 'dimgrey' ## for plotting a black outline of tip circles
    so_func=lambda k: tipsize*5 if k.traits["country"] in colour_dict else 0 #plots the uk tips over the grey ones
    zo_func=lambda k: 99
    # outline_func = lambda k: None
    outline_colour_func = lambda k: colour_dict[k.traits["country"]] if k.traits["country"] in colour_dict else 'dimgrey'
    zb_func=lambda k: 98
    zt_func=lambda k: 97
    font_size = 25
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
        if "country" in k.traits:
            if k.traits["country"] in colour_dict:
                name=k.traits["display"] ## get travelling segment(s)
            
                x=x_attr(k)
                y=y_attr(k)
            
                height = My_Tree.treeHeight+offset
                
                ax2.text(tallest_height+space_offset+space_offset, y, name, size=font_size, ha="left", va="center", fontweight="bold")
                ax2.plot([x+space_offset,tallest_height+space_offset],[y,y],ls='--',lw=0.5,color=l_func(k))

    ax2.spines['top'].set_visible(False) ## make axes invisible
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    
    ax2.set_xlim(-space_offset,absolute_x_axis_size)
    ax2.set_ylim(min_y,max_y)

    plt.yticks([])
    plt.xticks([])
    plt.title(intro)

    fig2.tight_layout()

def sort_trees_index(tree_dir):
    b_list = []
    d_list = []
    for r,d,f in os.walk(tree_dir):
        for thing in f:
            if "tree" in thing:
                a = thing.split(".")
                b = a[0].strip("UK")
                b_list.append(b)
        
    c = sorted(b_list, key=int)
    for i in c:
        d = "UK" + i
        d_list.append(d)
        
    return d_list


def make_pie(intro, df_counts):

    fig,ax1 = plt.subplots(figsize=(2,2), subplot_kw=dict(aspect="equal"), dpi=100)

    pie_data = list(df_counts[intro])
    pie_labels='England','Scotland','Wales','Nothern_Ireland'
    
    colour_dict = {"wales":"darkseagreen",
            "england":"indianred",
            "scotland":"steelblue",
            "northern_ireland":"skyblue"}

    pie_colours = []
    for i in df_counts.index:
        pie_colours.append(colour_dict[i])
        
    total = sum(pie_data)
    wedges,texts,autotexts = ax1.pie(pie_data,colors=pie_colours, autopct=lambda p:'{:.0f}'.format(p*total/100) if p > 0 else '', shadow=False, startangle=140)

    plt.setp(autotexts, size=10)
    plt.title(intro)

    #plt.show()
    
 
def make_all_of_the_trees(input_dir, df_counts, taxon_dict, min_uk_taxa=3):
    #outdir = outdir.rstrip("/")
    intro_list = []
    
    tallest_height = find_tallest_tree(input_dir)
    
    colour_dict = {"Wales":"darkseagreen",
               "England":"indianred",
               "Scotland":"steelblue",
               "Northern_Ireland":"skyblue"}

    lst = sort_trees_index(input_dir)

    for fn in lst:
        lineage = fn
        treename = "uk_lineage+" + lineage + ".tree"
        
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
                    make_scaled_tree_without_legend(tree,len(tips),colour_dict, tallest_height, lineage, taxon_dict)     
                except ValueError:
                    pass




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
                
