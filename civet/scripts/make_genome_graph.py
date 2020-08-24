#!/usr/bin/env python3
import os
import argparse
import matplotlib as mpl
from matplotlib import pyplot as plt
import csv
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from itertools import cycle

colour_list = ["#cc8180","#a0a997","#cbaca4","#cadfbc"]
colour_cycle = cycle(colour_list)

def parse_args():
    parser = argparse.ArgumentParser(description='Find sequences relative to Wuhan4 reference.')

    parser.add_argument("--input", action="store", type=str, dest="input")
    parser.add_argument("--output", action="store", type=str, dest="output")

    return parser.parse_args()

def next_colour():
    return next(colour_cycle)


def make_graph():

    args = parse_args()

    y_position = 0
    ref_vars = {}
    fig, ax = plt.subplots(1,1, figsize=(10,2.5), dpi=250)

    

    with open(args.input,newline="") as f:
        reader=csv.DictReader(f, delimiter='\t')
        for row in reader:
            y_position +=1
            snps = row["snps"].split(";")
            x = []
            y = []
            col = next_colour()
            rect = patches.Rectangle((0,y_position-0.5), 29903, 1 ,alpha=0.3, fill=True, edgecolor='none',facecolor=col)
            ax.add_patch(rect)
            
            for snp in snps:
                x_position = int(snp[:-2])
                base = snp[-1]
                ref = snp[-2]
                ref_vars[x_position]=ref
                
                ax.text(x_position, y_position, base, size=8, ha="center", va="center")
            
            ax.text(-20, y_position, row["name"], size=7, ha="right", va="center")
            
    rect = patches.Rectangle((0,-1), 29903, 1 ,alpha=0.2, fill=True, edgecolor='none',facecolor="dimgrey")
    ax.add_patch(rect)

    ax.text(-20, -0.5, "Reference", size=7, ha="right", va="center")

    for var in ref_vars:
    #     ax.scatter([var],[0], color="#924242", s = 9, marker="v")
        ax.text(var, -0.5, ref_vars[var], size=8, ha="center", va="center") 
        

    ax.spines['top'].set_visible(False) ## make axes invisible
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.yticks([])
            
    ax.set_xlim(0,29903)
    ax.set_ylim(-1,y_position+1)
    ax.tick_params(axis='x', labelsize=8)
    plt.xlabel("Genome position (base)", fontsize=9)
    plt.tight_layout()
    plt.savefig(args.output)

if __name__ == '__main__':

    make_graph()