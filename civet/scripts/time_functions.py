import datetime as dt 
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.ticker as plticker
import math


def summarise_dates(query_dict):

    overall_dates = []

    

    for tax in query_dict.values():
        overall_dates.extend(list(tax.date_dict.values()))

    if overall_dates != []:

        min_overall_date = min(overall_dates)
        max_overall_date = max(overall_dates)

        min_string = min_overall_date.strftime("%Y-%m-%d")
        max_string = max_overall_date.strftime("%Y-%m-%d")

        dates_present = True

    else:
        dates_present = False


    return dates_present, overall_dates, max_overall_date, min_overall_date,  max_string, min_string

def display_name(tax, custom_tip_fields):
    
    name = tax.name
            
    date = tax.sample_date
    display_name = f"{name}|{date}"
    
    if "adm2" in tax.attribute_dict.keys():
        adm2 = tax.attribute_dict["adm2"]
        display_name = f"{name}|{adm2}|{date}"

    if len(custom_tip_fields) > 0: 
        for label_element in custom_tip_fields:
            display_name = display_name + "|" + tax.attribute_dict[label_element]
    
    return display_name

def find_colour_dict(tree_list):
    
    colour_dict = {}
    count = 0

    cmap = cm.get_cmap("Paired")
    
    colors = cmap(np.linspace(0, 1, len(tree_list)))
    
    for option in sorted(tree_list):
        colour_dict[option] = colors[count]
        count += 1

    return colour_dict

def plot_time_series(query_dict, overall_max_date, overall_min_date, tree_list, custom_tip_fields):

    time_len = (overall_max_date - overall_min_date).days

    tick_loc_base = float(math.ceil(time_len/5))
    
    loc = plticker.MultipleLocator(base=tick_loc_base)

    fig, ax1 = plt.subplots(1,1, figsize=(20,15))
    ax2 = ax1.twinx()
    
    if time_len > 10:
        offset = dt.timedelta(time_len/10)
    else:
        offset = dt.timedelta(time_len/3)

    colour_dict = find_colour_dict(tree_list)

    count = 0

    for tax in query_dict.values():
        if tax.date_dict != {}:
        
            first_date = min(list(tax.date_dict.values()))
            last_date = max(list(tax.date_dict.values()))
            
            label = display_name(tax, custom_tip_fields)
            
            x = [first_date, last_date]
            y = [count, count]
            
            ax1.scatter(first_date, count, color=colour_dict[tax.tree], s=200, zorder=2)
            ax1.scatter(last_date, count, color=colour_dict[tax.tree], s=200, zorder=2)
            
            if first_date != last_date:
                ax1.plot(x,y, color=colour_dict[tax.tree],zorder=1)

            # ax2.plot([last_date, overall_max_date+offset],y,ls='dotted',lw=1,color="dimgrey")
            if x != overall_max_date:
                ax1.plot([last_date, overall_max_date+offset],y,ls='--',lw=1,color="dimgrey", zorder=1)

            ax2.text(overall_max_date+offset, count, label,size=15)
            
            count += 1
        
    ylim = ax1.get_ylim()
    ax2.set_ylim(ylim)
        
    ax1.spines['top'].set_visible(False) ## make axes invisible
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_yticks([])
    ax2.spines['top'].set_visible(False) ## make axes invisible
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_yticks([])

    ax1.tick_params(labelsize=20, rotation=90)
    ax1.xaxis.set_major_locator(loc)


    fig.tight_layout()






