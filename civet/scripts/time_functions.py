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

    else:
        return False


    return overall_dates, max_overall_date, min_overall_date,  max_string, min_string

def display_name(tax, custom_tip_fields):
    
    name = tax.name
            
    date = tax.sample_date
    display_name = f"{name}|{date}"
    
    if "adm2" in tax.attribute_dict.keys():
        adm2 = tax.attribute_dict["adm2"]
        display_name = f"{name}|{adm2}|{date}"

    count = 0
    if len(custom_tip_fields) > 0: 
        for label_element in custom_tip_fields:
            if count == 0:
                display_name = tax.attribute_dict[label_element] 
            else:   
                display_name = display_name + "|" + tax.attribute_dict[label_element]
            count += 1
    
    return display_name

def find_colour_dict(date_fields):
    
    colour_dict = {}
    count = 0

    cmap = cm.get_cmap("Paired")
    
    colors = cmap(np.linspace(0, 1, len(date_fields)))
    
    for option in sorted(date_fields):
        colour_dict[option] = colors[count]
        count += 1

    return colour_dict

def plot_time_series(tips, query_dict, overall_max_date, overall_min_date, date_fields, custom_tip_fields):

    colour_dict = find_colour_dict(date_fields)    

    time_len = (overall_max_date - overall_min_date).days

    height = math.sqrt(len(tips))*2 + 1

    if time_len > 20:
        tick_loc_base = float(math.ceil(time_len/5))
    else:
        tick_loc_base = 1.0

    loc = plticker.MultipleLocator(base=tick_loc_base) #Sets a tick on each integer multiple of a base within the view interval

    fig, ax1 = plt.subplots(1,1, figsize=(20,height))
    ax2 = ax1.twinx()
    
    if time_len > 10:
        offset = dt.timedelta(time_len/10)
    else:
        offset = dt.timedelta(time_len/3)


    count = 1

    for tax in tips:
        if tax.date_dict != {} and tax in query_dict.values():
        
            first_date_type = min(tax.date_dict.keys(), key=lambda k: tax.date_dict[k])
            last_date_type = max(tax.date_dict.keys(), key=lambda k: tax.date_dict[k])

            first_date = tax.date_dict[first_date_type]
            last_date = tax.date_dict[last_date_type]

            other_dates = {}
            for date_type, date in tax.date_dict.items():
                if date != first_date and date != last_date:
                    other_dates[date_type] = date
            
            label = display_name(tax, custom_tip_fields)
            
            x = [first_date, last_date]
            y = [count, count]
            
            ax1.scatter(first_date, count, color=colour_dict[first_date_type], s=200, zorder=2, label=first_date_type)
            ax1.scatter(last_date, count, color=colour_dict[last_date_type], s=200, zorder=2, label=last_date_type)

            for date_option, date in other_dates.items():
                ax1.scatter(date, count, color=colour_dict[date_option], s=200, zorder=2, label=date_option)
            
            if first_date != last_date:
                ax1.plot(x,y, color="dimgrey",zorder=1)

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

    handles, labels = ax1.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax1.legend(*zip(*unique))

    fig.tight_layout()






