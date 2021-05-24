import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt

## No uk-specific stuff in here if at all possible - the local lineages should just compare to whatever column they give
## Should just take any map file - should have a check that some of the sequences in match to 


##check that shape file and column to plot for local lineages is there

##not sure where the actual jsons should be made

def parse_map_options(metadata, map_queries, map_background, uk, query_map_column, query_map_colour, background_map_column,background_map_date_window, background_map_date_start,background_map_date_end, config):

    


 if args.map_queries or args.map_background:
        if args.uk: #or, this could be if it's run on CLIMB or some other proxy for it being UK
            from civet.output_options import uk_maps
            uk_maps.uk_specific_map_actions()

def map_checks()

def test_date_restrictions():



def do_date_restrictions():



def check_shape_file():


def check__map_col(): #presence of arg and also in the csv metadata
    #this can be general for any type of map, run it independently for local lineages or mapping the dots
    # if it's UK, then this will be provided already so the check can be the same

def check_presence_shapefile():



