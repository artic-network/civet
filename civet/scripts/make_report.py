#!/usr/bin/env python3
import os
from pweave import weave
import argparse
import shutil
import sys


thisdir = os.path.abspath(os.path.dirname(__file__))

def get_report_arguments(arg_file):

    args = []

    with open(arg_file) as f:
        for l in f:
            arg = l.strip("\n")
            arg = arg.replace("-", "_")
            args.append(arg)

    return args

def make_report(inputs, report_args_file):
    
    arguments = get_report_arguments(report_args_file)

    arg_dict = {}
    for name, value in zip(arguments,inputs):

        if name == "outfile":
            outfile = value
        elif name == "report_template":  #might be able to get rid of this, just have one common template that we use different configs for
            md_template = value
        else:
            if name == "cog_metadata":
                name = "full_metadata_file"
            else:
                name = name
            
            arg_dict[name] = value

    name_stem = ".".join(outfile.split(".")[:-1])
    arg_dict["name_stem_input"] = name_stem
                        
    tree_name_stems = []
    for r,d,f in os.walk(arg_dict["tree_dir"]):
        for fn in f:
            if fn.endswith(".tree"):
                basename = ".".join(fn.split(".")[:-1])
                stem = "_".join(basename.split("_")[:-1])
                tree_name_stems.append(stem)
    tree_name_stems = list(set(tree_name_stems))
    
    if len(tree_name_stems) > 1:
        sys.stderr.write("Error: Multiple tree names found")
        sys.exit(-1)
    elif len(tree_name_stems) == 0:
        sys.stderr.write("Error: No trees found in tree directory. Note, tree name much end with .tree")
        sys.exit(-1)
    else:
        tree_name_stem = tree_name_stems[0]

    arg_dict["tree_name_stem"] = tree_name_stem

    with open(outfile, 'w') as pmd_file:
    
        summary_dir = os.path.join(arg_dict["outdir"], "summary_files")

        change_line_dict = {}

        for key, value in arg_dict.items():
            new_key = key
            new_value = f'{key} = "{value}"\n'
            change_line_dict[key] = new_value
        
        if arg_dict["add_bars"]:
            change_line_dict["add_bars"] = f'add_bars = "{add_bars}"\n'
        else:
            change_line_dict["add_bars"] = 'add_bars = ""'
        
        with open(md_template) as f:
            for l in f:
                if "##CHANGE" in l:
                    for key in change_line_dict:
                        if key in l:
                            new_l = change_line_dict[key]
                else:
                    new_l = l

                pmd_file.write(new_l)
    
    weave(outfile, doctype = "pandoc", figdir=arg_dict["figdir"])

def main():
    parser = argparse.ArgumentParser(description="Report generator script")
    
    parser.add_argument("--report-args", required=True, help="report arguments text file", dest="report_args")
    
    parser.add_argument("-i","--input-csv", required=False, help="path to input file",dest="input_csv")
    parser.add_argument("--filtered-cog-metadata", required=False, help="path to combined metadata file",dest="filtered_cog_metadata")
    parser.add_argument("--cog-metadata", required=True, help="path to full COG metadata file",dest="cog_metadata")

    parser.add_argument("-o","--outfile", default="civet_report.pmd", help="output name stem as a string",dest="outfile")
    parser.add_argument("--report-template", help="report template file",dest="report_template")
    parser.add_argument("--outdir", help="output directory",dest="outdir")
    parser.add_argument("-t","--treedir", required=False, default="", help="path to tree directory",dest="treedir")
    parser.add_argument("--figdir", help="output directory",dest="figdir")

    
    parser.add_argument("-f", "--tree-fields-input",default="", help="desired fields to colour trees by in report. Default=UK country",dest="tree_fields_input")
    parser.add_argument("-gd", "--graphic-dict-input", default="", help="fields to colour by rather than display text. Add colour scheme optionally", dest="graphic_dict_input")
    parser.add_argument("-l", "--label-fields-input", default="", help="fields to add into labels in report trees. Default is adm2 and date", dest='label_fields_input')
    parser.add_argument("--date-fields-input", help="column headers containing date information as a a comma separated string.", dest="date_fields_input")
    parser.add_argument("--node-summary-option", action="store", help="field to summarise collapsed nodes by. Default=lineage", dest="node_summary_option")

    parser.add_argument("-sc", "--seq-centre",default="", help="Sequencing centre", dest="seq_centre")
    parser.add_argument("--failed-seqs", required=False, default="", help="csv of seqs that fail qc and the reason why",dest="failed_seqs")  
    parser.add_argument("--add-bars", action="store_true",dest="add_bars",default=False)  


    parser.add_argument("--clean-locs-file", required=True, help="CSV for cleaning adm2 regions in metadata", dest="clean_locs_file")
    parser.add_argument("--pc-file", required=True, help="file containing outer postcode to centroid mapping", dest="pc_file")

    parser.add_argument("--uk-map", required=True, help="shape file for uk counties", dest="uk_map")
    parser.add_argument("--channels-map", required=True, help="shape file for channel islands", dest="channels_map")
    parser.add_argument("--ni-map", required=True, help="shape file for northern irish counties", dest="ni_map")
    parser.add_argument("--urban-centres", default="", help="geojson for plotting urban centres", dest="urban_centres")

    parser.add_argument("--local-lineages", default="", action='store_true',help="List of rendered .png files for local lineage analysis", dest="local_lineages")
    parser.add_argument("--local-lin-maps", default="", action='store',help="List of rendered .png files for local lineage analysis", dest="local_lin_maps")
    parser.add_argument("--local-lin-tables", default="", action='store', help="List of .md tables for local lineage analysis", dest="local_lin_tables")

    parser.add_argument("--map-sequences", required=True, help="Bool for whether mapping of sequences by trait is required", dest="map_sequences")
    parser.add_argument("--map-cols", default="", help="either column names in input csv which contains x coords and y coords for mapping as a comma separated string OR column name containing outer postcode", dest="map_inputs")
    parser.add_argument("--input-crs", default="", help="coordinate reference system that x and y inputs are in", dest="input_crs")
    parser.add_argument("--mapping-trait", default="", help="trait to map sequences by", dest="mapping_trait")

    
    args = parser.parse_args()

    arg_list = get_arg_list(args)
    make_report(arg_list, args.report_args)


def get_arg_list(args):

    arg_list = [
        args.cog_metadata, args.input_csv, args.filtered_cog_metadata, #metadata inputs
        args.outfile, args.report_template, args.outdir, args.tree_dir, args.figdir, #directories and file sorting
        args.tree_fields_input, args.graphic_dict_input, args.label_fields_input, args.date_fields_input, args.node_summary_option, #display options for tree and timeline
        args.failed_seqs, args.seq_centre, args.add_bars, #misc options
        args.clean_locs_file, args.pc_file, args.uk_map, args.channels_map, args.ni_map, args.urban_centres #mapping files
        args.local_lineages, args.local_lin_maps, args.local_lin_tables, #background lineage mapping
        args.map_sequences, args.map_inputs, args.input_crs, args.mapping_trait] #options for mapping sequences in query

    return arg_list

if __name__ == "__main__":
    main()
