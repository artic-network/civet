#!/usr/bin/env python3
import os
from pweave import weave
import argparse
import shutil
import sys


thisdir = os.path.abspath(os.path.dirname(__file__))

def get_report_arguments1(arg_file):

    args = []

    with open(arg_file) as f:
        for l in f:
            arg = l.strip("\n")
            arg = arg.replace("-", "_")
            args.append(arg)

    return args

def get_report_arguments2(arg_file):

    args = [] #find and replace stuff as well

    with open(arg_file) as f:
        for l in f:
            arg = l.strip("\n")
            arg = arg.replace("-", "_")
            arg2 = "arg." + arg
            args.append(arg2)

    return args


# def make_report(cog_metadata, input_csv, filtered_cog_metadata, outfile, outdir, treedir, figdir, colour_fields, label_fields, node_summary, report_template, failed_seqs, seq_centre, clean_locs, uk_map, channels_map, ni_map, pc_file, local_lineages, local_lin_maps, local_lin_tables,map_sequences,map_inputs, input_crs,mapping_trait,urban_centres,add_bars, graphic_dict, date_fields):
def make_report(inputs, report_args_file):
    
    arguments = get_report_arguments1(report_args_file)

    arg_dict = {}
    for name, value in zip(arguments,inputs):
        if name == "cog_metadata":
            name = "full_metadata_file"
        else:
            name = name
        arg_dict[name] = value

    outfile = arg_dict["outfile"]
    outdir = arg_dict["outdir"]
    report_template = arg_dict["report_template"]
    figdir = arg_dict["figdir"]

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
    
        md_template = report_template
        summary_dir = os.path.join(outdir, "summary_files")

        change_line_dict = {}

        for key, value in arg_dict.items():
            new_key = key
            new_value = f'{key} = "{value}"\n'
            change_line_dict[key] = new_value

        # change_line_dict = {
        #                     "output_directory": f'output_directory = "{arg_dict[outdir]}"\n',
        #                     "name_stem_input": f'name_stem_input = "{arg_dict[name_stem]}"\n',
        #                     "full_metadata_file": f'full_metadata_file = "{cog_metadata}"\n',
        #                     "filtered_cog_metadata": f'filtered_cog_metadata = "{filtered_cog_metadata}"\n',
        #                     "input_csv": f'input_csv = "{input_csv}"\n',
        #                     "tree_fields_input": f'tree_fields_input = "{tree_fields_input}"\n',
        #                     "label_fields_input": f'label_fields_input = "{label_fields_input}"\n',
        #                     "date_fields_input":f'date_fields_input = "{date_fields_input}"\n',
        #                     "node_summary_option": f'node_summary_option = "{node_summary_option}"\n',
        #                     "figdir": f'figdir = "{figdir}"\n',
        #                     "tree_dir": f'tree_dir = "{tree_dir}"\n',
        #                     "tree_name_stem": f'tree_name_stem = "{tree_name_stem}"\n',
        #                     "summary_dir": f'summary_dir = "{summary_dir}"\n',
        #                     "failed_seqs": f'failed_seqs = "{failed_seqs}"\n',
        #                     "seq_centre": f'seq_centre = "{seq_centre}"\n',
        #                     "clean_locs_file": f'clean_locs_file = "{clean_locs_file}"\n',
        #                     "uk_map": f'uk_map = "{uk_map}"\n',
        #                     "channels_map": f'channels_map = "{channels_map}"\n',
        #                     "ni_map": f'ni_map = "{ni_map}"\n',
        #                     "pc_file": f'pc_file = "{pc_file}"\n',
        #                     "local_lineages":f'local_lineages = "{local_lineages}"\n',
        #                     "local_lin_maps" : f'local_lin_maps = "{local_lin_maps}"\n',
        #                     "local_lin_tables" : f'local_lin_tables = "{local_lin_tables}"\n',
        #                     "map_sequences":f'map_sequences = "{map_sequences}"\n', 
        #                     "map_inputs":f'map_inputs = "{map_inputs}"\n',
        #                     "mapping_trait":f'mapping_trait = "{mapping_trait}"\n',
        #                     "input_crs":f'input_crs = "{input_crs}"\n',
        #                     "urban_centres":f'urban_centres = "{urban_centres}"\n',
        #                     "graphic_dict_input":f'graphic_dict_input = "{graphic_dict_input}"\n'
        #                     }
        
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
    
    weave(outfile, doctype = "pandoc", figdir=figdir)

def main():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("-i","--input-csv", required=False, help="path to input file",dest="input_csv")
    parser.add_argument("--report-args", required=True, help="report arguments text file", dest="report_args")
    parser.add_argument("-f", "--tree-fields-input",default="", help="desired fields to colour trees by in report. Default=UK country",dest="tree_fields_input")
    parser.add_argument("-l", "--label-fields-input", default="", help="fields to add into labels in report trees. Default is adm2 and date", dest='label_fields_input')
    parser.add_argument("-gd", "--graphic-dict-input", default="", help="fields to colour by rather than display text. Add colour scheme optionally", dest="graphic_dict_input")
    parser.add_argument("--date-fields-input", help="column headers containing date information as a a comma separated string.", dest="date_fields_input")
    parser.add_argument("--node-summary-option", action="store", help="field to summarise collapsed nodes by. Default=lineage", dest="node_summary_option")

    parser.add_argument("-sc", "--seq-centre",default="", help="Sequencing centre", dest="seq_centre")

    parser.add_argument("--filtered-cog-metadata", required=False, help="path to combined metadata file",dest="filtered_cog_metadata")
    parser.add_argument("--cog-metadata", required=True, help="path to full COG metadata file",dest="cog_metadata")
    
    parser.add_argument("--failed-seqs", required=False, default="", help="csv of seqs that fail qc and the reason why",dest="failed_seqs")    
    parser.add_argument("-t","--tree-dir", required=False, default="", help="path to tree directory",dest="tree_dir")
    parser.add_argument("--report-template", help="report template file",dest="report_template")

    parser.add_argument("-o","--outfile", default="civet_report.pmd", help="output name stem as a string",dest="outfile")
    parser.add_argument("--outdir", help="output directory",dest="outdir")
    parser.add_argument("--figdir", help="output directory",dest="figdir")

    parser.add_argument("--clean-locs-file", required=True, help="CSV for cleaning adm2 regions in metadata", dest="clean_locs_file")
    parser.add_argument("--uk-map", required=True, help="shape file for uk counties", dest="uk_map")
    parser.add_argument("--channels-map", required=True, help="shape file for channel islands", dest="channels_map")
    parser.add_argument("--ni-map", required=True, help="shape file for northern irish counties", dest="ni_map")
    parser.add_argument("--pc-file", required=True, help="file containing outer postcode to centroid mapping", dest="pc_file")
    parser.add_argument("--add-bars", action="store_true",dest="add_bars",default=False)

    parser.add_argument("--map-sequences", required=True, help="Bool for whether mapping of sequences by trait is required", dest="map_sequences")
    parser.add_argument("--map-cols", default="", help="either column names in input csv which contains x coords and y coords for mapping as a comma separated string OR column name containing outer postcode", dest="map_inputs")
    parser.add_argument("--input-crs", default="", help="coordinate reference system that x and y inputs are in", dest="input_crs")
    parser.add_argument("--mapping-trait", default="", help="trait to map sequences by", dest="mapping_trait")
    parser.add_argument("--urban-centres", default="", help="geojson for plotting urban centres", dest="urban_centres")

    parser.add_argument("--local-lineages", default="", action='store_true',help="List of rendered .png files for local lineage analysis", dest="local_lineages")
    parser.add_argument("--local-lin-maps", default="", action='store',help="List of rendered .png files for local lineage analysis", dest="local_lin_maps")
    parser.add_argument("--local-lin-tables", default="", action='store', help="List of .md tables for local lineage analysis", dest="local_lin_tables")

    args = parser.parse_args()

    # arg_list = get_report_arguments2(args.report_args)

    arg_list = get_arg_list(args)
    make_report(arg_list, args.report_args)


def get_arg_list(args):

    arg_list = [args.cog_metadata, args.input_csv, args.filtered_cog_metadata, 
    args.outfile, args.outdir, args.tree_dir, 
    args.figdir, args.tree_fields_input, args.label_fields_input, 
    args.node_summary_option, args.report_template, args.failed_seqs, 
    args.seq_centre, args.clean_locs_file, args.uk_map, args.channels_map,
    args.ni_map, args.pc_file, args.local_lineages, args.local_lin_maps, 
    args.local_lin_tables,args.map_sequences, args.map_inputs, args.input_crs, args.mapping_trait, 
    args.urban_centres,args.add_bars, args.graphic_dict_input, args.date_fields_input]

    return arg_list

if __name__ == "__main__":
    main()
