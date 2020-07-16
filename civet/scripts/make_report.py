#!/usr/bin/env python3
import os
from pweave import *
import argparse
import shutil


thisdir = os.path.abspath(os.path.dirname(__file__))

def make_report(cog_metadata, input_csv, filtered_cog_metadata, outfile, outdir, treedir, figdir, fields, report_template, failed_seqs, no_seq, seq_centre, clean_locs, uk_map, channels_map, ni_map):

    name_stem = ".".join(outfile.split(".")[:-1])
    with open(outfile, 'w') as pmd_file:
    
        md_template = report_template
        summary_dir = os.path.join(outdir, "summary_files")
        with open(md_template) as f:
            for l in f:
                if "##CHANGE" in l:
                    if "output_directory" in l:
                        new_l = f'output_directory = "{outdir}"\n'
                    elif "name_stem_input" in l:
                        new_l = f'name_stem_input = "{name_stem}"\n'
                    elif "full_metadata_file" in l:
                        new_l = f'full_metadata_file = "{cog_metadata}"\n'
                    elif "filtered_cog_metadata" in l:
                        new_l = f'filtered_cog_metadata = "{filtered_cog_metadata}"\n'
                    elif "input_csv" in l:
                        new_l = f'input_csv = "{input_csv}"\n'
                    elif "input_directory" in l:
                        new_l = f'input_directory = "{treedir}"\n'
                    elif "desired_fields_input" in l:
                        new_l = f'desired_fields_input = "{fields}"\n'
                    elif "figdir" in l:
                        new_l = f'figdir = "{figdir}"\n'
                    elif "tree_dir" in l:
                        new_l = f'tree_dir = "{treedir}"\n'
                    elif "summary_dir" in l:
                        new_l = f'summary_dir = "{summary_dir}"\n'
                    elif "QC_fail_file" in l:
                        new_l = f'QC_fail_file = "{failed_seqs}"\n'
                    elif "missing_seq_file" in l:
                        new_l = f'missing_seq_file = "{no_seq}"\n'
                    elif "sequencing_centre" in l:
                        new_l = f'sequencing_centre = "{seq_centre}"\n'
                    elif "clean_locs_file" in l:
                        new_l = f'clean_locs_file = "{clean_locs}"\n'
                    elif "uk_map" in l:
                        new_l = f'uk_map = "{uk_map}"\n'
                    elif "channels_map" in l:
                        new_l = f'channels_map = "{channels_map}"\n'
                    elif "ni_map" in l:
                        new_l = f'ni_map = "{ni_map}"\n'
                    elif "local_lineages" in l:
                        new_l = f'local_lineages = "{local_lineages}"\n'
                    elif "local_lin_maps" in l:
                        new_l = f'local_lin_maps = "{local_lin_maps}"\n'
                    elif "local_lin_tables" in l:
                        new_l = f'local_lin_tables = "{local_lin_tables}"\n'
                else:
                    new_l = l

                pmd_file.write(new_l)

    weave(outfile, doctype = "pandoc", figdir=figdir)

def main():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("-i","--input-csv", required=False, help="path to input file",dest="input_csv")
    parser.add_argument("-f", "--fields",default="", help="desired fields for report. Default=date and UK country",dest="fields")
    parser.add_argument("-sc", "--sequencing-centre",default="", help="Sequencing centre", dest="sc")

    parser.add_argument("--filtered-cog-metadata", required=False, help="path to combined metadata file",dest="filtered_cog_metadata")
    parser.add_argument("--cog-metadata", required=True, help="path to full COG metadata file",dest="cog_metadata")
    
    parser.add_argument("--failed-seqs", required=False, default="", help="csv of seqs that fail qc and the reason why",dest="failed_seqs")
    parser.add_argument("--no-seq-provided", required=False, default="", help="file of seqs that weren't in cog and didn't have a sequence provided",dest="no_seq")
    
    parser.add_argument("-t","--treedir", required=False, default="", help="path to tree directory",dest="treedir")
    parser.add_argument("--report-template", help="report template file",dest="report_template")

    parser.add_argument("-o","--outfile", default="civet_report.pmd", help="output name stem as a string",dest="outfile")
    parser.add_argument("--outdir", help="output directory",dest="outdir")
    parser.add_argument("--figdir", help="output directory",dest="figdir")

    parser.add_argument("--clean-locs", required=True, help="CSV for cleaning adm2 regions in metadata", dest="clean_locs")
    parser.add_argument("--uk-map", required=True, help="shape file for uk counties", dest="uk_map")
    parser.add_argument("--channels-map", required=True, help="shape file for channel islands", dest="channels_map")
    parser.add_argument("--ni-map", required=True, help="shape file for northern irish counties", dest="ni_map")

    parser.add_argument("--local_lineages", default='False', action='store_true',help="List of rendered .png files for local lineage analysis. ';' delimited.", dest="local_lineages")
    parser.add_argument("--local_lin_maps", default='None', action='store',help="List of rendered .png files for local lineage analysis. ';' delimited.", dest="local_lin_maps")
    parser.add_argument("--local_lin_tables", default='None', action='store', help="List of .md tables for local lineage analysis. ';' delimited. ", dest="local_lin_tables")


    args = parser.parse_args()

    make_report(args.cog_metadata, args.input_csv, args.filtered_cog_metadata, args.outfile, args.outdir, args.treedir, args.figdir,args.fields, args.report_template, args.failed_seqs,args.no_seq, args.sc, args.clean_locs, args.uk_map, args.channels_map, args.ni_map)


if __name__ == "__main__":
    main()