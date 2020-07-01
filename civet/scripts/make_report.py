#!/usr/bin/env python3
import os
from pweave import *
import argparse
import shutil


thisdir = os.path.abspath(os.path.dirname(__file__))

def make_report(cog_metadata, input_csv, filtered_cog_metadata, outfile, outdir, treedir, figdir, fields, report_template, failed_seqs, no_seq, seq_centre):

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
                    ### Added
                    elif "local_lineage" in l:
                        new_l = f'local_lineage = "{local_lin}"\n'
                    elif "date_restriction" in l:
                        new_l = f'date_restriction = "{date_restrict}"\n'
                    elif "date_start" in l:
                        new_l = f'date_start = "{date_start}"\n'
                    elif "date_end" in l:
                        new_l = f'date_end = "{date_end}"\n'
                    elif "date_window_size" in l:
                        new_l = f'date_window_size = "{date_window}"\n'
                    
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
    ##Added
    parser.add_argument("--local-lin", default="False", help="Whether to produce localised lineages",dest="local_lin")
    parser.add_argument("--date-restrict", default="False", help="Whether to restrict dates",dest="date_restrict")
    parser.add_argument("--date-start", help="Start for date filter",dest="date_start")
    parser.add_argument("--date-end", help="End for date filter",dest="date_end")
    parser.add_argument("--date-window", help="Window for date filtering",dest="date_window")

    args = parser.parse_args()

    make_report(args.cog_metadata, args.input_csv, args.filtered_cog_metadata, args.outfile, args.outdir, args.treedir, args.figdir,args.fields, args.report_template, args.failed_seqs,args.no_seq, args.sc)


if __name__ == "__main__":
    main()