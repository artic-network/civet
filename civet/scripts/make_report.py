
import os
from pweave import *
import argparse
import shutil


thisdir = os.path.abspath(os.path.dirname(__file__))

def make_report(cog_metadata, input_csv, filtered_cog_metadata, outfile, outdir, treedir, figdir, fields, report_template, font_file):

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
                    elif "desired_fields" in l:
                        new_l = f'desired_fields = "{fields}"\n'
                    elif "font_file" in l:
                        new_l = f'font_file = "{font_file}"\n'
                    elif "figdir" in l:
                        new_l = f'figdir = "{figdir}"\n'
                    elif "tree_dir" in l:
                        new_l = f'tree_dir = "{treedir}"\n'
                    elif "summary_dir" in l:
                        new_l = f'summary_dir = "{summary_dir}"\n'
                    
                else:
                    new_l = l

                pmd_file.write(new_l)

    weave(outfile, doctype = "pandoc", figdir=figdir)

def main():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("-i","--input-csv", required=False, help="path to input file",dest="input_csv")
    parser.add_argument("-f", "--fields",default="", help="desired fields for report. Default=date and UK country",dest="fields")

    parser.add_argument("--filtered-cog-metadata", required=False, help="path to combined metadata file",dest="filtered_cog_metadata")
    parser.add_argument("--cog-metadata", required=True, help="path to full COG metadata file",dest="cog_metadata")
    
    parser.add_argument("-t","--treedir", required=False, default="", help="path to tree directory",dest="treedir")
    parser.add_argument("--report-template", help="report template file",dest="report_template")
    parser.add_argument("--font-file", help="custom font",dest="font_file")

    parser.add_argument("-o","--outfile", default="civet_report.pmd", help="output name stem as a string",dest="outfile")
    parser.add_argument("--outdir", help="output directory",dest="outdir")
    parser.add_argument("--figdir", help="output directory",dest="figdir")

    args = parser.parse_args()

    make_report(args.cog_metadata, args.input_csv, args.filtered_cog_metadata, args.outfile, args.outdir, args.treedir, args.figdir,args.fields, args.report_template, args.font_file)


if __name__ == "__main__":
    main()