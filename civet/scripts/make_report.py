
import os
from pweave import *
import argparse


thisdir = os.path.abspath(os.path.dirname(__file__))

def make_report(cog_metadata, input_csv, outfile, outdir, tree_dir, fields, report_template, font_file):

    fd = os.path.join(outdir, "figures")
    name_stem = ".".join(outfile.split(".")[:-1])

    with open(outfile, 'w') as pmd_file:
    
        md_template = report_template
        
        with open(md_template) as f:
            for l in f:
                if "##CHANGE" in l:
                    if "output_directory" in l:
                        new_l = f'output_directory = "{outdir}"\n'
                    elif "name_stem" in l:
                        new_l = f'name_stem = "{name_stem}"\n'
                    elif "big_metadata" in l:
                        new_l = f'big_metadata = "{cog_metadata}"\n'
                    elif "their_metadata" in l:
                        new_l = f'their_metadata = "{input_csv}"\n'
                    elif "input_directory" in l:
                        new_l = f'input_directory = "{tree_dir}"\n'
                    elif "desired_fields" in l:
                        new_l = f'desired_fields = "{fields}"\n'
                    elif "font_file" in l:
                        new_l = f'font_file = "{font_file}"\n'
                else:
                    new_l = l

                pmd_file.write(new_l)

    weave(outfile, doctype = "pandoc", figdir=fd)

def main():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("-i","--input-csv", required=False, help="path to input file",dest="input_csv")
    parser.add_argument("-f", "--fields",default="", help="desired fields for report. Default=date and UK country",dest="fields")

    parser.add_argument("--filtered-cog-metadata", required=False, help="path to combined metadata file",dest="filtered_cog_metadata")
    parser.add_argument("--cog-metadata", required=True, help="path to full COG metadata file",dest="cog_metadata")
    
    parser.add_argument("-t","--tree-dir", required=False, default="", help="path to tree directory",dest="tree_dir")
    parser.add_argument("--report-template", help="report template file",dest="report_template")
    parser.add_argument("--font-file", help="custom font",dest="font_file")
    
    parser.add_argument("-o","--outfile", default="civet_report.pmd", help="output name stem as a string",dest="outfile")
    parser.add_argument("--outdir", help="output directory",dest="outdir")

    args = parser.parse_args()

    make_report(args.cog_metadata, args.input_csv, args.outfile, args.outdir, args.tree_dir, args.fields, args.report_template, args.font_file)


if __name__ == "__main__":
    main()