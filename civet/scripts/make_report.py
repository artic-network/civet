
import os
from pweave import *
import argparse

def make_report(big_metadata, input_metadata, name_stem, output_directory, input_directory, desired_fields):
    fd = os.path.join(output_directory, "figures")

    outfile = os.path.join(output_directory, f"{name_stem}.pmd")
    with open(outfile, 'w') as pmd_file:
    
        pmd_string = name_stem + ".pmd"

        md_template = 'civet_template.pmd'
        
        with open(md_template) as f:
            for l in f:
                if "##CHANGE" in l:
                    if "output_directory" in l:
                        new_l = 'output_directory = "' + str(output_directory) + '"\n'
                    elif "name_stem" in l:
                        new_l = 'name_stem = "' + str(name_stem) + '"\n'
                    elif "big_metadata" in l:
                        new_l = 'big_metadata = "' + str(big_metadata) + '"\n'
                    elif "their_metadata" in l:
                        new_l = 'their_metadata = "' + str(input_metadata) + '"\n'
                    elif "input_directory" in l:
                        new_l = 'input_directory = "' + str(input_directory) + '"\n'
                    elif "desired_fields" in l:
                        new_l = 'desired_fields = "' + str(desired_fields) + '"\n'

                else:
                    new_l = l

                pmd_file.write(new_l)

    weave(pmd_string, doctype = "pandoc", figdir=fd)

def main():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("-i","--input-path", required=False, default="", help="path to tree and annotations inputs as a string",dest="input_path")
    parser.add_argument("-bm","--cog-metadata", required=True, help="path to full COG metadata file",dest="cog_metadata")
    parser.add_argument("-im","--input-metadata", required=False, help="path to their metadata file",dest="input_metadata")
    parser.add_argument("-s","--stem", default="civet_report", help="output name stem as a string",dest="stem")
    parser.add_argument("-od","--outdir", help="output directory",dest="outdir")
    parser.add_argument("-df", "--desired-fields",default=[], help="desired fields for report. Default=date and UK country",dest="desired_fields")


    args = parser.parse_args()

    make_report(args.cog_metadata, args.input_metadata, args.stem, args.outdir, args.input_path, args.desired_fields)


if __name__ == "__main__":
    main()