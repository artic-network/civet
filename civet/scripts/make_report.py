
import os
from pweave import *
import argparse

def make_report(big_metadata, input_metadata, name_stem, output_directory, input_directory, desired_fields):
    fd = os.path.join(output_directory, "figures")
    pmd_file = open(name_stem + ".pmd", 'w')
    
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


    pmd_file.close()


    weave(pmd_string, doctype = "pandoc", figdir=fd)

def main():
    parser = argparse.ArgumentParser(description="Report generator script")
    parser.add_argument("--i", required=False, default="", help="path to tree inputs as a string")
    parser.add_argument("--bm", required=True, help="path to full COG metadata file")
    parser.add_argument("--im", required=False, default="", help="path to their metadata file")
    parser.add_argument("--s", default="civet_report", help="output name stem as a string")
    parser.add_argument("--od", default="./", help="output directory, default is working directory")
    parser.add_argument("--df", default=[], help="desired fields that you want to be in the analysis, otherwise will just have date and UK country")


    args = parser.parse_args()

    make_report(args.bm, args.im, args.s, args.od, args.i, args.df)


if __name__ == "__main__":
    main()