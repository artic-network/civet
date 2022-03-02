#!/usr/bin/env python3
import subprocess
import os
import sys
from civet.utils.log_colours import green,cyan
import importlib

def which(dependency):
    try:
        subprocess.check_output(["which", dependency])
        return True
    except subprocess.CalledProcessError:
        return False

def check_module(module, missing):
    try:
        importlib.import_module(module)
    except ImportError:
        missing.append(module)

def check_this_dependency(dependency,missing):
    check = which(dependency)

    if not check:
        missing.append(dependency)

def check_dependencies(dependency_list, module_list):

    missing = []

    for dependency in dependency_list:
        check_this_dependency(dependency, missing)

    for module in module_list:
        check_module(module, missing)

    if missing:
        if len(missing)==1:
            sys.stderr.write(cyan(f'Error: Missing dependency `{missing[0]}`.')+'\nPlease update your civet environment.\n')
            sys.exit(-1)
        else:
            dependencies = ""
            for i in missing:
                dependencies+=f"\t- {i}\n"

            sys.stderr.write(cyan(f'Error: Missing dependencies.')+f'\n{dependencies}Please update your civet environment.\n')
            sys.exit(-1)

def check_scorpio_mutations(comma_separated_mutations):
    mutations = " ".join(comma_separated_mutations.split(","))
    test_fasta = os.path.join("../tests/action_test_data/civet_test.fa")
    command = "scorpio haplotype \
                  -i %s \
                  --mutations %s \
                  --append-genotypes \
                  -n mutations \
                  --dry-run" % (test_fasta, mutations)
    completed_process = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                                           universal_newlines=True)
    if completed_process.returncode != 0:
        sys.stderr.write(cyan(
            f'Error: Scorpio does not like your mutations\n%s' % completed_process.stderr))
        sys.exit(-1)
    else:
        print(green("Scorpio mutations are in approved format."))

# check_dependencies()
