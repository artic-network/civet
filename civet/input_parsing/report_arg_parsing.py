import sys
import csv
from civet.utils.log_colours import green,cyan
from civet.utils import misc
import datetime as dt

def qc_report_content(config):
    reports = config["report_content"]
    to_generate = []
    to_run = []
    for report_options in reports:
        report_options = report_options.split(",")
        try:
            report_options = [int(i) for i in report_options]
        except:
            sys.stderr.write(cyan(f'Error: -rc/ --report-content should be a comma separated numerical string.\n'))
            sys.exit(-1)
        
        for i in report_options:
            if i in range(1,7):
                to_run.append(i)
            else:
                sys.stderr.write(cyan(f'Error: {i} is an invalid -rc/ --report-content option. Options 1..6 inclusive.\n'))
                sys.exit(-1)

        report_options = sorted(list(set(report_options)))
        to_generate.append(report_options)
    to_run = sorted(list(set(to_run)))
    print(to_run, to_generate)
    config["report_content"] = to_run
    config["reports"] = to_generate

#then at some point we need to update the treefile with these display names using jclusterfunk

def report_group_parsing(report_content,anonymise,config):
    """
    parses the report group arguments 
    --report-content (Default 1,2,3)
    --anonymise (Default False)
    """

    # if command line arg, overwrite config value
    misc.add_arg_to_config("report_content",report_content,config)
    misc.add_arg_to_config("anonymise",anonymise,config)

    qc_report_content(config)
