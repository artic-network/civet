#!/usr/bin/env python3
from civet.utils.log_colours import green,cyan
from civet.utils import misc

import sys
import os

from datetime import datetime 
from datetime import date
import tempfile
import sys




def output_group_parsing(outdir,output_prefix,output_data,tempdir,no_temp,config):
    
    misc.add_path_to_config("outdir",outdir,config)
    misc.add_arg_to_config("output_prefix",output_prefix,config)
    misc.add_arg_to_config("output_data",output_data,config)
    misc.add_path_to_config("tempdir",tempdir,config)
    misc.add_arg_to_config("no_temp",no_temp,config)

