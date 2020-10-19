from civet import __version__

import setuptools
import argparse
import os.path
import snakemake
import sys
import tempfile
import csv
import os
import yaml
from datetime import datetime
from datetime import date
from Bio import SeqIO
import pkg_resources

from reportfunk.funks import io_functions as qcfunk
from reportfunk.funks import report_functions as rfunk
from reportfunk.funks import custom_logger as custom_logger
from reportfunk.funks import log_handler_handle as lh
from civet.scripts import civetfunks as cfunk


import pytest

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))
print("Running tests")

def test_outdir_default():
    config = cfunk.get_defaults()
    qcfunk.get_outdir(None,None,cwd,config)
    outdir = config["outdir"].split("/")[-1]
    assert outdir.startswith("civet")

def test_outdir_arg():
    config = cfunk.get_defaults()
    qcfunk.get_outdir("input_arg",None,cwd,config)
    outdir = config["outdir"].split("/")[-1]
    assert outdir.startswith("input_arg")

def test_output_prefix():
    config = cfunk.get_defaults()
    qcfunk.get_outdir(None,"prefix",cwd,config)
    outdir = config["outdir"].split("/")[-1]
    assert outdir.startswith("prefix")

def test_outdir_and_prefix():
    config = cfunk.get_defaults()
    qcfunk.get_outdir("outdir_arg","prefix",cwd,config)
    outdir = config["outdir"].split("/")[-1]
    assert outdir.startswith("outdir_arg")

def test_configure_update():
    config = cfunk.get_defaults()
    config["from_metadata"]="sample_date=2020-10-04"
    cfunk.configure_update(True,config)
    assert config["colour_by"]=="new:Paired"
    assert config["tree_fields"]=="new"
    assert config["table_fields"]=="central_sample_id,sequence_name,sample_date,uk_lineage,phylotype,tree,new"
    assert config["tempdir"]==""
    assert config["outdir"]==""

def test_configure_cluster():
    config = cfunk.get_defaults()
    cfunk.configure_cluster(config)
    assert config["down_distance"]==100

def test_get_sequencing_centre_header():
    config = cfunk.get_defaults()
    qcfunk.get_outdir(None,None,cwd,config)
    config["sequencing_centre"] = "EDIN"
    cfunk.get_sequencing_centre_header(config)
    assert config["sequencing_centre_file"] == "./figures/EDIN.png"

def test_min_len():
    config = cfunk.get_defaults()
    print(cwd)
    config["query"]=os.path.join(thisdir,"test.csv")
    config["fasta"]=os.path.join(thisdir,"test.fasta")
    config["outdir"] = thisdir
    config["background_metadata"] = os.path.join(thisdir,"test_background.csv")
    num_seqs = qcfunk.input_file_qc(40000,None,config)
    assert num_seqs == 0

def test_max_ambig():
    config = cfunk.get_defaults()
    print(cwd)
    config["query"]=os.path.join(thisdir,"test.csv")
    config["fasta"]=os.path.join(thisdir,"test.fasta")
    config["outdir"] = thisdir
    config["background_metadata"] = os.path.join(thisdir,"test_background.csv")
    num_seqs = qcfunk.input_file_qc(None,0,config)
    assert num_seqs == 2

def test_input_file_qc():
    config = cfunk.get_defaults()
    print(cwd)
    config["query"]=os.path.join(thisdir,"test.csv")
    config["fasta"]=os.path.join(thisdir,"test.fasta")
    config["outdir"] = thisdir
    config["background_metadata"] = os.path.join(thisdir,"test_background.csv")
    num_seqs = qcfunk.input_file_qc(None,None,config)
    assert num_seqs == 2
    test_min_len()
    test_max_ambig()

def test_background_metadata_header():
    config = cfunk.get_defaults()
    config["background_metadata"] = os.path.join(thisdir,"test_background.csv")
    test_header = "central_sample_id,biosample_source_id,sequence_name,secondary_identifier,sample_date,epi_week,country,adm1,is_surveillance,is_community,is_hcw,is_travel_history,travel_history,lineage,lineage_support,uk_lineage,acc_lineage,del_lineage,phylotype,adm2".split(",")
    qcfunk.check_metadata_for_search_columns(config)
    assert test_header == config["background_metadata_header"]
    
def test_generate_query_from_metadata():
    config = cfunk.get_defaults()
    config["background_metadata"] = os.path.join(thisdir,"test_background.csv")
    from_metadata = "adm2=Edinburgh sample_date=2020-10-01"

    query = qcfunk.generate_query_from_metadata("test_from_metadata.csv",from_metadata, config["background_metadata"], config)
    c = 0
    with open(query, "r") as f:
        for l in f:
            c +=1
    assert c == 5

    from_metadata = "adm2=Edinburgh sample_date=2020-10-01:2020-10-10"

    query = qcfunk.generate_query_from_metadata("test_from_metadata.csv",from_metadata, config["background_metadata"], config)
    c = 0
    with open(query, "r") as f:
        for l in f:
            c +=1
    assert c == 12

def test_check_query_for_input_column():
    config = cfunk.get_defaults()
    config["query"]=os.path.join(thisdir,"test.csv")
    qcfunk.check_query_for_input_column(config)
    test_query_header = "name,display,HCW_status,care_home,sample_date,date_2,date_3,adm2,test,test2,outer_postcode,x_col,y_col,adm1_blah".split(",")
    assert config["query_metadata_header"] == test_query_header

def test_map_sequences_config():
    config = cfunk.get_defaults()
    config["background_metadata"] = os.path.join(thisdir,"test_background.csv")
    qcfunk.check_metadata_for_search_columns(config)
    config["query_metadata_header"] =  "name,display,HCW_status,care_home,sample_date,date_2,date_3,adm2,test,test2,outer_postcode,x_col,y_col,adm1_blah".split(",")
    config["map_info"] = "outer_postcode"
    config["map_sequences"] = True
    cfunk.map_sequences_config(config)
    assert config["input_crs"] == "EPSG:4326"

def test_check_date_format():
    date_string1 = "2020-10-01"
    check_date = qcfunk.check_date_format(date_string1)
    print(check_date)
    assert str(check_date) == "2020-10-01"

    date_string2 = "10-01-2020"
    try:
        check_date = qcfunk.check_date_format(date_string2)
    except:
        check_date = -1
    assert check_date == -1

def main():
    test_outdir_default()

    test_outdir_arg()

    test_output_prefix()

    test_outdir_and_prefix()

    test_configure_update()

    test_configure_cluster()

    test_get_sequencing_centre_header()

    test_input_file_qc()

    test_background_metadata_header()

    test_generate_query_from_metadata()

    test_check_query_for_input_column()

    test_map_sequences_config()

    test_check_date_format()

if __name__ == '__main__':
    main()