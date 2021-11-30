#!/usr/bin/env python3
import os
import csv
import sys
import datetime as dt

from civet.utils.log_colours import green,cyan
from civet.utils.config import *

def add_col_to_metadata(new_column_name, new_column_dict, metadata, new_metadata, match_column, config): 
    #dictionary currently is key=sequence name and value=new col value
    print(green("Adding column to master metadata:"), new_column_name)
    with open(new_metadata, 'w') as fw:
        
        with open(metadata,"r") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            header.append(new_column_name)
            config[KEY_QUERY_CSV_HEADER] = header
            
            writer = csv.DictWriter(fw, fieldnames=config[KEY_QUERY_CSV_HEADER],lineterminator='\n')
            writer.writeheader()

            for row in reader:
                new_row = row

                if new_row[match_column] in new_column_dict:
                    new_row[new_column_name] = new_column_dict[new_row[match_column]]
                else:
                    new_row[new_column_name] = ""

                writer.writerow(new_row)

def check_date_format(to_check, line_count, header):

    date_format = "%Y-%m-%d"
    
    try:
        dt.datetime.strptime(to_check, date_format).date()
    except:
        sys.stderr.write(cyan(f'Date {to_check} on line {line_count} in column {header} in incorrect format. Please use YYYY-MM-DD\n'))
        sys.exit(-1)
                
def add_arg_to_config(key,arg,config):
    if arg:
        config[key] = arg

def add_file_to_config(key,arg,config):
    if arg:
        path_to_file = os.path.abspath(config["cwd"])
        full_path = os.path.join(path_to_file,arg)
        config[key]=full_path

def add_path_to_config(key,arg,config):
    if arg:
        expanded_path = os.path.expanduser(arg)
        path_to_cwd = os.path.abspath(config["cwd"])
        full_path = os.path.join(path_to_cwd,expanded_path)
        config[key]=full_path

def header(v):
    print(green("""\n
                                    __              __    
                              ____ |__|__  __ _____/  |_ 
                             / ___\|  \  \/ // __ \   __|
                            \  \___|  |\   /\  ___/|  |  
                             \____/ __| \_/  \____/ __|  

                **** Cluster Investigation & Virus Epidemiology Tool ****
                """)+green(f"""
                                        {v}""")+green("""
                        ****************************************
                                                                
                         Aine O'Toole, Verity Hill, Ben Jackson, 
                       JT McCrone, Rachel Colquhoun, Stefan Rooke, 
                                    Andrew Rambaut        
                                 Edinburgh University          
\n"""))

def preamble(v):
    header(v)
    

def funding():
    print(green("""
                    Funding:                
                                                                
                                    ARTIC Network               
                        Wellcome Trust Collaborators Award      
                                    206298/Z/17/Z               
                                                                
                            COVID-19 Genomics UK Consortium     
                        UK Department of Health and Social Care 
                            UK Research and Innovation          
                                                                
                                    ReservoirDOCs               
                    European Research Council Consolidator Grant
                                    ERC-2016-COG                
                                                             
"""))

def acknowledgements():
    print(green("""
                    Code contributors:           
                                                            
                        Ben Jackson         gofasta       
                        Rachel Colquhoun    background data
                        JT McCrone          figtree.js     
                        Stefan Rooke        local map 
                        Andrew Rambaut      jclusterfunk    
                                                            
                """))
def full_acknowledgements():
    funding()
    print(green("""
                    Acknowledgements:            

                    Code contributors:           
                                                            
                        Ben Jackson         gofasta       
                        Rachel Colquhoun    background data
                        JT McCrone          figtree.js     
                        Stefan Rooke        local map 
                        Andrew Rambaut      jclusterfunk   
                                                            
                    We thank the following for helpful suggestions, 
                    comments, beta-testing, feature requests and
                    patience.                
                                                            
                        :nickloman:         :mattloose:     
                        :mattbashton:       :tomconnor:     
                        :rebeccadewar:      :martinmchugh:    
                        :richardmyers:      :meerachand:    
                        :samnicholls:       :radpoplawski:   
                        :davidaanensen:     :benlindsey:    
                        :jeffbarrett:       :derekfairley:   
                        :josephhughes:      :davidrobertson:  
                        :richardorton:      :mattholden:
                        :ulfschaefer:       :nataliegroves:   
                        :nikosmanesis:      :jaynaraghwani:   
"""))

def be_arty():
    logo()

def logo():
    print("""
                                       &@                                       
                           *@@,,,,,,,,,,,,,,,,,,,,,@@/                          
                      %@,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@%                     
                   @,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@                  
                @,,,,,,,,,,,,,,,,%(**(,,,,,,,,,,,,,,,,,,,,,,,,,,@               
              @,,,,,,,,,,,,,,,%(((((%**%((((((%,,,,,,,,,,,,,,,,,,,@             
            @,,,,,,,,,,,,,,,%((((((((((%**%((((%,,,,,,,,,,,,,,,,,,,,@           
          @,,,,,,,,,,,,,,,%,,,,,,%(((######((((%,,,,,,,,,,,,,,,,,,,,,,@         
         @,,,,,,,,,%(((%(****#,,,,###****##((((#,,,,,,,,,,,,,,,,,,,,,,,@        
        @,,,,,,,,,,,,,,,((((%(((((##***%#%,,,,,#,,,,,,,,,,,,,,,,,,,,,,,,@       
       @,,,,,,,,,,,,,,#(((((((((((%((*%(((((#,,*,,,,,,,,,,,,,,,,,,,,,,,,,@      
      /*,,,,,,,,,,,,%((((((((((((((((##%(*********%,,,,,,,,,,,,,,,,,,,,,,*      
      @,,,,,,,,,,,%((((((((((((((((((((####  #(%*****%,,,,,,,,,,,,,,,,,,,,@     
      @,,,,,,,,,,(((((((((((((%##%%########(**%**,,,,,,,,,,,,,,,,,,,,,,,,,@     
      @,,,,/**#((((%%***#%#(%***********##,/##%/**%,,,,,,,,,,,,,,,,,,,,,,,@     
      @,,,,,****(((((((((((%*****(, ,,,.*/***%(*/###%,,,,,,,,,,,,,,,,,,,,,@     
      @,,,,,%****((((((((((((((((#*****,,,,,,,,,,,*%%,%###%%%/,**,,,,,,,,,@     
       @,,,,%******((((((((((((((((((%***%,,,,,,,,,,###,,,,,,,/#**%,,,,,,@      
       (*,,/*********((((((((((((((((%*,,,,,,,,,,,###,,,,/*,,,*(****,,,,/.      
        (%****(%#****/(((((((((((((,,,,,,,,,,,,,,##(,,,,/(,,,%,,****,,,*#       
          @*******##((((((((((((((((((((*,,,,,,,%#(,,,,,,,,,/,,***%,(%@         
           @***************************%,,,,,,,,%#(,,,,,,,%,,%*,,,%,,@          
             @********************%/*******,,,,,,,%,,,*#/(,*,,,,,,,@            
               @********************/******%,,,,%##%,,,,,,,,,,,,,@              
                 @@*******************###((((%#####,,,,,,,,,,,@@                
                    *@************%##((***((((####%,,,,,,,,@*                   
                         @@(***/#%*****%(((((%###,,,,,@@                        
                               @@@@/*((((((((%@@@@                              

""")

