#!/usr/bin/env python3
import os
from civet.utils.log_colours import green,cyan


def add_col_to_metadata(col_header, dictionary, metadata): #dictionary currently is key=sequence name and value=new col value

    new_metadata = f'{metadata.strip(".csv")}_new.csv'

    with open(new_metadata, 'w') as fw:
        
        with open(metadata) as f:
            read_data = csv.DictReader(f)
            fieldnames = read_data.fieldnames
            write_fieldnames = fieldnames.append(col_header)
            write_obj = csv.DictWriter(fw, fieldnames=write_fieldnames)
            write_obj.writeheader()

            for line in read_data:
                write_dict = {}
                for field in fieldnames:
                    write_dict[field]
                if line["sample_name"] in dictionary:
                    write_dict[col_header] = dictionary[line["sample_name"]]
                else:
                    write_dict[col_header] = ""

                write_obj.writerow(write_dict)


    os.remove(metadata)
    os.rename(new_metadata, metadata)
                    


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
                                                                
                      Aine O'Toole & Verity Hill & Rachel Colquhoun       
                                    Rambaut Group              
                                Edinburgh University          
\n"""))

def preamble(v):
    header(v)
    funding()
    acknowledgements()

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
                        JT McCrone          clusterfunk     
                        Stefan Rooke        local map 
                        Andrew Rambaut      jclusterfunk    
                                                            
                    Acknowledgements:            
                                                            
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

