#!/usr/bin/env python3
from civet.utils import log_colours as colour


def header(v):
    print(colour.green("""\n
                                    __              __    
                              ____ |__|__  __ _____/  |_ 
                             / ___\|  \  \/ // __ \   __|
                            \  \___|  |\   /\  ___/|  |  
                             \____/ __| \_/  \____/ __|  

                **** Cluster Investigation & Virus Epidemiology Tool ****
                """)+colour.green(f"""
                                        {v}""")+colour.green("""
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
    print(colour.green("""
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
    print(colour.green("""
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
