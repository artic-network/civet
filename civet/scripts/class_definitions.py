


def convert_date(date_string):
    bits = date_string.split("-")
    date_dt = dt.date(int(bits[0]),int(bits[1]), int(bits[2]))
    
    return date_dt

class taxon():

    def __init__(self, name, global_lin, uk_lin, phylotype, label_fields, tree_fields):

        self.name = name

        self.sample_date = "NA"

        self.date_dict = {}
        
        if global_lin == "":
            self.global_lin = "NA"
        else:
            self.global_lin = global_lin
        
        if uk_lin == "":
            self.uk_lin = "NA"
        else:
            self.uk_lin = uk_lin
       
        if phylotype == "":
            self.phylotype = "NA"
        else:
            self.phylotype = phylotype
       
        self.in_cog = False
        
        self.attribute_dict = {}
        self.attribute_dict["adm1"] = "NA"
        
        for i in label_fields:
            self.attribute_dict[i] = "NA"
        for i in tree_fields:
            self.attribute_dict[i] = "NA"
        
        
        self.tree = "NA"

        self.closest_distance = "NA"
        self.snps = "NA"

class lineage():
    
    def __init__(self, name, taxa):
        
        self.name = name
        self.taxa = taxa
        self.dates = []
        self.global_lins = set()
        
        for tax in taxa:
            if tax.sample_date != "NA":
                tax.date_dt = convert_date(tax.sample_date)
                self.dates.append(tax.date_dt)
            self.global_lins.add(tax.global_lin)
                
        if self.dates == []:
            self.first_date = "NA"
        else:
            self.first_date = min(self.dates)