
def timeline_checking(metadata, timeline_dates, config):  

    misc.add_arg_to_config("timeline_dates",timeline_dates,config) #default is None

    if config["timeline_dates"]:

        date_header_list = config["timeline_dates"].split(",")

        with open(metadata) as f:
            data = csv.DictReader(f)
            for header in date_header_list:
                if header not in data.fieldnames:
                    sys.stderr.write(cyan(f"Error: {header} (specified for use in timeline plot) not found in metadata file.\n"))
                    sys.exit(-1)

            line_count = 0
            for line in data:
                line_count += 1
                for header in date_header_list:
                    if line[header] != "":
                        misc.check_date_format(line[header], line_count, header)


    return config



def make_timeline_json(config):

    date_cols = config["timeline_dates"].split(",")

    overall = defaultdict(dict)
    overall['catchments'] = defaultdict(dict)

    with open(config["query_metadata"]) as f:
        data = csv.DictReader(f)
        for l in data:
            if l['query_boolean'] == "TRUE":
                if l['catchment'] in overall['catchments']:
                    dict_list = overall['catchments'][l['catchment']]
                else:
                    dict_list = []
                
                new_dict = {}
                new_dict["sequence_name"] = l['display_name'] 
                for col in date_cols:
                    new_dict[col] = l[col]
                
                dict_list.append(new_dict)
                overall['catchments'][l['catchment']] = dict_list

    with open(os.path.join(config["tempdir"],'timeline_data.json'), 'w') as outfile:
        json.dump(overall, outfile)