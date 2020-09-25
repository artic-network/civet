def get_defaults():
    default_dict = {"maxambig":0.5,
                    "minlen":10000,
                    "datadir":"civet-cat",
                    "input_column":"name",
                    "data_column":"central_sample_id",
                    "display_name":None,
                    "private":False,
                    "distance":2,
                    "up_distance":2,
                    "down_distance":2,
                    "threshold":1,
                    "add_bars":False,
                    "tree_fields":"adm1",
                    "date_range_start":False,
                    "date_range_end":False,
                    "date_window":7,
                    "global_search":False,
                    "label_fields":"NONE",
                    "date_fields":"NONE",
                    "graphic_dict":"adm1",
                    "date_restriction":False,
                    "local_lineages":False,
                    "map_sequences":False,
                    "delay_collapse":False,
                    "threads":1,
                    "force":True,
                    "trim_start":265,   # where to pad to using datafunk
                    "trim_end":29674
                    }
    return default_dict

def define_seq_db(config, default_dict):
    if config["global_search"] == True:
        config["seq_db"] = config["cog_global_seqs"]
    else:
        config["seq_db"] = config["cog_seqs"]

