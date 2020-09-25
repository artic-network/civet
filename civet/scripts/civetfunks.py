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
                    "date_range_start":False,
                    "date_range_end":False,
                    "date_window":7,
                    "global_search":False,
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
