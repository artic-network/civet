# QUERY KEYS
KEY_IDS="ids"
KEY_INPUT_METADATA="input_metadata"

# INPUT OPTION KEYS
KEY_FROM_METADATA="from_metadata"
KEY_MAX_QUERIES="max_queries"

# INPUT COLUMN KEYS
KEY_INPUT_ID_COLUMN="input_id_column"
KEY_INPUT_DATE_COLUMN="input_date_column"
KEY_INPUT_DISPLAY_COLUMN="input_display_column"

# INPUT SEQ OPTION KEYS
KEY_INPUT_SEQUENCES="input_sequences"
KEY_NUM_SEQS="num_seqs"
KEY_MAX_AMBIGUITY="max_ambiguity"
KEY_MIN_LENGTH="min_length"
KEY_TRIM_START="trim_start"
KEY_TRIM_END="trim_end"

# BACKRGOUND VARIABLES KEYS
CIVET_DATADIR="datadir"
KEY_BACKGROUND_METADATA="background_metadata"
KEY_BACKGROUND_SEQUENCES="background_sequences"
KEY_BACKGROUND_SNPS="background_snps"
KEY_BACKGROUND_TREE="background_tree"

# BACKGROUND COLUMN KEYS
KEY_BACKGROUND_ID_COLUMN="background_id_column"


# DIRECTORIES
KEY_TEMPDIR="tempdir"
KEY_INPUT_PATH="input_path"
KEY_CWD="cwd"
KEY_OUTDIR="outdir"

# QUERY MAP KEYS
KEY_QUERY_MAP_FILE="query_map_file"
KEY_LATITUDE_COLUMN="latitude_column"
KEY_LONGITUDE_COLUMN="longitude_column"
KEY_BACKGROUND_MAP_DATE_RANGE="background_map_date_range"

# MISC KEYS
KEY_ANONYMISE="anonymise"
KEY_VERBOSE="verbose"
KEY_THREADS="threads"
KEY_LOG_API="log_api"
KEY_LOG_STRING="log_string"
KEY_QUIET="quiet"
KEY_CIVET_MODE="civet_mode"
KEY_DATE="date"
KEY_AUTHORS="authors"

RESOURCE_KEY_FILENAME="filename"
RESOURCE_KEY_DIRECTORY="directory"
RESOURCE_KEY="KEY"

# HEADER KEYS
KEY_QUERY_CSV_HEADER="query_csv_header"

PROTECTED_COL_NAMES = ["hash","catchment","query_boolean","qc_status","source","seq_N_content","seq_length","in_tree"]

# DEPENDENCIES AND RESOURCES TO CHECK
dependency_list = ["gofasta","minimap2","snakemake","iqtree","jclusterfunk"]
module_list = ["mako","Bio"]

resources = [
        {RESOURCE_KEY:"reference_sequence",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"reference.fasta"},
        {RESOURCE_KEY:"outgroup_fasta",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"outgroup.fasta"},
        {RESOURCE_KEY:"report_template",
        RESOURCE_KEY_DIRECTORY:"data/report_modules",
        RESOURCE_KEY_FILENAME:"report_template.mako"},
        {RESOURCE_KEY:"mako_query_table",
        RESOURCE_KEY_DIRECTORY:"data/report_modules",
        RESOURCE_KEY_FILENAME:"query_table.txt"},
        {RESOURCE_KEY:"mako_timeline",
        RESOURCE_KEY_DIRECTORY:"data/report_modules",
        RESOURCE_KEY_FILENAME:"timeline.txt"},
        {RESOURCE_KEY:"html_colours",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"html_colours.csv"},
        {RESOURCE_KEY:"global_acceptable_values",
        RESOURCE_KEY_DIRECTORY:"data/map_data",
        RESOURCE_KEY_FILENAME:"global_acceptable_values.tsv"},
        {RESOURCE_KEY:"uk_acceptable_values",
        RESOURCE_KEY_DIRECTORY:"data/map_data",
        RESOURCE_KEY_FILENAME:"UK_acceptable_values.tsv"},
        {RESOURCE_KEY:"adm0_centroids",
        RESOURCE_KEY_DIRECTORY:"data/map_data",
        RESOURCE_KEY_FILENAME:"adm0_centroids.csv"},
        {RESOURCE_KEY:"adm1_centroids",
        RESOURCE_KEY_DIRECTORY:"data/map_data",
        RESOURCE_KEY_FILENAME:"adm1_centroids.csv"},
        {RESOURCE_KEY:"uk_centroids",
        RESOURCE_KEY_DIRECTORY:"data/map_data",
        RESOURCE_KEY_FILENAME:"uk_centroids.csv"}
    ]