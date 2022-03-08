# QUERY KEYS
KEY_IDS="ids"
KEY_INPUT_METADATA="input_metadata"
KEY_QUERY_METADATA="query_metadata"

# INPUT OPTION KEYS
KEY_FROM_METADATA="from_metadata"
KEY_MAX_QUERIES="max_queries"

# INPUT COLUMN KEYS
KEY_INPUT_ID_COLUMN="input_id_column"
KEY_INPUT_DATE_COLUMN="input_date_column"
KEY_INPUT_DISPLAY_COLUMN="input_display_column"
KEY_GLOBAL_SNIPIT="global_snipit"
KEY_FOCAL_ALIGNMENT="focal_alignment"
KEY_REFERENCE_NAME="reference_name"

# INPUT SEQ OPTION KEYS
KEY_INPUT_SEQUENCES="input_sequences"
KEY_NUM_SEQS="num_seqs"
KEY_MAX_AMBIGUITY="max_ambiguity"
KEY_MIN_LENGTH="min_length"
KEY_TRIM_START="trim_start"
KEY_TRIM_END="trim_end"

# SEQ RELATED KEYS
KEY_QUERY_FASTA="query_fasta"
KEY_MATCHED_FASTA="matched_fasta"

# COLUMN NAME KEYS
KEY_HASH="hash"
KEY_QUERY_BOOLEAN="query_boolean"
KEY_QUERY="query"
KEY_IN_TREE="in_tree"
KEY_CATCHMENT="catchment"

# BACKRGOUND VARIABLES KEYS
CIVET_DATADIR="datadir"
KEY_BACKGROUND_METADATA="background_metadata"
KEY_BACKGROUND_SEQUENCES="background_sequences"
KEY_BACKGROUND_SNPS="background_snps"
KEY_BACKGROUND_TREE="background_tree"

# BACKGROUND COLUMN KEYS
KEY_BACKGROUND_ID_COLUMN="background_id_column"
KEY_SEQUENCE_ID_COLUMN="sequence_id_column"
KEY_BACKGROUND_DATE_COLUMN="background_date_column"
KEY_BACKGROUND_LOCATION_COLUMN="background_location_column"

# OUTPUT AND DIRECTORIES
KEY_TEMPDIR="tempdir"
KEY_NO_TEMP="no_temp"
KEY_INPUT_PATH="input_path"
KEY_CWD="cwd"
KEY_OUTDIR="outdir"
KEY_DATA_OUTDIR="data_outdir"
KEY_OUTPUT_PREFIX="output_prefix"
KEY_OUTPUT_DATA="output_data"
KEY_DATESTAMP="datestamp"
KEY_OVERWRITE="overwrite"

# BACKGROUND CURATION KEYS
KEY_UNALIGNED_SEQUENCES="unaligned_sequences"
KEY_BACKGROUND_DATA_OUTDIR="background_data_outdir"
KEY_PRIMARY_FIELD_DELIMTER="primary_field_delimiter"
KEY_PRIMARY_METADATA_FIELDS="primary_metadata_fields"
KEY_SECONDARY_FIELDS="secondary_fields"
KEY_SECONDARY_FIELD_DELIMTER="secondary_field_delimiter"
KEY_SECONDARY_FIELD_LOCATION="secondary_field_location"
KEY_SECONDARY_METADATA_FIELDS="secondary_metadata_fields"


# CATCHMENT KEYS
KEY_CATCHMENT_COUNT="catchment_count"
KEY_FIGURE_CATCHMENTS="figure_catchments"
KEY_CATCHMENT_BACKGROUND_SIZE="catchment_background_size"
KEY_SNP_DISTANCE="snp_distance"
KEY_SNP_DISTANCE_UP="snp_distance_up"
KEY_SNP_DISTANCE_DOWN="snp_distance_down"
KEY_SNP_DISTANCE_SIDE="snp_distance_side"
KEY_COLLAPSE_THRESHOLD="collapse_threshold"

KEY_DOWNSAMPLE="downsample"
KEY_MODE="mode"
KEY_FACTOR="factor"
KEY_DOWNSAMPLE_FIELD="downsample_field"
KEY_DOWNSAMPLE_COLUMN="downsample_column"

# REPORT KEYS
KEY_COLOUR_THEME="colour_theme"
KEY_COLOUR_MAP="colour_map"
KEY_REPORT_CONTENT="report_content"
KEY_OUTPUT_REPORTS="output_reports"
KEY_REPORTS="reports"
KEY_REPORT_PRESET="report_preset"
KEY_REPORT_TITLE="report_title"
KEY_ANONYMISE="anonymise"
KEY_HTML_COLOURS="html_colours"
KEY_REPORT_TEMPLATE="report_template"

KEY_QUERY_SUMMARY_DATA="query_summary_data"
KEY_FASTA_SUMMARY_PASS="fasta_summary_pass"
KEY_FASTA_SUMMARY_FAIL="fasta_summary_fail"

# TABLE KEYS
KEY_MUTATIONS="mutations"
KEY_QUERY_TABLE_CONTENT="query_table_content"
KEY_FASTA_TABLE_CONTENT="fasta_table_content"

# TREE KEYS
KEY_TREE_ANNOTATIONS="tree_annotations"
KEY_MAX_TREE_SIZE="max_tree_size"

# TIMELINE KEYS
KEY_TIMELINE_DATES="timeline_dates"
KEY_TIMELINE_GROUP_COLUMN="timeline_group_column"

# BACKGROUND MAP KEYS
KEY_BACKGROUND_MAP_COLUMN="background_map_column"
KEY_BACKGROUND_MAP_LOCATION="background_map_location"
KEY_BACKGROUND_MAP_FILE="background_map_file"
KEY_CENTROID_FILE="centroid_file"
KEY_BACKGROUND_MAP_COLOURS="background_map_colours"
KEY_BACKGROUND_MAP_OTHER_COLOURS="background_map_other_colours"

# QUERY MAP KEYS
KEY_QUERY_MAP_FILE="query_map_file"
KEY_LATITUDE_COLUMN="latitude_column"
KEY_LONGITUDE_COLUMN="longitude_column"
KEY_BACKGROUND_MAP_DATE_RANGE="background_map_date_range"

# QUERY TIME SERIES
KEY_SERIES_COLOUR_FACTOR="series_colour_factor"

# GENERATE BACKGROUND DATA
KEY_GENERATE_CIVET_BACKGROUND_DATA = "generate_civet_background_data"
INPUT_SEQ_OPTIONS = ["gisaid","auspice_source_fasta","fasta"]
INPUT_METADATA_OPTIONS = ["auspice_source_tsv","csv"]

GISAID_HEADER_METADATA_FIELDS = ["sequence_name","gisaid_id","date","country"]
AUSPICE_SOURCE_HEADER_METADATA_FIELDS = ["sequence_name","country"]
AUSPICE_TSV_FIELDS_TO_PARSE = ["strain","gisaid_epi_isl","date","region","country"]
AUSPICE_ACKNOWLEDGEMENTS = ["strain","originating_lab","submitting_lab","authors"]

# MISC KEYS
KEY_VERBOSE="verbose"
KEY_THREADS="threads"
KEY_LOG_API="log_api"
KEY_LOG_STRING="log_string"
KEY_QUIET="quiet"
KEY_CIVET_MODE="civet_mode"
KEY_DATE="date"
KEY_AUTHORS="authors"

KEY_REFERENCE_SEQUENCE="reference_sequence"

RESOURCE_KEY_FILENAME="filename"
RESOURCE_KEY_DIRECTORY="directory"
RESOURCE_KEY="KEY"

# HEADER KEYS
KEY_QUERY_CSV_HEADER="query_csv_header"

PROTECTED_COL_NAMES = ["hash","catchment","query_boolean","qc_status","source","seq_N_content","seq_length","in_tree"]

# background metadata fields
QC_FIELDS = ["N_count","proportion_N","seq_length","QC_status"]
VALUE_PRIMARY_METADATA_FIELDS = ["sequence_name","gisaid_id","sample_date"]
VALUE_SECONDARY_METADATA_FIELDS = ["virus","country","sequence_id","year"]

# DEPENDENCIES AND RESOURCES TO CHECK
dependency_list = ["gofasta","minimap2","snakemake","iqtree","jclusterfunk","scorpio","constellations"]
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