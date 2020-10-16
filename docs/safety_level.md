
## Data anonymisation

adm2 & COG-IDs, CLIMB and local civet. This section clarifies what data is where and what settings of civet let the user remain in line with the data access agreement.


### Running locally vs running on CLIMB

When runing civet with the `--CLIMB` flag, the user can access adm2 information. This information must not leave CLIMB and there is COG-ID masking that takes place with the default civet options. 

When running civet on a local machine with the `-r` flag, only publically available data is rsynced from CLIMB. This means that running civet on your local machine will not allow you to access COG-wide adm2. If you have adm2 for your own queries (or outer postcode or other more fine-grain geographic information), this can be supplied for your queries in the input.csv file.

### --safety-level and data masking

`--safety-level` is an integer value (0, 1 or 2) that configures the level of anonymisation needed in the report. Default: 1. When running on CLIMB specify either 1 or 2 to decide between accessing background COG IDs or background adm2 information. 

0) Anonymisation off. COG-IDs and adm2 metadata shown for all sequences if adm2 present in either query or background metadata. This information should not be copied or disseminated 
1) No COG-IDs shown for anything but query sequences specified with the `-i / --input` flag (i.e. not with the `--from-metadata` option). If run on CLIMB, adm2 present in the report. 
2) COG-IDs shown on query sequenecs and background sequences, but none show adm2 (equivalent to the 'publically available' data). 

The query sequences, if they are specified as a list of sequences or in a csv, will always show either their name or the specified `--display-name` regardless of safety level.

<strong>--label-fields</strong> 
A comma separated string of metadata headers containing information to be displayed in the labels in the phylogeny.

By default, the `input_column` and `sample_date` (if provided) will be used. Depending on the --safety-level and if adm2 is present in either metadata (for example if you ran on CLIMB) this will also be shown. 

If any fields are provided using the label fields option, these will be overwritten, and the display name and specified fields will be shown. 

For example, a query sequence has the following metadata:

| display name | adm1 | adm2 | date | care home | 
| --- | --- | --- | --- | --- |
| EDB1234 | Scotland | Edinburgh | 2020-09-23 | A | 

Under safety level 0 and safety level 1, the default label would be:
"EDB1234|Edinburgh|2020-09-23"

Under safety level 2, the default label would be:
"EDB1234|2020-09-23"

and using  
```--label-fields adm1,care_home``` 
regardless of safety level will show:
"EDB1234|Scotland|A"

Therefore once you specify *any* label fields, you must specify all those that you want other than the display name.

**Note: the label-fields overwrites the safety level for all query data. If you specify adm2, regardless of safety level, the query sequences will show adm2 in their labels**. However, the appropriate anonymisation of the COG IDs will still take place if safety level is 1 and the queries were found using `--from-metadata`, and background sequence labels are assigned by safety level.
