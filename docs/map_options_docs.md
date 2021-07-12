![](./doc_figures/website_header.png)


## Mapping options
 
There are two types of mapping in civet:
1. Plot background diversity of lineages in the area of interest (report option 6)
2. Plot where your queries are (report option 7) 
 
Both maps will show a default first view, but it is possible to alter that view by using the mouse wheel to zoom in and out, and by clicking and dragging on the map itself.
 
### Mapping queries (report option 6)
 
This option allows you to plot queries on a map using their coordinates.
 
Arguments:
--query-map-file/-qmfile : the underlying map file that the dots will be plotted on. Default: world map
--latitude-column/-lat : column in the metadata that contains latitude information for queries. Default: latitude
--longitude-column/-long : column in the metadata that contains longitude information for queries. Default: longitude
 
### Mapping background diversity (report option 7)
 
This option summarises the PANGO lineages as radial charts during a set time period at a specific geographical level. Only lineages that make up more than 5% of the overall sequence count are shown, with others shown under “other”.
 
For more information on PANGO lineages, see https://www.pango.network/
 
Arguments:
--background-map-file/-bmfile : the underlying map file that the radial charts will be displayed on. 
Default: world map
--centroid-file : file containing centroids of locations that lineages are summarised by, used to place radial plots on the map. Only necessary if a custom background map file is provided.
--background-map-date-range/-bmrange : Date range to restrict summary to. This can either be in terms of two dates (inclusive) in the format YYYY-MM-DD separated by a “:” (eg) or a single integer specifying the number of days either side of the date range of the queries. Eg --background-map-date-range 2020-09-10:2020-09-30 or --background-map-date-range 10 
Default: 2019-11-01 to present day
--background-map-column/-bmcol : column in the background dataset containing geographical information to summarise PANGO lineages by. Must correspond to locations in the centroid file and background map file. 
Default: --location-column (country by default)
--background-map-location/-bmloc : comma separated list of locations to show background lineage diversity for. Must be valid locations matching the map file and the centroid file. 
Default: all valid locations in the background metadata file in the appropriate column.
 
NB if the environment variable civet_mode is set to CLIMB, the default for the background-map-column is “suggested_adm2_grouping” instead, with the corresponding map and centroid files. This is an output of the COG-UK data processing pipeline, and is a version of the administrative level 2 (roughly county) of where the sequence was sampled. See https://github.com/COG-UK/docs/blob/master/geography_cleaning.md for more information.


### [Next: Acknowledgements](./acknowledgements.md)
