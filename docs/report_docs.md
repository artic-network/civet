![](./doc_figures/website_header.png)

# Report options
How to customise the civet report
 
## Report content and report preset
 
Each element of the civet report is associated with a number:
 
Report elements:
1. Tables describing the queries
2. Summarise the phylogenetic context of the queries in table form
3. Show phylogenetic trees of the queries and their surrounding phylogenetic context, known as "catchments"
4. Show "snpit" plot of each catchment, which visualises the SNPs found in each query sequence compared to the reference sequence.
5. Produces timelines of queries per catchment
6. Summarises Pango lineages found in the background dataset at appropriate geographical level (eg country) and displays them on a map.
7. Plots queries as dots on a map
 
Provide combinations of these letters as a comma separated string to define the elements present in the report.
 
You can also use presets using --report-preset/-rp:
 
hold_the_sauce: 1,2 (just summary tables, no figures)
the_usual : 1,2,3,4,5 (All figures except maps) 
the_works : 1,2,3,4,5,7 (All figures except background diversity map)
the_whole_shebang : 1,2,3,4,5,6,7 (everything!)
 
Default: 1,2,3,4,5
 
Usage: 
--report-content {$LIST} or -rc {$LIST}
--report-preset {$PRESET} or -rp {$PRESET}
 
## Report title
 
Title to show at the top of the report.
 
Default: civet report
Usage: --report-title {$TITLE} or -rt {$TITLE}
 
## Colour theme
 
Report theme colour
 
Default: #7178bc
Usage: --colour-theme/-ct {$HEX_CODE}
 
## Colour map
 
Comma separated string of HEX codes or HTML-support names to colour factors in different figures in the report by.
Default: civet logo colours
Usage: --colour-map/-cmap {$LIST}
 
## Report column and anonymise
 
Controls how queries are referred to in the report tables and figures.
 
Either, you may provide a column in the metadata to refer to queries by using --report-column/-rcol (eg can be useful for referring to patient IDs).
 
Alternatively, the --anonymise flag will provide an arbitrary number to refer to the sequence by.
 
Default: False
Usage: 
--report-column{$COLUMN_NAME} or -rcol {$COLUMN_NAME}
--anonymise
 
## location column
 
Column containing location information in the background metadata csv. If present, it will be displayed by default in query and catchment summary tables. 
Default: country
Usage: --location-column/-lcol {$COLUMN_NAME}
 
 
## date column and background date column
Columns in the input query csv and background csv respectively that contain date information. If present, it will be displayed by default in query and catchment summary tables.
 
Default: sample_date
Usage: 
--date-column{$COLUMN_NAME} or -dcol {$COLUMN_NAME}
--background-date-column{$COLUMN_NAME} or -bdate {$COLUMN_NAME}
 
## table content and catchment table
 
Comma separated strings to define what is present in the query summary tables (report option 1) and the catchment summary tables (report option 2).
 
Defaults: 
Table content: --background_column,--background_date_column,source,lineage,country,catchment
Catchment table:
count,country,lineage
Usage: 
--table-content {$LIST}
--catchment-table {$LIST}
 
## Timeline dates and timeline colours (report option 5)
 
Provide a comma separated list of column headers present in the query metadata or background metadata to plot in the timeline, and hex codes or names to colour each dot by. If providing colours, please provide the same number of colours as the number of fields displayed in the timeline.
 
Defaults: 
sample_date or --date-column if provided
 
Usage:
--timeline-dates {$LIST_HEADERS} or -td {$LIST_HEADERS}
--timeline-colours {$LIST_HEXES}







### [Next: Map options and usage instructions](./map_option_docs.md)
