![](./doc_figures/website_header.png)

# Input background data

<strong>civet</strong> summarises information around a set of sequences of interest. It relies on the user providing a background tree, alignment and metadata file. 

The data files <strong>civet</strong>  looks for are:
```
cog_global_2020-XX-YY_alignment.fasta
cog_global_2020-XX-YY_metadata.csv
cog_global_2020-XX-YY_tree.nexus
```

### --CLIMB

For SARS-CoV-2, this data is hosted on CLIMB as part of COG-UK. To run <strong>civet</strong>  on CLIMB with the latest data, either
1) Use the ``--CLIMB`` flag 
or
2) Specify ``CLIMB: True`` in the config.yaml file 

This provides <strong>civet</strong>  with the path to the latest data on CLIMB and allows the user to access adm2 information. 

### -r / --remote

Alternatively, run <strong>civet</strong>  remotely from CLIMB with 
1) The ``-r / --remote`` flag 
or
2) By adding ``remote: True`` to the config file. 

If SSH keys are configured, simply run:

```
civet -i input.csv -r 
```
Otherwise, provide a climb username with ``-uun / --username``:
```
civet -i input.csv -r -uun climb-covid19-smithj
```

Notes:
- This data will access a version of the COG-UK data that is publically available (does not contain adm2 information)
- To access CLIMB in this way, you must have a valid COG-UK CLIMB username and be in the UK

By default, the data will be pulled down to a directory called ``civet-cat`` in the current working directory. 

### -d / --datadir

The user can specify a custom background data directory with the ``-d / --datadir`` flag. 

This can be used with the `remote` option to rsync to an alternative location or without the without the remote flag, <strong>civet</strong>  can just accept the data in that directory as input background data. 

This can also be run on CLIMB without the --CLIMB flag to specify an older version of the dataset. 

```
civet -i input.csv -d path/to/data_directory 
```
### --background-metadata

By default, civet will look for a csv containing background data in the data directory. However, to provide custom background data, use 
1) the ``--background-metadata`` flag
or
2) add `background_metadata: path/to/metadata.csv` to the config file


### Background metadata requirements


The following fields must be **always** present in this background metadata, or civet will not run:

- **sequence_name** containing names of every sequence
- **country** containing the country of sampling
- A field to match the input data with containing COG IDs. The default header for this column is set to **central_sample_id**, but this can be changed by altering the ``--data-column`` argument
- A date column containing the date of sampling. The default header for this column is set to **sample_date**, but can be changed by altering the ``--database-sample-date-column`` argument.

The other fields depend on what you have supplied as optional in the report:
- **node_summary:** if you have altered this argument, the column header in the background metadata that this has been set to must be supplied.
- If the ``--local-lineages`` analysis is being undertaken, **adm2** must be provided.


Some data can be provided in either the query csv or the background metadata:

- Any field you include in ``--table-fields``, ``--label-fields`` or ``--tree-fields``. 
- If the default options for the ``--colour-by`` argument are used, **adm1** must be provided in either the background metadata or the query csv.
- If the default options for ``--table-fields`` argument are used, **uk_lineage, lineage, phylotype** must be included in the background metadata or the query csv.


### [Next: Usage](./usage.md)
