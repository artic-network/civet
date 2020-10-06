
<section id="banner">
    <div class="content">
      <header>
        <h2>Background data</h2>
        <p>How to access global background data</p>
      </header>
    </div>
    <span class="image object">
        <img src="./figures/civet_logo.png" alt="" style="max-width:150px"/>
        </span>
</section>


### Input background data

<strong>civet</strong> summarises information around a set of sequences of interest. It relies on the user providing a background tree, alignment and metadata file. 

The data files <strong>civet</strong>  looks for are:
```
cog_global_alignment.fasta
cog_global_metadata.csv
cog_global_tree.nexus
```

<strong>--CLIMB</strong>

For SARS-CoV-2, this data is hosted on CLIMB as part of COG-UK. To run <strong>civet</strong>  on CLIMB with the latest data, use the ``--CLIMB`` flag (or specify ``CLIMB: True`` in the config file) command to <strong>civet</strong> . This provides <strong>civet</strong>  with the path to the latest data on CLIMB. 

<strong>-r / --remote</strong>

Alternatively, run <strong>civet</strong>  remotely from CLIMB with the ``-r / --remote`` flag or by adding ``remote: True`` to the. If SSH keys are configured, simply run:

```
civet -i input.csv -r 
```
Otherwise, provide a climb username with ``-uun / --username``:
```
civet -i input.csv -r -uun climb-covid19-smithj
```

By default, the data will be pulled down to a directory called ``civet-cat`` in the current working directory. 
<strong>-d / --datadir</strong>

The user can specify a custom directory with the ``-d / --datadir`` flag. 

If the user has background data locally (or would like to run on CLIMB with a different version of the data), by specifying ``--datadir`` without the remote flag, <strong>civet</strong>  can just accept the data in that directory as input background data.

```
civet -i input.csv -d path/to/data_directory 
```

By default, civet will look for a csv containing background data in the data directory. However, to provide custom background data, use the ``--background-metadata`` flag. 

The following fields must be **always** present in this background metadata, or civet will not run properly:

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





