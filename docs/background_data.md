
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