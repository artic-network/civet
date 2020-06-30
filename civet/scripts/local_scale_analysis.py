import geopandas as gp
import pandas as pd
import sjoin
from libpysal.weights import Queen, KNN, attach_islands, DistanceBand, set_operations
import altair as alt
import maup
import datetime as dtime


##adm1.5 regions as NHS boards/spatial units
england_boards="https://github.com/SRooke/GeoFiles/raw/master/NHS_England_Region_localOffices.geojson"
scotland_boards="https://github.com/SRooke/GeoFiles/raw/master/ScottishHealthBoards.json"
wales_boards="https://github.com/SRooke/GeoFiles/raw/master/WalesHealthBoards.geojson"
NI_all="http://osni-spatialni.opendata.arcgis.com/datasets/c55271de58434719b3248e4e3c2e9bb0_0.geojson"
##other county boundaries
england_counties="https://opendata.arcgis.com/datasets/3fb7653392cc4d0fafecb56d96d22109_0.geojson"
england_metroCounties="https://opendata.arcgis.com/datasets/389f538f35ef4eeb84965dfd7c0a0b47_0.geojson"
UK_localauthorities="https://opendata.arcgis.com/datasets/b6d2e15801de45328b760a4f55d74318_0.geojson"
NI_counties_data="http://osni-spatialni.opendata.arcgis.com/datasets/f1f716b3f1e542e39d36e45438d5e240_1.geojson"





def adm2cleaning(data_cog):
    data_cog.dropna(subset=['adm2'],inplace=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^MOTHERWELL$','NORTH LANARKSHIRE', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^GREATER_LONDON$','GREATER LONDON', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^BORDERS$','SCOTTISH BORDERS', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^PAISLEY$','RENFREWSHIRE', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^SUTTON-IN-ASHFIELD$','NOTTINGHAMSHIRE', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^LONDON$','GREATER LONDON', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^KILMARNOCK$','EAST AYRSHIRE', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^PERTH$','PERTHSHIRE AND KINROSS', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^ISLE OF ANGLESEY$','ANGLESEY', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^SHETLAND$','SHETLAND ISLANDS', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^GREATER LONDONDERRY$','LONDONDERRY', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^STAFFORD$','STAFFORDSHIRE', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^INVERNESS$','HIGHLAND', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^KINGS NORTON$','WEST MIDLANDS', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^GLASGOW CITY$','GLASGOW', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^CITY OF LONDON$','GREATER LONDON', regex=True)
    data_cog['adm2']=data_cog['adm2'].str.replace('^RHONDDA CYNON TAF$','RHONDDA, CYNON, TAFF', regex=True)
    return data_cog


def adm1_5_assignment(data_boundaries_URL, data_adm):
    data_NHSboard = gp.read_file(data_boundaries_URL)
    data_NHSboard = data_NHSboard.to_crs("EPSG:4326")
    data_adm = data_adm.to_crs("EPSG:4326")
    data_adm = data_adm.set_geometry(col=data_adm.geometry.centroid)
    adm1_5 = sjoin.sjoin_nearest(data_adm, data_NHSboard, how='left', search_radius=0.2)
    return adm1_5
cog_meta=pd.read_csv("~/GitKraken/civet/workspace/cog_global_2020-06-19_metadata.csv")
cog_meta['central_sample_id']=cog_meta['sequence_name'].str.split('/', expand=True)[1]
cog_meta=cog_meta.loc[cog_meta['country'] =='UK']
cog_noadm2=cog_meta.loc[~cog_meta['adm2'].notna()]
##this will remove adm2 NAs from DF
cog_meta_clean=adm2cleaning(cog_meta)

###Generate the correct spatial joins for adm2 centroids to 'adm1.5'; NHS health boards for England and Scotland and Wales
england_adm2toboards=adm1_5_assignment(england_boards, gadm3.loc[gadm3['NAME_1'] == "England"])
england_MCstoboards=adm1_5_assignment(england_boards, england_metro_counties)
scotland_adm2toboards=adm1_5_assignment(scotland_boards, gadm3.loc[gadm3['NAME_1'] == 'Scotland'])
###Generate simple translation dictionary for adm2 -> respective health boards
england_boardTranslation= dict(dict(zip(england_adm2toboards.NAME_2, england_adm2toboards.nhsrlo19nm)),**dict(zip(england_MCstoboards.mcty18nm,england_MCstoboards.nhsrlo19nm)))
scotland_boardTranslation = dict(zip(scotland_adm2toboards.NAME_2, scotland_adm2toboards.HBName))

mainland_boardTranslation = dict(england_boardTranslation, ** scotland_boardTranslation)
mainland_boardTranslation = {k.upper(): v for k, v in mainland_boardTranslation.items()}
###Preparing cog_meta for spatial filtering & processing
cog_meta_mainland=cog_meta_clean.loc[cog_meta_clean['adm1'].isin(['UK-SCT', 'UK-ENG', 'UK-WAL'])]
cog_meta_mainland['HBName']=cog_meta_mainland.loc[:, 'adm2'].map(mainland_boardTranslation)
cog_meta_mainland['HBCode']=cog_meta_mainland.loc[:, 'HBName'].map(HBCode_translation)

