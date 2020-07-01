import geopandas as gp
import pandas as pd
import sjoin
from libpysal.weights import Queen, attach_islands, DistanceBand, set_operations
import maup
import datetime as dtime
import json



##GADM data
gadmpkg="https://github.com/SRooke/GeoFiles/raw/master/gadm36_GBR.gpkg"
##adm1.5 regions as NHS boards/spatial units
england_boards="https://github.com/SRooke/GeoFiles/raw/master/NHS_England_Region_localOffices.geojson"
scotland_boards="https://github.com/SRooke/GeoFiles/raw/master/ScottishHealthBoards.json"
wales_boards="https://github.com/SRooke/GeoFiles/raw/master/WalesHealthBoards.geojson"
NI_all="https://github.com/SRooke/GeoFiles/raw/master/NIall.geojson"
##other county boundaries
england_counties="https://github.com/SRooke/GeoFiles/raw/master/England_counties.geojson"
england_metroCounties="https://raw.githubusercontent.com/SRooke/GeoFiles/master/England_metrocounties.geojson"
NI_counties_data="https://github.com/SRooke/GeoFiles/raw/master/NI_historicCounties.geojson"

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


def dateRestriction(DFin, start_date, end_date):
    DFin['sample_date']=pd.to_datetime(DFin['sample_date'], format="%Y-%m-%d")
    start_date=dtime.datetime.strptime(start_date, "%Y-%m-%d")
    end_date=dtime.datetime.strptime(end_date, "%Y-%m-%d")
    datemask = (DFin['sample_date'] > start_date) & (DFin['sample_date'] <= end_date)
    DFout=DFin.loc[datemask]
    return DFout

def uk_lineage_json(HBCODE,cog_mainland_daterestricted):
    lineage_uk=cog_mainland_daterestricted.loc[cog_mainland_daterestricted['HBCode'] == HBCODE]['uk_lineage'].value_counts()[:10].to_frame()
    lineage_global=cog_mainland_daterestricted.loc[cog_mainland_daterestricted['HBCode'] == HBCODE]['lineage'].value_counts()[:10].to_frame()
    lineage_uk=lineage_uk.reset_index()
    lineage_global=lineage_global.reset_index()
    lineage_uk.columns=['UK_Lineage','Count']
    #lothian_lineage_uk['UK_Lineage']=lothian_lineage_uk['UK_Lineage'].str.replace('UK','')
    lineage_global.columns=['UK_Lineage','Count']
    lineage_uk_json=lineage_uk.to_json(orient='records')
    lineage_global_json=lineage_global.to_json(orient='records')
    return lineage_uk_json

def update_adm15(combinedMAP):
    ##Push all names to single column from mismatched regional shapes
    combinedMAP['HBName'].update(combinedMAP['nhsrlo19nm'])
    combinedMAP['HBName'].update(combinedMAP['lhb19nm'])
    combinedMAP['HBCode'].update(combinedMAP['nhsrlo19cd'])
    combinedMAP['HBCode'].update(combinedMAP['lhb19cd'])
    combinedMAP=combinedMAP.drop(columns=['nhsrlo19cd','nhsrlo19nm','lhb19cd','lhb19nm','lhb19nmw','Shape_Leng','Shape_Area', 'bng_e','bng_n', 'objectid','st_areashape','st_lengthshape'])
    return combinedMAP

def subMapExtractor(loc_code, neighborW, wholeMap):
    ###Takes central location ID, the weighted neighbors, and the full map to extract from
    ## need to extend the neighbors to include self-referencing code
    loc_codeList=[loc_code]
    loc_codeList.extend(neighborW.neighbors[loc_code])
    subMap=wholeMap.query('HBCode in @loc_codeList')
    return subMap

def getIslands(neighborW):
    ## Pull islands for specific joining to the mainland as a whole reference unit
    islandList=neighborW.islands
    return islandList


def getForthBridge(DFin, weightedMatrix):
    ##Fix to ensure Fife and Lothian link via bridge
    landbridge=DistanceBand.from_dataframe(DFin.loc[DFin['HBCode'].str.startswith('S')], 0.35, ids='HBCode')
    w_fixed_setUnion=set_operations.w_union(weightedMatrix,landbridge)
    return w_fixed_setUnion

def mappie():
    multimappie_blank = """{
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "width": 1000,
    "height": 1500,
    "projection": {"type":"mercator"},
    "layer": [
      {
          "data": {
          "values":""
          },

        "mark": {"type":"geoshape", "fill": "#D6D6D6", "stroke": "white"}
    },
      {
          "data": {
          "values":""
          },

        "mark": {"type":"geoshape", "fill": "#112B31", "stroke": "white"}
    }

    ],
    "resolve": {"scale":{"color":"independent", "radius":"independent", "scale":"independent", "theta":"independent"} }
    }
    """
    return multimappie_blank

def multimappie():
    multimappie_blank = """{
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "width": 1000,
    "height": 1500,
    "projection": {"type":"mercator"},
    "layer": [
      {
          "data": {
          "values":""
          },

        "mark": {"type":"geoshape", "fill": "#112B31", "stroke": "white"}
    },
      {
          "data": {
          "values":""
          },

        "mark": {"type":"geoshape", "fill": "#D6D6D6", "stroke": "white"}
    }

    ],
    "resolve": {"scale":{"color":"independent", "radius":"independent", "scale":"independent", "theta":"independent"} }
    }
    """
    return multimappie_blank

def blankJSON():
    blank_json = """{
      "data": {
        "values": ""
      },

      "layer": [
        {
        "mark": {"type": "arc", "innerRadius": 20, "cornerRadius":0.5,"padAngle":0.05, "stroke": "white","strokeWidth":0.5, "thetaOffset":-1.5},
        "encoding": { 
        "text": {"field": "UK_Lineage", "type": "nominal"}
            }
        },
            {
        "mark": {"type": "text", "radiusOffset": 30, "radius":30, "fill":"black", "fontWeight":"bold", "fontSize":15, "font":"Helvetica Neue","angle":0,"align":"center", "thetaOffset":-1.5, "stroke":"white", "strokeWidth":0.1},
        "encoding": { 
          "text": {"field": "UK_Lineage", "type": "nominal"}
            }
        }
    ],
      "encoding": {
        "theta": {"field": "Count", "type": "quantitative", "sort":"descending", "stack": true},
        "radius": {"field": "Count", "type": "quantitative", "scale": {"type": "sqrt", "zero": true, "range": [10, 100]}},
        "color": {"field": "UK_Lineage", "sort":"ascending", "type": "nominal", "scale":{"range": ["#344966", "#E0BE15", "#7B517B", "#5374A2", "#5B8E7D", "#B34248", "#023C40", "#0D1821", "#69626D","#C2C2C2"]}},
          "latitude" : {"datum":""},
          "longitude":{"datum":""}
            }
    }"""
    return blank_json

def tabulateLins(HBCODE, DF):
    tableList=[]
    for each in HBCODE:
        if len(DF.loc[DF['HBCode']==each]['uk_lineage'].value_counts().to_frame()) > 20:
            uk_lin_DF=DF.loc[DF['HBCode'] == each]['uk_lineage'].value_counts[:20].to_frame()
            uk_lin_DF=uk_lin_DF.reset_index()
            uk_lin_DF.columns=['UK Lineage', 'Count']
        else:
            uk_lin_DF=DF.loc[DF['HBCode'] == each]['uk_lineage'].value_counts.to_frame()
            uk_lin_DF=uk_lin_DF.reset_index()
            uk_lin_DF.columns=['UK Lineage', 'Count']
        tableList
        if len(DF.loc[DF['HBCode']==each]['lineage'].value_counts().to_frame()) > 10:
            uk_lin_DF=DF.loc[DF['HBCode'] == each]['lineage'].value_counts[:10].to_frame()
            uk_lin_DF=uk_lin_DF.reset_index()
            uk_lin_DF.columns=['Global Lineage', 'Count']
        else:
            uk_lin_DF=DF.loc[DF['HBCode'] == each]['lineage'].value_counts.to_frame()
            uk_lin_DF=uk_lin_DF.reset_index()
            uk_lin_DF.columns=['Global Lineage', 'Count']


def mapProduce(HBSet, DF, region):
    mapPIE=json.loads(mappie())
    multimapPIE=json.loads(multimappie())
    if len(HBSet) == 1:
        code=HBSet['HBCode'][0]
        lat=HBSet.iloc[0].centroid.y.item()
        long=HBSet.iloc[0].centroid.x.item()
        tempJSON=blankJSON()
        tempJSON=json.loads(tempJSON)
        data=uk_lineage_json(code, DF)
        tempJSON['data']['values']=data
        tempJSON['encoding']['longitude']['datum']=long
        tempJSON['encoding']['latitude']['datum']=lat
        mapPIE['layer'].append(tempJSON)
        mapPIE['layer'][0]['data']['values']=region
        mapPIE['layer'][1]['data']['values']=HBSet['features']
    elif len(HBSet) > 1:
        for row, frame in HBSet.iterrows():
            code=frame['HBCode']
            lat = HBSet.loc[HBSet.index == row].centroid.y.item()
            long = HBSet.loc[HBSet.index == row].centroid.x.item()
            tempjson = blankJSON()
            tempjson = json.loads(tempjson)
            data = uk_lineage_json(code, DF)
            tempjson['data']['values'] = data
            tempjson['encoding']['longitude']['datum'] = long
            tempjson['encoding']['latitude']['datum'] = lat
            multimapPIE['layer'].append(tempjson)
        mapPIE['layer'][0]['data']['values'] = HBSet['features']
        mapPIE['layer'][1]['data']['values'] = region

def getSampleData_final(MetadataLoc, dates=None):
    cog_meta=pd.read_csv(f'{MetadataLoc}')
    cog_meta['central_sample_id'] = cog_meta['sequence_name'].str.split('/', expand=True)[1]
    cog_meta = cog_meta.loc[cog_meta['country'] == 'UK']
    cog_meta_clean=adm2cleaning(cog_meta)
    if dates is not None:
        start=dates['start']
        end=dates['end']
        cog_meta_clean = dateRestriction(cog_meta_clean, start, end)
    return cog_meta_clean

def central_surrounding_regions(HBCODE, neighborW, MAP):
    subregion=subMapExtractor(HBCODE,neighborW,MAP)
    central=subregion.loc[subregion['HBCode'] == HBCODE]
    neighbors=subregion.loc[subregion['HBCode'] != HBCODE]
    return central, neighbors, subregion

##### Start Code Execution #####

######## Handle Maps  ########

##GADM definitions - importing from geopackages
gadm3=gp.read_file(gadmpkg, layer="gadm36_GBR_3")
gadm3['NAME_2']=gadm3['NAME_2'].str.upper()
gadm3['NAME_3']=gadm3['NAME_3'].str.upper()

##English Metropolitan counties
england_metro_counties=gp.read_file(england_metroCounties)

##Historic NI counties
NI_counties=gp.read_file(NI_counties_data)

##Generate list of regions
gadm3_names=gadm3['NAME_3'].str.upper().to_list()
NI_county_names=NI_counties['CountyName'].to_list()
england_MCs_names=england_metro_counties['MCTY18NM'].str.upper().to_list()

###Generate the correct spatial joins for adm2 centroids to 'adm1.5'; NHS health boards for England and Scotland and Wales
england_adm2toboards=adm1_5_assignment(england_boards, gadm3.loc[gadm3['NAME_1'] == "England"])
england_MCstoboards=adm1_5_assignment(england_boards, england_metro_counties)
scotland_adm2toboards=adm1_5_assignment(scotland_boards, gadm3.loc[gadm3['NAME_1'] == 'Scotland'])

###Generate simple translation dictionary for adm2 -> respective health boards
england_boardTranslation=dict(dict(zip(england_adm2toboards.NAME_2,england_adm2toboards.nhsrlo19nm)),**dict(zip(england_MCstoboards.MCTY18NM,england_MCstoboards.nhsrlo19nm)))
scotland_boardTranslation=dict(zip(scotland_adm2toboards.NAME_2, scotland_adm2toboards.HBName))

##Ensuring all locations are upper-case for consistency
mainland_boardTranslation=dict(england_boardTranslation, ** scotland_boardTranslation)
mainland_boardTranslation={k.upper():v for k,v in mainland_boardTranslation.items()}

### Import and assign consistent CRS to all map sub-units
scot_boards=gp.read_file(scotland_boards)
scot_boards=scot_boards.to_crs("EPSG:4326")
eng_boards=gp.read_file(england_boards)
eng_boards=eng_boards.to_crs("EPSG:4326")
NI_region=gp.read_file(NI_all)
NI_region=NI_region.to_crs("EPSG:4326")
wal_boards=gp.read_file(wales_boards)
wal_boards=wal_boards.to_crs("EPSG:4326")

##Concat all regional boards together
mainlandBoards=gp.GeoDataFrame(pd.concat([scot_boards,eng_boards,wal_boards])).reset_index()

### Elminate small-scale gaps between polygons to ensure contiguity calculations work
gapclosedMainland=maup.close_gaps(mainlandBoards['geometry'])

###Get final regional boards map
mainland_boards=gp.GeoDataFrame(mainlandBoards)
mainland_boards.geometry=gapclosedMainland
mainland_boards=update_adm15(mainland_boards)

###Get contiguity neighbors for mainland
mainland_boards_W=Queen.from_dataframe(mainland_boards, idVariable='HBCode')
##Hacky fix to link Fife and Lothian
mainland_boards_final_W=getForthBridge(mainland_boards,mainland_boards_W)


HBCode_translation=dict(zip(mainland_boards.HBName, mainland_boards.HBCode))




cog_restricted=getSampleData_final(DATALOC)

###Preparing cog_meta for spatial filtering & processing
cog_meta_mainland=cog_restricted.loc[cog_restricted['adm1'].isin(['UK-SCT', 'UK-ENG', 'UK-WAL'])]
cog_meta_mainland['HBName']=cog_meta_mainland.loc[:, 'adm2'].map(mainland_boardTranslation)
cog_meta_mainland['HBCode']=cog_meta_mainland.loc[:, 'HBName'].map(HBCode_translation)


central, neighboring, submap = central_surrounding_regions(HB_code, mainland_boards_final_W, mainland_boards)



