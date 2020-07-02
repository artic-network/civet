if __name__ == "__main__":
    main()

import geopandas as gp
import pandas as pd
from libpysal.weights import Queen, attach_islands, DistanceBand, set_operations
import datetime as dtime
import json
import pickle
from vega import VegaLite
from IPython.display import display,display_png
from IPython.display import display_html
def VegaLite(spec):
    bundle = {}
    bundle['application/vnd.vegalite.v3+json'] = spec
    display(bundle, raw=True)

def adm2cleaning(data_cog):
    data_cog2=data_cog.dropna(subset=['adm2']).copy()
    data_cog2['adm2']=data_cog['adm2'].str.replace('^MOTHERWELL$','NORTH LANARKSHIRE', regex=True).copy()
    data_cog2['adm2']=data_cog['adm2'].str.replace('^GREATER_LONDON$','GREATER LONDON', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^BORDERS$','SCOTTISH BORDERS', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^PAISLEY$','RENFREWSHIRE', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^SUTTON-IN-ASHFIELD$','NOTTINGHAMSHIRE', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^LONDON$','GREATER LONDON', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^KILMARNOCK$','EAST AYRSHIRE', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^PERTH$','PERTHSHIRE AND KINROSS', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^ISLE OF ANGLESEY$','ANGLESEY', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^SHETLAND$','SHETLAND ISLANDS', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^GREATER LONDONDERRY$','LONDONDERRY', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^STAFFORD$','STAFFORDSHIRE', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^INVERNESS$','HIGHLAND', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^KINGS NORTON$','WEST MIDLANDS', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^GLASGOW CITY$','GLASGOW', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^CITY OF LONDON$','GREATER LONDON', regex=True)
    data_cog2['adm2']=data_cog['adm2'].str.replace('^RHONDDA CYNON TAF$','RHONDDA, CYNON, TAFF', regex=True)
    return data_cog2


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

singlemappie_blank = '''{
  "$schema": "https://vega.github.io/schema/vega-lite/v4.json",
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
    '''

multimappie_blank="""{
"$schema": "https://vega.github.io/schema/vega-lite/v4.json",
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


blank_json="""{
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

def tabulateLins(HBCODE, DF, HBNAME):
    #lineage_uk=cog_mainland_daterestricted.loc[cog_mainland_daterestricted['HBCode'] == HBCODE]['uk_lineage'].value_counts()[:10].to_frame()
    for each in [HBCODE]:
        if len(DF.loc[DF['HBCode']==each]['uk_lineage'].value_counts().to_frame()) > 20:
            uk_lin_DF=DF.loc[DF['HBCode'] == each]['uk_lineage'].value_counts()[:20].to_frame()
            uk_lin_DF=uk_lin_DF.reset_index()
            uk_lin_DF.columns=['UK Lineage', 'Count']
        else:
            uk_lin_DF=DF.loc[DF['HBCode'] == each]['uk_lineage'].value_counts().to_frame()
            uk_lin_DF=uk_lin_DF.reset_index()
            uk_lin_DF.columns=['UK Lineage', 'Count']
        if len(DF.loc[DF['HBCode']==each]['lineage'].value_counts().to_frame()) > 10:
            glob_lin_DF=DF.loc[DF['HBCode'] == each]['lineage'].value_counts()[:10].to_frame()
            glob_lin_DF=glob_lin_DF.reset_index()
            glob_lin_DF.columns=['Global Lineage', 'Count']
        else:
            glob_lin_DF=DF.loc[DF['HBCode'] == each]['lineage'].value_counts().to_frame()
            glob_lin_DF=glob_lin_DF.reset_index()
            glob_lin_DF.columns=['Global Lineage', 'Count']
    tableList=[uk_lin_DF, glob_lin_DF]
    print(HBNAME)
    print(uk_lin_DF.to_markdown())
    print(glob_lin_DF.to_markdown())

    return tableList


def mapProduce(HBSet, DF, region):
    mapPIE=json.loads(multimappie_blank)
    #multimapPIE=json.loads(multimappie())
    HBSet_json = HBSet.to_json()
    HBSet_json = json.loads(HBSet_json)
    region_json = region.to_json()
    region_json = json.loads(region_json)
    #print(HBSet.features)
    if len(HBSet) == 1:
        for row, frame in HBSet.iterrows():
            code=frame['HBCode']
            lat=HBSet.loc[HBSet.index == row].centroid.y.item()
            long=HBSet.loc[HBSet.index == row].centroid.x.item()
            tempJSON=json.loads(blank_json)
            data=uk_lineage_json(code, DF)
            #print(data)
            tempJSON['data']['values']=data
            tempJSON['encoding']['longitude']['datum']=long
            tempJSON['encoding']['latitude']['datum']=lat
            mapPIE['layer'].append(tempJSON)
        mapPIE['layer'][0]['data']['values']=region_json['features']
        mapPIE['layer'][1]['data']['values']=HBSet_json['features']
        VegaLite(mapPIE)
        return mapPIE
    elif len(HBSet) > 1:
        for row, frame in HBSet.iterrows():
            code=frame['HBCode']
            lat = HBSet.loc[HBSet.index == row].centroid.y.item()
            long = HBSet.loc[HBSet.index == row].centroid.x.item()
            tempJSON=json.loads(blank_json)
            data=uk_lineage_json(code, DF)
            tempJSON['data']['values']=data
            tempJSON['encoding']['longitude']['datum']=long
            tempJSON['encoding']['latitude']['datum']=lat
            mapPIE['layer'].append(tempJSON)
        mapPIE['layer'][0]['data']['values'] = HBSet_json['features']
        mapPIE['layer'][1]['data']['values'] = region_json['features']
        return mapPIE

def getSampleData(metadataLoc):
    cog_meta = pd.read_csv(metadataLoc)
    cog_meta['sample_date'] = pd.to_datetime(cog_meta['sample_date'], format="%Y-%m-%d")
    return cog_meta


def getSampleData_final(MetadataDF, dates=None):
    cog_meta=MetadataDF
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

def adm2_to_centralHBCode(adm2, translation_dict, HbtoCode):
    HBs=[]
    #print(adm2)
    for each in adm2:
        if each in translation_dict:
            HBs.append(translation_dict[each])
    if len(HBs) > 0:
        centralHB=max(HBs,key=HBs.count)
        centralHBCode=HbtoCode[centralHB]
        #[k for k, v in a.items() if v == b]
        print(centralHBCode)
        return centralHBCode
    elif len(HBs) == 0:
        return None

def defineDateRestriction(samplesDF, specifiedRange):
    datedSamples=samplesDF.drop_na(subset=['Collection_Date'])
    if len(datedSamples) == 0:
        print('No collection dates specified, will revert to using all available data for local lineage analysis.')
        return None
    if  len(datedSamples) == 1:
        daterange=''
        return daterange






##### Start Code Execution #####
#mainland_boards=gp.read_file(f'{data}/maps/Mainland_HBs_gapclosed_simplified.geojson')
#mainland_boards=update_adm15(mainland_boards)
###Get contiguity neighbors for mainland
#mainland_boards_W=Queen.from_dataframe(mainland_boards, idVariable='HBCode')
##Hacky fix to link Fife and Lothian
#mainland_boards_final_W=getForthBridge(mainland_boards,mainland_boards_W)
###Create translation dict for board codes
#HBCode_translation=dict(zip(mainland_boards.HBName, mainland_boards.HBCode))
##########Data loading###########
##Get HB translation dict##
#HBTranslation=pickle.load(f'{data}/maps/HB_Translation.pkl')
###Checking user defined dates###
#if DATERANGE is not None:
#    cog_restricted = getSampleData_final(f'{metadata}', DATERANGE)
#else:
#    cog_restricted=getSampleData_final(f'{metadata}')
###Preparing cog_meta for spatial filtering & processing
#cog_meta_mainland=cog_restricted.loc[cog_restricted['adm1'].isin(['UK-SCT', 'UK-ENG', 'UK-WAL'])]
#cog_meta_mainland['HBName']=cog_meta_mainland.loc[:, 'adm2'].map(mainland_boardTranslation)
#cog_meta_mainland['HBCode']=cog_meta_mainland.loc[:, 'HBName'].map(HBCode_translation)
### Get submitted sample metadata ###
#sampleDataDF=pd.read_csv(sampleCSV)
#HB_code=adm2_to_centralHBCode(sampleDataDF['adm2'].to_list())
#if HB_code is not None:
#    ## Get the localised regions ##
#    central, neighboring, submap = central_surrounding_regions(HB_code, mainland_boards_final_W, mainland_boards)
#    ## Generate tabular data for each region ##
#    for each in [central, neighboring]:
#        for row, frame in each.iterrows():
#            print(frame['HBName'])
#            tabulateLins(frame['HBCode'], cog_meta_mainland)
#    ## Generated Mapping ##
#    for each, region in [central, neighboring], ['Central Health Board', 'Neighboring Health Boards']:
#        print (region)
#        mapProduce(each, cog_meta_mainland, submap)
