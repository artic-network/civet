#!/usr/bin/env python3
import os
import warnings
import pickle
import geopandas as gp
import pandas as pd
from libpysal.weights import Queen, attach_islands, DistanceBand, set_operations
import json
#from vega import VegaLite
import argparse
warnings.filterwarnings("ignore")
print("running local scale analysis")
# ~~~~~~~ Define input variables ~~~~~~~~~~

parser = argparse.ArgumentParser(description='Parse barcode info and minimap paf file, create report.')
parser.add_argument("--date-restriction", default="False", action="store", type=str, dest="date_restriction")
parser.add_argument("--date-pair-start", default="", action="store", type=str, dest="date_pair_start")
parser.add_argument("--date-pair-end", default="", action="store", type=str, dest="date_pair_end")
parser.add_argument("--cog-meta-global", action="store", type=str, dest="cog_metadata")
parser.add_argument("--user-sample-data", action="store", type=str, dest="user_sample_data")
parser.add_argument("--date-window", action="store",required=False, type=int, dest="date_window")
parser.add_argument("--output-base-dir", action="store", type=str, dest="output_base_dir")
parser.add_argument("--output-temp-dir", action="store", type=str, dest="output_temp_dir")
parser.add_argument("--hb-translation", action="store", type=str, dest="hb_translation")
parser.add_argument("--uk-map", action="store", type=str, dest="uk_map")
#parser.add_argument("--civet-cat", action="store", type=str, dest="civet_cat_dir")

argsIN=parser.parse_args()

currentDir=os.getcwd()
#####    NOTE ------- NEEEDS TO KNOW WHAT THE GENERATED OUTPUT DIR IS
outDIR=argsIN.output_base_dir
#os.path.join(currentDir, 'GENERATEDOUTPUTDIR', 'figures', 'Mapping')

## Needed temporarily to point to correct additional data files (map, pkl)
#civet_dir = "/mnt/e/Users/Stefan/GitKraken/civet"
# Needed to pull metadata
#civet_cat_dir=os.path.join(currentDir, 'civet_cat')

### needs importation from installation of civet
translator=argsIN.hb_translation
mapfile=argsIN.uk_map

date_pair=[]
for each in [argsIN.date_pair_start, argsIN.date_pair_end]:
    if each != "None":
        date_pair.append(each)


# ~~~~~~~~ Defining functions ~~~~~~~~~~
#def VegaLite(spec):
#    bundle = {}
 #   bundle['application/vnd.vegalite.v3+json'] = spec
 #   display(bundle, raw=True)


def adm2cleaning(data_cog, samplecsv=False):
    data_cog2 = ""
    if samplecsv == False:
        data_cog2 = data_cog.dropna(subset=['adm2']).copy()
    if samplecsv == True:
        data_cog2 = data_cog.copy()
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^MOTHERWELL$', 'NORTH LANARKSHIRE', regex=True).copy()
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^GREATER_LONDON$', 'GREATER LONDON', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^BORDERS$', 'SCOTTISH BORDERS', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^PAISLEY$', 'RENFREWSHIRE', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^SUTTON-IN-ASHFIELD$', 'NOTTINGHAMSHIRE', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^LONDON$', 'GREATER LONDON', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^KILMARNOCK$', 'EAST AYRSHIRE', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^PERTH$', 'PERTHSHIRE AND KINROSS', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^ISLE OF ANGLESEY$', 'ANGLESEY', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^SHETLAND$', 'SHETLAND ISLANDS', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^GREATER LONDONDERRY$', 'LONDONDERRY', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^STAFFORD$', 'STAFFORDSHIRE', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^INVERNESS$', 'HIGHLAND', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^KINGS NORTON$', 'WEST MIDLANDS', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^GLASGOW CITY$', 'GLASGOW', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^CITY OF LONDON$', 'GREATER LONDON', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^RHONDDA CYNON TAF$', 'RHONDDA, CYNON, TAFF', regex=True)
    return data_cog2


def dateRestriction(DFin, dateDict):
    DFin['sample_date'] = pd.to_datetime(DFin['sample_date'], format="%Y-%m-%d")
    datemask = (DFin['sample_date'] > dateDict['start_date']) & (DFin['sample_date'] <= dateDict['end_date'])
    DFout = DFin.loc[datemask]
    return DFout


def uk_lineage_json(HBCODE, cog_mainland_daterestricted):
    lineage_uk = cog_mainland_daterestricted.loc[cog_mainland_daterestricted['HBCode'] == HBCODE][
                     'uk_lineage'].value_counts()[:10].to_frame()
    lineage_global = cog_mainland_daterestricted.loc[cog_mainland_daterestricted['HBCode'] == HBCODE][
                         'lineage'].value_counts()[:10].to_frame()
    lineage_uk = lineage_uk.reset_index()
    lineage_global = lineage_global.reset_index()
    lineage_uk.columns = ['UK_Lineage', 'Count']
    # lothian_lineage_uk['UK_Lineage']=lothian_lineage_uk['UK_Lineage'].str.replace('UK','')
    lineage_global.columns = ['UK_Lineage', 'Count']
    lineage_uk_json = lineage_uk.to_json(orient='records')
    lineage_global_json = lineage_global.to_json(orient='records')
    return lineage_uk_json


def update_adm15(combinedMAP):
    ##Push all names to single column from mismatched regional shapes
    combinedMAP['HBName'].update(combinedMAP['nhsrlo19nm'])
    combinedMAP['HBName'].update(combinedMAP['lhb19nm'])
    combinedMAP['HBCode'].update(combinedMAP['nhsrlo19cd'])
    combinedMAP['HBCode'].update(combinedMAP['lhb19cd'])
    combinedMAP = combinedMAP.drop(
        columns=['nhsrlo19cd', 'nhsrlo19nm', 'lhb19cd', 'lhb19nm', 'lhb19nmw', 'Shape_Leng', 'Shape_Area', 'bng_e',
                 'bng_n', 'objectid', 'st_areashape', 'st_lengthshape'])
    return combinedMAP


def subMapExtractor(loc_code, neighborW, wholeMap):
    ###Takes central location ID, the weighted neighbors, and the full map to extract from
    ## need to extend the neighbors to include self-referencing code
    loc_codeList = [loc_code]
    loc_codeList.extend(neighborW.neighbors[loc_code])
    subMap = wholeMap.query('HBCode in @loc_codeList')
    return subMap


def getIslands(neighborW):
    ## Pull islands for specific joining to the mainland as a whole reference unit
    islandList = neighborW.islands
    return islandList


def getForthBridge(DFin, weightedMatrix):
    ##Fix to ensure Fife and Lothian link via bridge
    landbridge = DistanceBand.from_dataframe(DFin.loc[DFin['HBCode'].str.startswith('S')], 0.35, ids='HBCode')
    w_fixed_setUnion = set_operations.w_union(weightedMatrix, landbridge)
    return w_fixed_setUnion

def tabulateLins(HBCODE, DF, HBNAME):
    for each in [HBCODE]:
        if len(DF.loc[DF['HBCode'] == each]['uk_lineage'].value_counts().to_frame()) > 20:
            uk_lin_DF = DF.loc[DF['HBCode'] == each]['uk_lineage'].value_counts()[:20].to_frame()
            uk_lin_DF = uk_lin_DF.reset_index()
            uk_lin_DF.columns = ['UK Lineage', 'Count']
            uk_lin_DF['Count'] = uk_lin_DF['Count'].astype(int)
        else:
            uk_lin_DF = DF.loc[DF['HBCode'] == each]['uk_lineage'].value_counts().to_frame()
            uk_lin_DF = uk_lin_DF.reset_index()
            uk_lin_DF.columns = ['UK Lineage', 'Count']
            uk_lin_DF['Count'] = uk_lin_DF['Count'].astype(int)

        if len(DF.loc[DF['HBCode'] == each]['lineage'].value_counts().to_frame()) > 10:
            glob_lin_DF = DF.loc[DF['HBCode'] == each]['lineage'].value_counts()[:10].to_frame()
            glob_lin_DF = glob_lin_DF.reset_index()
            glob_lin_DF.columns = ['Global Lineage', 'Count']
            glob_lin_DF['Count'] = glob_lin_DF['Count'].astype(int)
        else:
            glob_lin_DF = DF.loc[DF['HBCode'] == each]['lineage'].value_counts().to_frame()
            glob_lin_DF = glob_lin_DF.reset_index()
            glob_lin_DF.columns = ['Global Lineage', 'Count']
            glob_lin_DF['Count'] = glob_lin_DF['Count'].astype(int)

    tableList = [uk_lin_DF, glob_lin_DF]
    # print(HBNAME)
    # [print(each.to_markdown()) for each in tableList]
    MDout = pd.concat(tableList, axis=1).fillna("").to_markdown(showindex=False)
    return HBNAME, MDout


def mapProduce(HBSet, DF, region, centralCode=None):
    singlemapPIE = json.loads(singlemappie_blank)
    multimapPIE = json.loads(multimappie_blank)
    # print(type(mapPIE))
    # multimapPIE=json.loads(multimappie())
    HBSet_json = HBSet.to_json()
    HBSet_json = json.loads(HBSet_json)
    region_json = region.to_json()
    region_json = json.loads(region_json)
    # print(HBSet.features)
    if len(HBSet) == 1:
        for row, frame in HBSet.iterrows():
            code = frame['HBCode']
            lat = HBSet.loc[HBSet.index == row].centroid.y.item()
            long = HBSet.loc[HBSet.index == row].centroid.x.item()
            tempJSON = json.loads(blank_json)
            data = uk_lineage_json(code, DF)
            # print(data)
            tempJSON['data']['values'] = data
            tempJSON['encoding']['longitude']['datum'] = long
            tempJSON['encoding']['latitude']['datum'] = lat
            singlemapPIE['layer'].append(tempJSON)
        singlemapPIE['layer'][0]['data']['values'] = region_json['features']
        singlemapPIE['layer'][1]['data']['values'] = HBSet_json['features']
        #VegaLite(singlemapPIE)
        return singlemapPIE
    elif len(HBSet) > 1:
        for row, frame in HBSet.iterrows():
            code = frame['HBCode']
            lat = HBSet.loc[HBSet.index == row].centroid.y.item()
            long = HBSet.loc[HBSet.index == row].centroid.x.item()
            tempJSON = json.loads(blank_json)
            data = uk_lineage_json(code, DF)
            tempJSON['data']['values'] = data
            tempJSON['encoding']['longitude']['datum'] = long
            tempJSON['encoding']['latitude']['datum'] = lat
            multimapPIE['layer'].append(tempJSON)
        multimapPIE['layer'][0]['data']['values'] = region_json['features']
        multimapPIE['layer'][1]['data']['values'] = HBSet_json['features']
        #VegaLite(multimapPIE)
        return multimapPIE
    if centralCode is not None:
        newDF1 = HBSet.loc[HBSet['HBCode' == centralCode]].copy()
        newDF2 = HBSet.loc[HBSet['HBCode'] != centralCode].copy()
        newDF3 = pd.concat([newDF1, newDF2])
        for row, frame in newDF3.iterrows():
            code = frame['HBCode']
            lat = newDF3.loc[newDF3.index == row].centroid.y.item()
            long = newDF3.loc[newDF3.index == row].centroid.x.item()
            tempJSON = json.loads(blank_json)
            data = uk_lineage_json(code, DF)
            tempJSON['data']['values'] = data
            tempJSON['encoding']['longitude']['datum'] = long
            tempJSON['encoding']['latitude']['datum'] = lat
            multimapPIE['layer'].append(tempJSON)
        multimapPIE['layer'][0]['data']['values'] = region_json['features']
        multimapPIE['layer'][1]['data']['values'] = HBSet_json['features']
        #VegaLite(multimapPIE)
        return multimapPIE


def getSampleData_final(MetadataDF, HBTranslation, HBCode_translation):
    cog_meta = MetadataDF
    cog_meta['central_sample_id'] = cog_meta['sequence_name'].str.split('/', expand=True)[1]
    cog_meta = cog_meta.loc[cog_meta['adm1'].isin(['UK-SCT', 'UK-ENG', 'UK-WAL'])]
    cog_meta['HBName'] = cog_meta.loc[:, 'adm2'].map(HBTranslation)
    cog_meta['HBCode'] = cog_meta.loc[:, 'HBName'].map(HBCode_translation)
    cog_meta_clean = adm2cleaning(cog_meta)
    return cog_meta_clean


def central_surrounding_regions(HBCODE, neighborW, MAP):
    subregion = subMapExtractor(HBCODE, neighborW, MAP)
    central = subregion.loc[subregion['HBCode'] == HBCODE]
    neighbors = subregion.loc[subregion['HBCode'] != HBCODE]
    return central, neighbors, subregion


def adm2_to_centralHBCode(sampleframe, translation_dict, HbtoCode):
    HBs = []
    sampleframe['adm2'] = sampleframe['adm2'].str.upper()
    sampleframe = adm2cleaning(sampleframe, samplecsv=True)
    # print(sampleframe)
    adm2 = sampleframe['adm2'].to_list()
    # print(adm2)
    for each in adm2:
        if each in translation_dict:
            HBs.append(translation_dict[each])
    if len(HBs) > 0:
        centralHB = max(HBs, key=HBs.count)
        centralHBCode = HbtoCode[centralHB]
        # [k for k, v in a.items() if v == b]
        # print(centralHBCode)
        return centralHBCode
    elif len(HBs) == 0:
        return None


def defineDateRestriction(samplesDF, windowSize):
    datedSamples = samplesDF.dropna(subset=['Collection_Date'])
    if len(datedSamples) == 0:
        print('No collection dates specified, will revert to using all available data for local lineage analysis.')
        return None
    if len(datedSamples) > 0:
        datedSamples['Collection_Date'] = pd.to_datetime(datedSamples['Collection_Date'], format="%Y-%m-%d")
        datedSamples.sort_values(by=['Collection_Date'], inplace=True)
        start = datedSamples.iloc[[0]]['Collection_Date'].item() + pd.DateOffset(days=-windowSize)
        end = datedSamples.iloc[[-1]]['Collection_Date'].item() + pd.DateOffset(days=windowSize)
        daterange = {'start_date': start, 'end_date': end}
        return daterange


def finaliseMapping(boardMAP):
    ###Get contiguity neighbors for mainland
    MAP_W = Queen.from_dataframe(boardMAP, idVariable='HBCode')
    ##Hacky fix to link Fife and Lothian
    final_W = getForthBridge(boardMAP, MAP_W)
    return final_W


def hbcode_hbname_translation(mapDF):
    ## Create translation dict for HB codes and names
    translation = dict(zip(mapDF.HBName, mapDF.HBCode))
    return translation


def do_date_restriction(cogDF, samplecsv, start, end, window=7, restriction_bool="False"):
    if restriction_bool == 'True':
        if start is not None and end is not None:
            date_list = {'start_date': date_start, 'end_date': date_end}
            cogOut = dateRestriction(cogDF, date_list)
            return cogOut
        else:
            date_list = defineDateRestriction(inputSamples, window)
            if date_list is not None:
                cogOut = dateRestriction(cogDF, date_list)
                return cogOut
            else:
                cogOut = cogDF
                return cogOut
    else:
        return cogDF


def tableget(subFrame, cogDF):
    return (tabulateLins(subFrame['HBCode'], cogDF, subFrame['HBName']))

# ~~~~~~~~~~  JSON definitions  ~~~~~~~~~~~~~

singlemappie_blank = '''{
  "$schema": "https://vega.github.io/schema/vega-lite/v4.13.1.json",
    "width": 1000,
    "height": 1500,
    "projection": {"type":"mercator"},
    "layer": [
      {
          "data": {
          "values":""
          },

        "mark": {"type":"geoshape", "fill": "#DCE1DE", "stroke": "white"}
    },
      {
          "data": {
          "values":""
          },

        "mark": {"type":"geoshape", "fill": "#9CC5A1", "stroke": "white"}
    }

    ],
    "resolve": {"scale":{"color":"independent", "radius":"independent", "scale":"independent", "theta":"independent"} }
    }
    '''

multimappie_blank = """{
"$schema": "https://vega.github.io/schema/vega-lite/v4.13.1.json",
"width": 1000,
"height": 1500,
"projection": {"type":"mercator"},
"layer": [
  {
      "data": {
      "values":""
      },

    "mark": {"type":"geoshape", "fill": "#DCE1DE", "stroke": "white"}
},
  {
      "data": {
      "values":""
      },

    "mark": {"type":"geoshape", "fill": "#BFD9C2", "stroke": "white"}
}
],
"resolve": {"scale":{"color":"shared", "radius":"shared", "scale":"independent", "theta":"independent"} }
}
"""

blank_json = """{
  "data": {
    "values": ""
  },

  "layer": [
    {
    "mark": {"type": "arc", "innerRadius": 15, "cornerRadius":0.5,"padAngle":0.05, "stroke": "white","strokeWidth":0.5, "thetaOffset":-1.5},
    "encoding": { 
    "text": {"field": "UK_Lineage", "type": "nominal"}
        }
    },
        {
    "mark": {"type": "text", "radiusOffset": 30, "radius":30, "fill":"black", "fontWeight":"bold", "fontSize":15, "font":"Helvetica Neue","angle":0,"align":"center", "thetaOffset":-1.5, "stroke":"white", "strokeWidth":0.01},
    "transform":[{"window": [{"op": "rank","as": "rank"}],
                  "sort": [{ "field": "Count", "order": "descending" }]
    }, 
    {"filter": "datum.rank <= 5"}
    ],
    "encoding": { 
      "text": {"field": "UK_Lineage", "type": "nominal"}
        }
    }
],
  "encoding": {
    "theta": {"field": "Count", "type": "quantitative", "sort":"descending", "stack": true},
    "radius": {"field": "Count", "type": "quantitative", "scale": {"type": "sqrt", "zero": true, "range": [10, 100]}},
    "color": {"field": "UK_Lineage", "sort":"", "type": "ordinal", "scale":{"range": ["#51A3A3", "#095D45", "#F58300", "#313B72", "#5B8E7D", "#504D56", "#605D65", "#6E6C73", "#7B7980","#87858C"]}},
      "latitude" : {"datum":""},
      "longitude":{"datum":""}
        }
}"""

# ~~~~~~~~~~~  Code execution ~~~~~~~~~~~~~~~
# ~~~~~~~~

HBTranslation=pickle.load(open(translator, 'rb'))
COGDATA=pd.read_csv(argsIN.cog_metadata)
inputSamples = pd.read_csv(argsIN.user_sample_data)
mainland_boards=gp.read_file(mapfile)
if len(date_pair) == 2:
    date_start=date_pair[0]
    date_start=pd.to_datetime(date_start, format="%Y-%m-%d")
    date_end=date_pair[1]
    date_end=pd.to_datetime(date_end, format="%Y-%m-%d")
if len(date_pair) == 1:
    date_start=date_pair[0]
    date_start=pd.to_datetime(date_start, format="%Y-%m-%d")
    date_end=pd.to_datetime('today')
if len(date_pair) == 0:
    date_start=None
    date_end=None

# ~~~~~~~
mainland_W=finaliseMapping(mainland_boards)
HBCode_name_translation=hbcode_hbname_translation(mainland_boards)

# ~~~~~~~
mainland_W=finaliseMapping(mainland_boards)
HBCode_name_translation=hbcode_hbname_translation(mainland_boards)

# ~~~~~~~
#mainland_boards=update_adm15(mainland_boards)
###Checking user defined dates###
if argsIN.date_restriction != 'True':
    cog_meta=do_date_restriction(COGDATA, inputSamples, date_start, date_end)
if argsIN.date_restriction == 'True':
    if date_window == '':
        cog_meta=do_date_restriction(COGDATA, inputSamples, date_start, date_end, restriction_bool='True')
    if date_window != '':
        cog_meta=do_date_restriction(COGDATA, inputSamples, date_start, date_end, restriction_bool='True', window=int(date_window))
### Final restricted meta
cog_final = getSampleData_final(cog_meta, HBTranslation, HBCode_name_translation)
## Proessing input csv ##
Central_HB_code=adm2_to_centralHBCode(inputSamples, HBTranslation, HBCode_name_translation)

# ~~~~~~~
if Central_HB_code is not None:
    ## Get the localised regions ##
    central, neighboring, submap = central_surrounding_regions(Central_HB_code, mainland_W, mainland_boards)
    ## Generate tabular data for each region ##
    for row, frame in central.iterrows():
        HB_name, MDTable = tableget(frame, cog_final)
        with open(os.path.join(outDIR, f'{HB_name}_central_lineageTable.md'), 'w') as f:
            f.write(f'### {HB_name}\n')
            f.write(f'{MDTable}\n\n')
    for row, frame in neighboring.iterrows():
        HB_name, MDTable = tableget(frame, cog_final)
        with open(os.path.join(outDIR, f'{HB_name}_neighboring_lineageTable.md'), 'w') as f:
            f.write(f'### {HB_name}\n')
            f.write(f'{MDTable}\n\n')

    for each in [central, neighboring]:
        for row, frame in each.iterrows():
            HB_name, MDTable = tableget(frame, cog_final)

            with open(os.path.join(outDIR, f'{HB_name}_lineageTable.md'), 'w') as f:
                f.write(f'### {HB_name}\n')
                f.write(f'{MDTable}\n\n')
# ~~~~~~~
####  TARGET FILE WOULD COME FROM HERE; 'central_map_ukLink.vl.json'
if Central_HB_code is not None:
    ## Get the localised regions ##
    central, neighboring, submap = central_surrounding_regions(Central_HB_code, mainland_W, mainland_boards)
    ## Generated Mapping ##
    centralmapOUT=mapProduce(central, cog_final, submap)
    neighboringmapOUT=mapProduce(neighboring, cog_final, submap)
    regionmapOUT=mapProduce(submap, cog_final, submap, Central_HB_code)
    for mapjson, location in zip([centralmapOUT,neighboringmapOUT,regionmapOUT],['central', 'neighboring', 'region']):
        with open(os.path.join(argsIN.output_temp_dir, f'{location}_map_ukLin.vl.json'), 'w') as f:
            json.dump(mapjson, f)
            ##### Shell commands required
            # subprocess.call(['npx', 'vl2png', f'{location}_map_ukLin.vl.json', f{location}_map_ukLin.png])