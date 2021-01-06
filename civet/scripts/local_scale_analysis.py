#!/usr/bin/env python3
import os
import warnings
import pickle
import geopandas as gp
import pandas as pd
from libpysal.weights import Queen, attach_islands, DistanceBand, set_operations
import json
import csv
import numpy as np
from collections import Counter
#from vega import VegaLite
import argparse
warnings.filterwarnings("ignore")
print("running local scale analysis")
# ~~~~~~~ Define input variables ~~~~~~~~~~

parser = argparse.ArgumentParser(description='Parse barcode info and minimap paf file, create report.')
parser.add_argument("--date-restriction", default=False, action="store", type=bool, dest="date_restriction")
parser.add_argument("--date-pair-start", default="", action="store", dest="date_pair_start")
parser.add_argument("--date-pair-end", default="", action="store", dest="date_pair_end")
parser.add_argument("--cog-meta-global", action="store", type=str, dest="cog_metadata")
parser.add_argument("--user-sample-data", action="store", type=str, dest="user_sample_data")
parser.add_argument("--combined-metadata",action="store",type=str,dest="combined_metadata")
parser.add_argument("--input-name", action="store", type=str, dest="input_name")
parser.add_argument("--date-window", action="store",required=False, type=int, dest="date_window")
parser.add_argument("--output-base-dir", action="store", type=str, dest="output_base_dir")
parser.add_argument("--output-temp-dir", action="store", type=str, dest="output_temp_dir")
parser.add_argument("--hb-translation", action="store", type=str, dest="hb_translation")
parser.add_argument("--uk-map", action="store", type=str, dest="uk_map")
#parser.add_argument("--civet-cat", action="store", type=str, dest="civet_cat_dir")

argsIN=parser.parse_args()

currentDir=os.getcwd()
outDIR=argsIN.output_base_dir

translator=argsIN.hb_translation
mapfile=argsIN.uk_map

date_pair=[]
for each in [argsIN.date_pair_start, argsIN.date_pair_end]:
    if each != "False":
        date_pair.append(each)


# ~~~~~~~~ Defining functions ~~~~~~~~~~

def adm2cleaning(data_cog, samplecsv=False):
    data_cog2 = ""
    if not samplecsv:
        data_cog2 = data_cog.dropna(subset=['adm2']).copy()
    else:
        data_cog2 = data_cog.copy()
    
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^MOTHERWELL$', 'NORTH_LANARKSHIRE', regex=True).copy()
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^GREATER_LONDON$', 'GREATER_LONDON', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^BORDERS$', 'SCOTTISH _ORDERS', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^PAISLEY$', 'RENFREWSHIRE', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^SUTTON-IN-ASHFIELD$', 'NOTTINGHAMSHIRE', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^LONDON$', 'GREATER_LONDON', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^KILMARNOCK$', 'EAST_AYRSHIRE', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^PERTH$', 'PERTHSHIRE_AND_KINROSS', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^ISLE OF ANGLESEY$', 'ANGLESEY', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^SHETLAND$', 'SHETLAND_ISLANDS', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^GREATER LONDONDERRY$', 'LONDONDERRY', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^STAFFORD$', 'STAFFORDSHIRE', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^INVERNESS$', 'HIGHLAND', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^KINGS NORTON$', 'WEST_MIDLANDS', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^GLASGOW CITY$', 'GLASGOW', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^CITY OF LONDON$', 'GREATER_LONDON', regex=True)
    data_cog2['adm2'] = data_cog2['adm2'].str.replace('^RHONDDA CYNON TAF$', 'RHONDDA, CYNON, TAFF', regex=True)
    return data_cog2


def dateRestriction(DFin, dateDict):
    DFin['sample_date'] = pd.to_datetime(DFin['sample_date'], format="%Y-%m-%d")
    datemask = (DFin['sample_date'] > dateDict['start_date']) & (DFin['sample_date'] <= dateDict['end_date'])
    DFout = DFin.loc[datemask]
    return DFout


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
      #the first part is never true, it never finds the HBCode being the same
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


def uk_lineage_json(HBCODE, cog_mainland_daterestricted):
    lineage_uk = cog_mainland_daterestricted.loc[cog_mainland_daterestricted['HBCode'] == HBCODE][
                     'uk_lineage'].value_counts()[:10].to_frame()
    lineage_global = cog_mainland_daterestricted.loc[cog_mainland_daterestricted['HBCode'] == HBCODE][
                         'lineage'].value_counts()[:10].to_frame()
    lineage_uk = lineage_uk.reset_index()
    lineage_global = lineage_global.reset_index()
    lineage_uk.columns = ['UK_Lineage', 'Count']
    lin_list = lineage_uk['UK_Lineage'].to_list()
    regionName = HBname_code_translation[HBCODE]
    lineage_uk['region'] = regionName
    # lothian_lineage_uk['UK_Lineage']=lothian_lineage_uk['UK_Lineage'].str.replace('UK','')
    lineage_global.columns = ['UK_Lineage', 'Count']
    lineage_uk_json = lineage_uk.to_json(orient='records')
    lineage_global_json = lineage_global.to_json(orient='records')
    return lineage_uk_json, lin_list


def lineageRanking(DFin, HBSet, centralCode):
    DF = DFin.copy()
    linsOut = ''
    colsOut = ''
    sortedcounts = DF.value_counts(subset=['HBCode', 'uk_lineage']).to_frame().reset_index().sort_values(by=0,
                                                                                                                ascending=False)
    if centralCode is not None:
        centralloc = [centralCode]
        noncentrallocs = HBSet.query('HBCode not in @centralloc')
        region = HBSet['HBCode'].to_list()
        filteredregionalcounts = sortedcounts.query('HBCode in @region')
        filteredcentralcounts = sortedcounts.query('HBCode in @centralloc')
        filterednoncentralcounts = sortedcounts.query('HBCode in @noncentrallocs')
        topLinlist = filteredcentralcounts['uk_lineage'][:10].to_list()
        restLinlist = list(filterednoncentralcounts.query('uk_lineage not in @topLinlist')['uk_lineage'].unique())
        regionalList = topLinlist + restLinlist
        filler = [graphcolourFiller for i in range(len(regionalList) - 10)]
        linsOut = topLinlist + regionalList
        colsOut = graphcoloursBase + filler
    else:
        region = HBSet['HBCode'].to_list()
        filteredregionalcounts = sortedcounts.query('HBCode in @region')
        topLinlist = filteredregionalcounts['uk_lineage'][:10].to_list()
        restLinlist = list(filteredregionalcounts.query('uk_lineage not in @topLinlist')['uk_lineage'].unique())
        regionalList = topLinlist + restLinlist
        filler = [graphcolourFiller for i in range(len(regionalList) - 10)]
        linsOut = topLinlist + regionalList
        colsOut = graphcoloursBase + filler
    return linsOut, colsOut


def mapProduce(HBSet, DF, region, centralCode=None):
    ##base schema
    singlecanvasJSON = json.loads(multimappie_blank)
    multicanvasJSON = json.loads(multimappie_blank)
    HBSet_json = HBSet.to_json()
    HBSet_json = json.loads(HBSet_json)
    region_json = region.to_json()
    region_json = json.loads(region_json)
    if len(HBSet) == 1:
        # holds JSON-ified data
        datacollated = []
        subfigs = []
        centralLins = ""
        for row, frame in HBSet.iterrows():
            code = frame['HBCode']
            regionName = HBname_code_translation[code]
            # gather dataset and easily processed desc. stats from pandas-based functions
            data, lin_list = uk_lineage_json(code, DF)
            # Append this dataset - allowing for multiple regions
            datacollated.append(data)
            # gather encoding for chart positioning
            lat = HBSet.loc[HBSet.index == row].centroid.y.item()
            long = HBSet.loc[HBSet.index == row].centroid.x.item()
            # load in the subgraph schema
            subfig = json.loads(blank_subfig_json)
            subfig['encoding']['longitude']['datum'] = long
            subfig['encoding']['latitude']['datum'] = lat
            # give subfig appropriate filter to correctly pull from primary dataset
            subfig['layer'][0]['transform'][0]['filter'] = f"datum.region == '{regionName}'"
            subfig['layer'][1]['transform'][0]['filter'] = f"datum.region == '{regionName}'"
            subfigs.append(json.dumps(subfig))
            # Determine lineages of central interest
            centralLins = lin_list
        # load in dataset base spec
        dataJSON = json.loads(blank_data_json)
        for each in datacollated:
            # json obj status lost with change in scope
            currentList = dataJSON['mainSet']
            updatedList = currentList + json.loads(each)
            dataJSON['mainSet'] = updatedList
        templayer = json.loads("{}")
        # gathering sub-graphs
        for each in subfigs:
            templayer.update(json.loads(each))
        # changing encoding schema; assign colour to focal lineages & allow for greyed out lineages outside this scope
        encodingJSON = json.loads(blank_data_encoding)
        encodingJSON['color']['scale']['domain'] = centralLins
        encodingJSON['color']['scale']['range'] = graphcoloursBase[:10]
        # amending elements on base schema
        singlecanvasJSON['layer'].append(templayer)
        singlecanvasJSON['datasets'] = dataJSON
        singlecanvasJSON['encoding'] = encodingJSON
        singlecanvasJSON['layer'][0]['data']['values'] = region_json['features']
        singlecanvasJSON['layer'][1]['data']['values'] = HBSet_json['features']
        # VegaLite(singlemapPIE)
        return singlecanvasJSON



    elif len(HBSet) > 1:
        # holds JSON-ified data
        datacollated = []
        subfigs = []
        listLins = ""
        for row, frame in HBSet.iterrows():
            code = frame['HBCode']
            regionName = HBname_code_translation[code]
            # gather dataset and easily processed desc. stats from pandas-based functions
            data, lin_list = uk_lineage_json(code, DF)
            # Append this dataset - allowing for multiple regions
            datacollated.append(data)
            # gather encoding for chart positioning
            lat = HBSet.loc[HBSet.index == row].centroid.y.item()
            long = HBSet.loc[HBSet.index == row].centroid.x.item()
            # load in the subgraph schema
            subfig = json.loads(blank_subfig_json)
            subfig['encoding']['longitude']['datum'] = long
            subfig['encoding']['latitude']['datum'] = lat
            # give subfig appropriate filter to correctly pull from primary dataset
            subfig['layer'][0]['transform'][0]['filter'] = f"datum.region == '{regionName}'"
            subfig['layer'][1]['transform'][0]['filter'] = f"datum.region == '{regionName}'"
            subfigs.append(json.dumps(subfig))
        # load in dataset base spec
        dataJSON = json.loads(blank_data_json)
        for each in datacollated:
            # json obj status lost with change in scope
            currentList = dataJSON['mainSet']
            updatedList = currentList + json.loads(each)
            dataJSON['mainSet'] = updatedList
        # gathering sub-graphs
        for each in subfigs:
            multicanvasJSON['layer'].append(json.loads(each))
        # changing encoding schema; assign colour to focal lineages & allow for greyed out lineages outside this scope
        encodingJSON = json.loads(blank_data_encoding)
        linmap, colmap = lineageRanking(DF, HBSet, centralCode)
        encodingJSON['color']['scale']['domain'] = linmap
        encodingJSON['color']['scale']['range'] = colmap
        multicanvasJSON['datasets'] = dataJSON
        multicanvasJSON['encoding'] = encodingJSON
        multicanvasJSON['layer'][0]['data']['values'] = region_json['features']
        multicanvasJSON['layer'][1]['data']['values'] = HBSet_json['features']
        # VegaLite(singlemapPIE)
        return multicanvasJSON

def decide_HB(metadata_df, HB_translation):

  possible_HBs = []
  extra_translation = {}

  for adm2 in metadata_df["adm2"]:
    if not pd.isnull(adm2):
      if "|" in adm2:
        HB = decide_single_HB(adm2, HBTranslation)  
        try:
          HB_translation[adm2] = HB
        except:
          HB_translation = ""
      elif "RHONDDA" in adm2:  
        HB_translation[adm2] = "Cwm Taf Morgannwg University Health Board"    

  return HB_translation

def decide_single_HB(adm2, HB_translation):

  possible_HBs = []
  possible_adm2s = adm2.split("|")

  for item in possible_adm2s:
    if "RHONDDA" in item:
      possible_HBs.append("Cwm Taf Morgannwg University Health Board")
    else:
      try:
        possible_HBs.append(HB_translation[item])
      except:
        pass

  HB_counts = Counter(possible_HBs)
  try:
    HB = HB_counts.most_common(1)[0][0]
  except:
    HB = ""

  return HB


def getSampleData_final(MetadataDF, HBTranslation, HBCode_translation):
    cog_meta = MetadataDF
    cog_meta['central_sample_id'] = cog_meta['sequence_name'].str.split('/', expand=True)[1]
    cog_meta = cog_meta.loc[cog_meta['adm1'].isin(['UK-SCT', 'UK-ENG', 'UK-WLS', 'UK-NIR', 'Scotland', 'Northern_Ireland', "England", "Wales"])]
    HB_translation = decide_HB(cog_meta, HBTranslation)
    cog_meta['HBName'] = cog_meta['adm2'].map(HBTranslation)
    cog_meta['HBCode'] = cog_meta['HBName'].map(HBCode_translation)
    cog_meta_clean = adm2cleaning(cog_meta)H
    return cog_meta_clean


def central_surrounding_regions(HBCODE, neighborW, MAP):
    subregion = subMapExtractor(HBCODE, neighborW, MAP)
    central = subregion.loc[subregion['HBCode'] == HBCODE]
    neighbors = subregion.loc[subregion['HBCode'] != HBCODE]
    return central, neighbors, subregion


def adm2_to_centralHBCode(sampleframe, translation_dict, HbtoCode):
    HBs = []
    sampleframe['adm2'] = sampleframe['adm2'].str.upper().replace(" ","_")
    sampleframe = adm2cleaning(sampleframe, samplecsv=True)
    # print(sampleframe)
    adm2 = sampleframe['adm2'].to_list()
    # print(adm2)
    for each in adm2:
        if not pd.isnull(each):
          if each in translation_dict:
            HB = translation_dict[each]
          elif "|" in each:
            HB = decide_single_HB(each, translation_dict)
          elif "RHONDDA" in each:
            HB = "Cwm Taf Morgannwg University Health Board"
          HBs.append(HB)
    if len(HBs) > 0:
        centralHB = max(HBs, key=HBs.count)
        centralHBCode = HbtoCode[centralHB]
        # [k for k, v in a.items() if v == b]
        # print(centralHBCode)
        return centralHBCode
    elif len(HBs) == 0:
        return None


def supplement_sample_csv(sample_df,combined_metadata_df,input_name):

  col_list = sample_df.columns
  potential_date_cols = []
  for i in col_list:
    if "date" in i.lower():
        potential_date_cols.append(i)

  if "adm2" not in col_list:
    testing_similarity = combined_metadata_df['query'][combined_metadata_df["query"] == combined_metadata_df["closest"]]

    if len(testing_similarity) > 0 and not combined_metadata_df["adm2"].isnull().all():
      sample_df["adm2"] = combined_metadata_df["adm2"]
    else:
      return False

  if len(potential_date_cols) == 0:
    testing_similarity = combined_metadata_df['query'][combined_metadata_df["query"] == combined_metadata_df["closest"]]

    if len(testing_similarity) > 0:
      if "sample_date" in combined_metadata_df:
        sample_df["sample_date"] = combined_metadata_df["sample_date"]


  return sample_df

def defineDateRestriction(samplesDF, windowSize):
  
  col_list = samplesDF.columns
  potential_date_cols = []
  for i in col_list:
    if "date" in i.lower():
      potential_date_cols.append(i)

  if len(potential_date_cols) == 0:
    print('No dates found, will revert to using all available data for local lineage analysis.')
    return None

  collection = False
  for j in potential_date_cols:
    if "collection" in j.lower():
      date_col = j
      collection = True
      break
    else:
      date_col = j
      
  if collection:
    print(f"Collection date not found, using {j} to restrict local lineages by")

  datedSamples = samplesDF.dropna(subset=[date_col])
  if len(datedSamples) == 0:
      print('No collection dates specified, will revert to using all available data for local lineage analysis.')
      return None
  if len(datedSamples) > 0:
      datedSamples[date_col] = pd.to_datetime(datedSamples[date_col], format="%Y-%m-%d")
      datedSamples.sort_values(by=[date_col], inplace=True)
      start = datedSamples.iloc[[0]][date_col].item() + pd.DateOffset(days=-windowSize)
      end = datedSamples.iloc[[-1]][date_col].item() + pd.DateOffset(days=windowSize)
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

def do_date_restriction(cogDF, samplecsv, start, end, window=7, restriction_bool=False):
    if restriction_bool:
        if start and end:
            date_list = {'start_date': date_start, 'end_date': date_end}
            cogOut = dateRestriction(cogDF, date_list)
            return cogOut
        else:
            date_list = defineDateRestriction(samplecsv, window)
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



def hbname_hbcode_translation(mapDF):
    ## Create translation dict for HB codes and names
    translation = dict(zip(mapDF.HBCode, mapDF.HBName))
    return translation

# ~~~~~~~~~~  JSON definitions  ~~~~~~~~~~~~~

singlemappie_blank = '''
{
  "$schema": "https://vega.github.io/schema/vega-lite/v4.13.1.json",
  "width": 1000,
  "height": 1500,
  "datasets": {},
  "encoding": {},
  "projection": {"type": "mercator"},
  "layer": [
    {
      "data": {"values": ""},
      "mark": {"type": "geoshape", "fill": "#DCE1DE", "stroke": "white"}
    },
    {
      "data": {"values": ""},
      "mark": {"type": "geoshape", "fill": "#BFD9C2", "stroke": "white"}
    }
  ],
  "resolve": {
    "scale": {
      "color": "shared",
      "radius": "shared",
      "theta": "independent"
    }
  }
}
'''

multimappie_blank = """{
  "$schema": "https://vega.github.io/schema/vega-lite/v4.13.1.json",
  "width": 1000,
  "height": 1500,
  "datasets": {},
  "encoding": {},
  "projection": {"type": "mercator"},
  "layer": [
    {
      "data": {"values": ""},
      "mark": {"type": "geoshape", "fill": "#BFD9C2", "stroke": "white"}
    },
    {
      "data": {"values": ""},
      "mark": {"type": "geoshape", "fill": "#DCE1DE", "stroke": "white"}
    }
  ],
  "resolve": {
    "scale": {
      "color": "shared",
      "radius": "shared",
      "theta": "independent"
    }
  }
}
"""
##move datasets to outside layers
blank_data_json = """
{
    "mainSet": [],
    "data": {"name": "mainSet"},
    "encoding":{"sort": {"field": "Count", "order": "descending"}}
  }
"""

blank_data_encoding = """
{
      "color": {
        "field": "UK_Lineage",
        "type": "ordinal",
        "scale": {
          "domain": ["ordered UK_lin list here"],
          "range": ["colour range + filler greys here"]
        }
      }
}

"""

blank_subfig_json = """
{
      "layer": [
        {
          "mark": {
            "type": "arc",
            "innerRadius": 15,
            "cornerRadius": 0.5,
            "padAngle": 0.05,
            "stroke": "white",
            "strokeWidth": 0.5,
            "thetaOffset": -1.5
          },
          "transform": [{"filter": "datum.Region == 'region goes here'"}],
          "data": {"name": "mainSet"}
        },
        {
          "data": {"name": "mainSet"},
          "mark": {
            "type": "text",
            "radiusOffset": 30,
            "radius": 30,
            "fill": "black",
            "fontWeight": "bold",
            "fontSize": 15,
            "font": "Roboto",
            "angle": 0,
            "align": "center",
            "thetaOffset": -1.5,
            "stroke": "white",
            "strokeWidth": 0.01
          },
          "transform": [
            {"filter":"datum.region == 'region goes here'"},
            {
              "window": [{"op": "rank", "as": "rank"}],
              "sort": [{"field": "Count", "order": "descending"}]
            },
            {"filter": "datum.rank <= 5"}
          ],
          "encoding": {"text": {"field": "UK_Lineage", "type": "nominal"}}
        }
      ],
      "encoding": {
        "theta": {
          "field": "Count",
          "type": "quantitative",
          "sort": "descending",
          "stack": true
        },
        "radius": {
          "field": "Count",
          "type": "quantitative",
          "scale": {"type": "sqrt", "zero": true, "range": [10, 100]}
        },
        "latitude": {
          "datum": ""
        },
        "longitude": {
          "datum": ""
        }
      },
      "resolve": {"scale":{"theta":"independent"}}
}      
"""

graphcoloursBase = ["#51A3A3", "#536B78", "#FFA552", "#313B72", "#5B8E7D", "#ACCBE1", "#504D56", "#605D65", "#6E6C73",
                    "#7B7980"]
graphcolourFiller = "#87858C"

# ~~~~~~~~~~~  Code execution ~~~~~~~~~~~~~~~
# ~~~~~~~~

HBTranslation = pickle.load(open(translator, 'rb'))
COGDATA = pd.read_csv(argsIN.cog_metadata)
inputSamples = pd.read_csv(argsIN.user_sample_data)
combined_metadata = pd.read_csv(argsIN.combined_metadata)
mainland_boards=gp.read_file(mapfile)

#Replace spaces with underscores in the HB translations
new_translator = {}
for k,v in HBTranslation.items():
  new_translator[k.replace(" ","_")] = v
HBTranslation = new_translator


#if the inputs don't have adm2 or dates, but if some of them are in COG
inputSamples = supplement_sample_csv(inputSamples, combined_metadata, argsIN.input_name)
if type(inputSamples) == bool:
  print("NO ADM2 PRESENT, CANNOT PRODUCE LOCAL LINEAGE MAPPING")
else:

  #sort out the date inputs
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
  # ~~~~~~~~
  mainland_W=finaliseMapping(mainland_boards)
  HBCode_name_translation=hbcode_hbname_translation(mainland_boards)
  HBname_code_translation=hbname_hbcode_translation(mainland_boards)

  # ~~~~~~~
  #mainland_boards=update_adm15(mainland_boards)
  ###Checking user defined dates###
  if not argsIN.date_restriction:
      cog_meta=do_date_restriction(COGDATA, inputSamples, date_start, date_end)
  else:
      cog_meta=do_date_restriction(COGDATA, inputSamples, date_start, date_end, restriction_bool=True, window=int(argsIN.date_window))
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
          HB_code = frame["HBCode"]
          with open(os.path.join(outDIR, f'{HB_code}_central_lineageTable.md'), 'w') as f:
              f.write(f'### {HB_name}\n')
              f.write(f'{MDTable}\n\n')
      for row, frame in neighboring.iterrows():
          HB_name, MDTable = tableget(frame, cog_final)
          # HB_name = HB_name.replace(" ","_")
          HB_code = frame["HBCode"]
          with open(os.path.join(outDIR, f'{HB_code}_neighboring_lineageTable.md'), 'w') as f:
              f.write(f'### {HB_name}\n')
              f.write(f'{MDTable}\n\n')

  # ~~~~~~~
  if Central_HB_code is not None:
      ## Get the localised regions ##
      central, neighboring, submap = central_surrounding_regions(Central_HB_code, mainland_W, mainland_boards)
      ## Generated Mapping ##
      centralmapOUT=mapProduce(central, cog_final, submap)
      neighboringmapOUT=mapProduce(neighboring, cog_final, submap)
      regionmapOUT=mapProduce(submap, cog_final, submap, Central_HB_code)
      for mapjson, location in zip([centralmapOUT,neighboringmapOUT,regionmapOUT],['central', 'neighboring', 'region']):
          with open(os.path.join(argsIN.output_temp_dir, f'{location}_map_ukLin.vl.json'), 'w') as f:
              json.dump(mapjson, f, indent=True)
              ##### Shell commands required
              # subprocess.call(['npx', 'vl2png', f'{location}_map_ukLin.vl.json', f{location}_map_ukLin.png])