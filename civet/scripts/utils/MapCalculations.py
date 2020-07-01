import sjoin
import maup
import sjoin



def adm1_5_assignment(data_boundaries_URL, data_adm):
    data_NHSboard = gp.read_file(data_boundaries_URL)
    data_NHSboard = data_NHSboard.to_crs("EPSG:4326")
    data_adm = data_adm.to_crs("EPSG:4326")
    data_adm = data_adm.set_geometry(col=data_adm.geometry.centroid)
    adm1_5 = sjoin.sjoin_nearest(data_adm, data_NHSboard, how='left', search_radius=0.2)
    return adm1_5



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
