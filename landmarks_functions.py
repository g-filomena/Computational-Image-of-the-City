import pandas as pd, numpy as np, geopandas as gpd, matplotlib.pyplot as plt
import functools
import math
from math import sqrt
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, mapping, MultiLineString
from shapely.ops import cascaded_union, linemerge
from scipy.sparse import linalg
import pysal as ps

from time import sleep
import sys
pd.set_option('precision', 10)

import utilities as uf

"""
This set of functions is designed for extracting the computational image of the city,
Nodes, paths and districts are extracted with street network analysis, employing the primal and the dual graph representations.
Landmarks are extracted via a salience assessment process.

While the use of the terms "nodes" and "edges" can be cause confusion between the graph component and the Lynch components, nodes and edges are here used instead of vertexes and links to be consistent with NetworkX definitions.

"""
        
	

def select_buildings(city_buildings, area_to_clip, height_field, base_field = None, area_obstructions = None):
    """
    
    The function downloads and creates a simplified OSMNx graph for a selected area.
    Afterwards, geopandas dataframes for nodes and edges are created, assigning new nodeID and streeID identifiers.
    osmid are indeed heavy and confusing.
        
    Parameters
    type_download: string, {shapefilePolygon', 'OSMpolygon', 'distance_from_address', 'shapefilePolygon'}
    place: string, name of cities or areas in OSM
    network_type: string,  {'walk', 'bike', 'drive', 'drive_service', 'all', 'all_private', 'none'}
        what type of street or other network to get - from OSMNx paramaters
    epsg: int
    distance: float, only yse if type_download = 'distance from address'
    project: boolean
    ----------
    
    Returns
    -------
    GeoDataFrames
    """
    
    buildings = city_buildings[city_buildings.geometry.within(area_to_clip.geometry.loc[0])]

    city_buildings['area'] = city_buildings['geometry'].area
    city_buildings['height'] = city_buildings[height_field]
    if (base is None): city_buildings['base'] = 0
    else: city_buildings['base'] = city_buildings[base]
        
    city_buildings = city_buildings[city_buildings['area'] > 199]
    city_buildings = city_buildings[city_buildings['height'] >= 1]
    city_buildings = city_buildings[['height', 'base','geometry', 'area']]
    buildings['buildingID'] = buildings.index.values.astype(int)
    
    if  (area_obstructions is None): area_obstructions = area_to_clip.geometry.loc[0].buffer(800)
    else: area_obstructions = area_obstructions.geometry.loc[0]
    
    obstructions = city_buildings[city_buildings.geometry.within(area_obstructions)]
    buildings = city_buildings[city_buildings.geometry.within(area_to_clip.geometry.loc[0])]
    buildings['buildingID'] = buildings.index.values.astype(int)
    buildings['r_height'] = buildings['height'] + buildings['base']
    return(buildings, obstructions)
   
def structural_properties(buildings_gdf, obstructions, street_gdf):
        """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    
    index_geometry = buildings_gdf.columns.get_loc("geometry")+1 
    sindex = obstructions.sindex
    street_network = street_gdf.geometry.unary_union
    
    for row in buildings_gdf.itertuples():
        
        g = row[index_geometry]
        b200 = g.buffer(200)
        b150 = g.buffer(150)
        t = g.envelope
        coords = mapping(t)['coordinates'][0]
        d = [(Point(coords[0])).distance(Point(coords[1])), (Point(coords[1])).distance(Point(coords[2]))]

        width = min(d)
        length = max(d)

        buildings_gdf.set_value(row[0], 'width', width)
        buildings_gdf.set_value(row[0], 'length', length)
        
        possible_matches_index = list(sindex.intersection(b200.bounds))
        possible_matches = obstructions.iloc[possible_matches_index]   
        if(len(possible_matches_index) < 1): continue 
        
        possible_neigh_index = list(sindex.intersection(b150.bounds))
        possible_neigh = obstructions.iloc[possible_neigh_index]
        precise_neigh = possible_neigh[possible_neigh.intersects(b150)]
        buildings_gdf.set_value(row[0], 'neigh', len(precise_neigh)) #neighbours
        
        dist = g.distance(street_network)
        buildings_gdf.set_value(row[0], 'road', dist) #distance road

    buildings_gdf['fac'] = buildings_gdf['height']*(buildings_gdf.width) #facade area
    buildings_gdf.drop(['width','length'], axis=1, inplace = True)
    
    return buildings_gdf

def advance_visibility_buildings(landmarks_gdf, obstructions_gdf):
            """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    
    visibility_polygons = landmarks_gdf[['buildingID', 'geometry']].copy()
    landmarks_gdf['a_vis'] = 0.0
    
    sindex = obstructions_gdf.sindex
    index_geometry = landmarks_gdf.columns.get_loc("geometry")+1
    counter = 0
    
    for row in landmarks_gdf.itertuples():
        sys.stdout.write('\r')
        sys.stdout.write("progress: "+str(int(counter/len(landmarks_gdf)*100))+" %")
        sys.stdout.flush()
        sleep(0.25)
        
        counter += 1 
        origin = row[index_geometry].centroid
        exteriors = list(row[index_geometry].exterior.coords)
        no_holes = Polygon(exteriors)
        
        possible_obstacles_index = list(sindex.intersection(origin.buffer(2000).bounds))
        possible_obstacles = obstructions_gdf.iloc[possible_obstacles_index]
        possible_obstacles = obstructions_gdf[obstructions_gdf.geometry != row[index_geometry]]
        possible_obstacles = obstructions_gdf[~obstructions_gdf.geometry.within(no_holes)]

        start = 0.5
        i = start
        list_lines = []
            
        while(i <= 360):

            coords = uf.get_coord_angle([origin.x, origin.y], distance = 500, angle = i)

            line = LineString([origin, Point(coords)])
            obstacles = possible_obstacles[possible_obstacles.crosses(line)]
            ob = cascaded_union(obstacles.geometry)

            if len(obstacles > 0):
                t = line.intersection(ob)

                try:
                    intersection = t[0].coords[0]
                except:
                    intersection = t.coords[0]
                
                lineNew = LineString([origin, Point(intersection)])

            else: lineNew = line

            list_lines.append(lineNew)
            i = i+1

        list_points = [Point(origin)]

        
        for i in list_lines: list_points.append(Point(i.coords[1]))
        list_points.append(Point(origin))
        poly = Polygon([[p.x, p.y] for p in list_points])
        
        try:
            poly_vis = poly.difference(row[index_geometry])
        except:
            pp = poly.buffer(0)
            poly_vis = pp.difference(row[index_geometry])
        
        
        landmarks_gdf.set_value(row[0],'a_vis', poly_vis.area)
        
        try:
            if len(poly_vis) > 1: #MultiPolygon
                for i in range(0, len(poly_vis)): 
                    if (poly_vis[i].area < 100): del poly_vis[i]
        except:
            poly_vis = poly_vis

        visibility_polygons.set_value(row[0],'geometry', poly_vis)
        
    return landmarks_gdf, visibility_polygons



def visibility(buildings_gdf, sight_lines):
        """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    
    sight_lines = sight_lines.copy()
    sight_lines.drop(['Shape_Leng'], axis = 1, inplace = True)
    sight_lines['length'] = sight_lines['geometry'].length
    sight_lines = (sight_lines[['buildingID', 'length', 'nodeID']].groupby(['buildingID', 'nodeID'], 
                                            as_index = False)['length'].max())
    
    avg = (sight_lines[['buildingID', 'length']].groupby(['buildingID'], 
                                            as_index = False)['length'].mean())
    avg.rename(columns={'length':'mean_length'}, inplace=True)
    
    count = (sight_lines[['buildingID', 'length']].groupby(['buildingID'], 
                                            as_index=False)['length'].count())
    count.rename(columns={'length':'n_slines'}, inplace=True)

    tmp = sight_lines.set_index('buildingID')
    distant = tmp.groupby('buildingID').agg(lambda copy: copy.values[copy['length'].values.argmax()])
    distant = distant.reset_index()

    visibility_tmp = pd.merge(distant, avg, left_on= 'buildingID', right_on = 'buildingID')
    visibility = pd.merge(visibility_tmp, count, left_on= 'buildingID', right_on = 'buildingID')                    
    visibility.drop(['DIST_ALONG', 'Visibility', 'geometry'], axis = 1, errors = 'ignore', inplace=True)
    visibility.rename(columns = {'length':'dist','mean_length':'m_dist'}, inplace=True) 
    
    tmp = pd.merge(buildings_gdf, visibility[['buildingID', 'dist', 'm_dist', 'n_slines']], on = 'buildingID', how= 'left')  
    tmp['dist'].fillna((tmp['dist'].min()), inplace = True)
    tmp['m_dist'].fillna((tmp['m_dist'].min()), inplace=True)
    tmp['n_slines'].fillna((tmp['n_slines'].min()), inplace=True)
    
    col = ['dist', 'm_dist', 'n_slines']                     
    for i in col: uf.scaling_columnDF(tmp, i)
    tmp['vis'] = tmp['dist_sc']*0.5+tmp['m_dist_sc']*0.25+tmp['n_slines_sc']*0.25
    return tmp

def cultural_meaning(buildings_gdf, cultural_elements, score = None):
        """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    sindex = cultural_elements.sindex 
    buildings_gdf['cult'] = 0
    
    index_geometry = buildings_gdf.columns.get_loc("geometry")+1 
    
    for row in buildings_gdf.itertuples():
        g = row[index_geometry] #geometry
        possible_matches_index = list(sindex.intersection(g.bounds))
        possible_matches = cultural_elements.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(g)]
        
        if (score == None): cm = len(precise_matches)
        elif len(precise_matches) == 0: cm = 0
        else: cm = precise_matches[score].sum()
        
        buildings_gdf.set_value(row[0], 'cult', cm) #cultural meaning
     
    return buildings_gdf

def classify_lu(buildings_gdf, land_use):
    
    university = ['university', 'college', 'research']
    commercial = ['bank', 'service',  'commercial',  'retail', 'Retail',  'pharmacy', 'commercial;educa', 'shop', 'Commercial',
                  'supermarket', 'offices', 'foundation', 'office', 'books', 'Commercial services', 'Commercial Land', 
                  'Mixed Use Res/Comm',  'Commercial Condo Unit']
    
    residential = [ 'apartments', None, 'NaN', 'residential','flats', 'no', 'houses', 'garage', 'garages', 'building', 
                  'roof', 'storage_tank', 'shed', 'silo',  'parking',  'toilets', 'picnic_site','hut', 'information', 'viewpoint',
                  'atm',   'canopy', 'smokestack', 'greenhouse', 'fuel', 'Residential Condo Unit', 'Apartments 4-6 Units', 
                  'Residential Two Family', 'Apartments 7 Units above', 'Residential Single Family', 'Condominium Parking', 
                  'Residential Three Family', 'Condominium Master', 'Residential Land']
    
    attractions = ['Attractions', 'museum',  'castle', 'cathedral', 'attraction','aquarium', 'monument',  'gatehouse',
                   'terrace', 'tower', 'Attraction And Leisure']
    hospitality = [ 'hotel',  'hostel', 'guest_house']
    eating_drinking = [ 'restaurant', 'fast_food', 'cafe', 'bar',  'pub', 'Accommodation, eating and drinking',]
    public = ['post_office', 'townhall', 'public_building',  'library','civic', 'courthouse', 'public', 'embassy',
              'Public infrastructure', 'community_centre', 'parking', 'dormitory', 'Exempt', 'Exempt 121A']
    library = ['library']
    sport = ['stadium', 'Sport and entertainment', 'Sports Or Exercise Facility']
    entertainment = [ 'exhibition_centr','theatre', 'cinema']
    education = ['school', 'kindergarten', 'Education', 'Education and health']
    religious = ['church', 'place_of_worship','convent', 'rectory', 'Religious Buildings']
    emergency_service = [ 'fire_station','police', 'Emergency Service']
    transport = [ 'station', 'train_station']
    medical_care = ['hospital', 'doctors', 'dentist','clinic','veterinary', 'Medical Care']
    industrial = [ 'industrial', 'factory', 'construction', 'Manufacturing and production',  'gasometer', 'data_center']
    cultural = ['club_house','gallery', 'arts_centre','Cultural Facility']
    military = ['general aviation', 'Barracks']
    transport = ['Transport', 'Road Transport', 'station', 'subway_entrance', 'bus_station']

    buildings_gdf[land_use] = buildings_gdf[land_use].map( lambda x: 'university' if x in university
                                                              else 'commercial' if x in commercial
                                                              else 'residential' if x in residential
                                                              else 'attractions' if x in attractions
                                                              else 'library' if x in library
                                                              else 'hospitality' if x in hospitality
                                                              else 'eating_drinking' if x in eating_drinking
                                                              else 'public' if x in public
                                                              else 'sport' if x in sport
                                                              else 'entertainment' if x in entertainment
                                                              else 'education' if x in education
                                                              else 'religious' if x in religious
                                                              else 'emergency_service' if x in emergency_service
                                                              else 'industrial' if x in industrial
                                                              else 'cultural' if x in cultural
                                                              else 'transport' if x in transport
                                                              else 'medical_care' if x in medical_care
                                                              else 'military' if x in military
                                                              else x)
    
     buildings_gdf[buildings_gdf[land_use].str.contains('residential') | 
                   buildings_gdf[land_use].str.contains('Condominium') |
                  buildings_gdf[land_use].str.contains('Residential')] = 'residential'
    
     buildings_gdf[buildings_gdf[land_use].str.contains('commercial') | 
                   buildings_gdf[land_use].str.contains('Commercial')] = 'residential'
    
    return(obstructions_gdf)




def land_use_from_polygons(buildings_gdf, other_gdf, column, land_use_field):
        """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    
    buildings_gdf[column] = 'NaN'
    sindex = other_gdf.sindex

    index_geometry = buildings_gdf.columns.get_loc("geometry")+1 
    index_geometryPol = other_gdf.columns.get_loc("geometry")+1 

    for row in buildings_gdf.itertuples():

        g = row[index_geometry] #geometry

        possible_matches_index = list(sindex.intersection(g.bounds))
        possible_matches = other_gdf.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(g)]
        precise_matches['area'] = 0.0

        
        if (len(precise_matches) == 0): continue

        for row_other in precise_matches.itertuples():
            t = row_other[index_geometryPol]
            try:
                area_intersec = t.intersection(g).area
            except: 
                 continue
                    
            precise_matches.set_value(row_other[0], 'area', area_intersec)
    
        pm = precise_matches.sort_values(by='area', ascending=False).reset_index()
        
        if (pm['area'].loc[0] > (g.area * 0.59)):
            main_use = pm[land_use_field].loc[0]
            buildings_gdf.set_value(row[0], column, main_use)
            
        else: continue
        
    return buildings_gdf


def land_use_from_points(buildings_gdf, other_gdf, column, land_use_field):
        """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
        
    other_gdf['nr'] = 1
    buildings_gdf[column] = 'NaN'
    sindex = other_gdf.sindex

    index_geometry = buildings_gdf.columns.get_loc("geometry")+1 
    index_geometryPol = other_gdf.columns.get_loc("geometry")+1 

    for row in buildings_gdf.itertuples():

        g = row[index_geometry] #geometry

        possible_matches_index = list(sindex.intersection(g.bounds))
        possible_matches = other_gdf.iloc[possible_matches_index]

        if (len(possible_matches)==0): continue
        else:
            use = possible_matches.groupby([land_use_field],as_index=False)['nr'].sum().sort_values(by='nr',
                                    ascending=False).reset_index()
        
        main_use = use[land_use_field].loc[0]
        buildings_gdf.set_value(row[0], column, main_use)
               
    return buildings_gdf


        
def pragmatic_meaning(buildings_gdf):
        """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
        
    buildings_gdf['nr'] = 1
    sindex = buildings_gdf.sindex
    buildings_gdf['prag']= 0.0
    index_geometry = buildings_gdf.columns.get_loc("geometry")+1
    index_land_use = buildings_gdf.columns.get_loc("land_use")+1 

    for row in buildings_gdf.itertuples():
        g = row[index_geometry] #geometry
        b = g.buffer(200)
        use = row[index_land_use]

        possible_matches_index = list(sindex.intersection(b.bounds))
        possible_matches = buildings_gdf.iloc[possible_matches_index]
        precise_matches = possible_matches [possible_matches.intersects(b)]

        neigh = precise_matches.groupby(['land_use'], as_index=True)['nr'].sum()
        Nj = neigh.loc[use]
        #Pj = Nj/N

        Pj = 1-(Nj/precise_matches['nr'].sum())
        buildings_gdf.set_value(row[0], 'prag', Pj) #pragmatic meaning
        
    return buildings_gdf
        
        
def compute_scores(buildings_gdf):
        """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """

    col = ['area', 'fac', 'height', 'prag', 'a_vis','cult']
    for i in col: uf.scaling_columnDF(buildings_gdf, i)

    col = ['neigh', 'road']
    for i in col: uf.scaling_columnDF(buildings_gdf, i, inverse = True)     

    buildings_gdf['vScore'] = buildings_gdf['fac_sc']*0.30 + buildings_gdf['height_sc']*0.20 + buildings_gdf['vis']*0.5
    buildings_gdf['sScore'] = buildings_gdf['area_sc']*0.30 + buildings_gdf['neigh_sc']*0.20 + buildings_gdf['a_vis_sc']*0.30     +buildings_gdf['road_sc']*0.20

    col = ['vScore', 'sScore']
    for i in col: uf.scaling_columnDF(buildings_gdf, i)

    buildings_gdf['gScore'] = (buildings_gdf['vScore_sc']*0.50 + buildings_gdf['sScore_sc']*0.30 + buildings_gdf['cult_sc']*0.10 
                       + buildings_gdf['prag_sc']*0.10)

    uf.scaling_columnDF(buildings_gdf, 'gScore')
    
    return buildings_gdf

# Landmark Routing

def local_scores(landmarks_gdf, buffer_extension):
        """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    
    landmarks_gdf = landmarks_gdf.copy()
    spatial_index = landmarks_gdf.sindex
    index_geometry = landmarks_gdf.columns.get_loc("geometry")+1
    landmarks_gdf['lScore'] = 0.0
    
    col = ['area', 'fac', 'height', 'prag', 'a_vis','cult', 'vis']
    col_inverse = ['neigh', 'road']
   
    for row in landmarks_gdf.itertuples():
        b = row[index_geometry].centroid.buffer(buffer_extension)
        
#         LL = landmarks_gdf[landmarks_gdf.geometry.within(row[index_geometry].centroid.buffer(buffer_extension))]
        possible_matches_index = list(spatial_index.intersection(b.bounds))
        possible_matches = landmarks_gdf.iloc[possible_matches_index]
        LL = possible_matches[possible_matches.intersects(b)]
        
        for i in col: uf.scaling_columnDF(LL, i)
        for i in col_inverse: uf.scaling_columnDF(LL, i, inverse = True)
        
        LL['vScore'] =  LL['fac_sc']*0.30 + LL['height_sc']*0.20 + LL['vis']*0.5
        LL['sScore'] =  LL['area_sc']*0.30 + LL['neigh_sc']*0.20 + LL['a_vis_sc']*0.3 + LL['road_sc']*0.20
        LL['lScore'] =  LL['vScore_sc']*0.20 + LL['sScore_sc']*0.30 + LL['cult_sc']*0.10 + LL['prag_sc']*0.40
        
        localScore = float("{0:.3f}".format(LL['lScore'].loc[row[0]]))
        landmarks_gdf.set_value(row[0], 'lScore', localScore)
    uf.scaling_columnDF(landmarks_gdf, 'lScore', inverse = False)
        
    return landmarks_gdf



