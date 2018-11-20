import osmnx as ox, networkx as nx, matplotlib.cm as cm, pandas as pd, numpy as np
import geopandas as gpd
import functools
import community
import math
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.sparse import linalg
import pysal as ps

from shapely.geometry import Point, LineString, Polygon, MultiPolygon, mapping, MultiLineString
from math import sqrt
import pandas as pd
from shapely.ops import cascaded_union, linemerge
pd.set_option('precision', 10)
from time import sleep
import sys
import utilities as uf
"""
This set of functions is designed for extracting the computational image of the city,
Nodes, paths and districts are extracted with street network analysis, employing the primal and the dual graph representations.
Landmarks are extracted via a salience assessment process.

While the use of the terms "nodes" and "edges" can be cause confusion between the graph component and the Lynch components, nodes and edges are here used instead of vertexes and links to be consistent with NetworkX definitions.

"""
        
	

def select_buildings(city_buildings, area_to_clip, height_field, base = None, area_obstructions = None):
    
    buildings = city_buildings[city_buildings.geometry.within(area_to_clip.geometry.loc[0])]

    city_buildings["area"] = city_buildings['geometry'].area
    city_buildings["height"] = city_buildings[height_field]
    if (base is None): city_buildings["base"] = 0
    else: city_buildings["base"] = city_buildings[base]
        
    city_buildings = city_buildings[city_buildings['area'] > 199]
    city_buildings = city_buildings[['height', 'base','geometry', 'area']]
    buildings['buildingID'] = buildings.index.values.astype(int)
    
    if  (area_obstructions is None): area_obstructions = area_to_clip.geometry.loc[0].buffer(800)
    else: area_obstructions = area_obstructions.geometry.loc[0]
    
    obstructions = city_buildings[city_buildings.geometry.within(area_obstructions)]
    buildings = city_buildings[city_buildings.geometry.within(area_to_clip.geometry.loc[0])]
    buildings['buildingID'] = buildings.index.values.astype(int)
    buildings['r_height'] = buildings['height'] + buildings['base']
    return(buildings, obstructions)

def get_obstructions(city_buildings, area_to_clip):
        obstructions = city_buildings[city_buildings.geometry.within(area_to_clip.geometry.loc[0])]
        return obstructions
    
def structural_properties(buildings_gdf, obstructions, street_gdf):
    
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


def visibility(buildings_gdf, sight_lines):
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


def land_use_from_polygons(buildings_gdf, other_gdf, column, land_use_field):
    
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


def local_salience(nodes_gdf, landmarks_gdf):
    
    spatial_index = landmarks_gdf.sindex
    nodes_gdf['loc_land'] = 'None'
    nodes_gdf['loc_scor'] = 'None'
    index_geometry = nodes_gdf.columns.get_loc("geometry")+1

    for row in nodes_gdf.itertuples():
       
        g = row[index_geometry] #geometry
        b = g.buffer(50)    
        
        possible_matches_index = list(spatial_index.intersection(b.bounds))
        possible_matches = landmarks_gdf.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(b)]
        
        if (len(precise_matches) == 0): continue
    
        precise_matches = precise_matches.sort_values(by = 'lScore_sc', ascending=False).reset_index()
        precise_matches = precise_matches.round({'lScore_sc':3})
        
        list_local = precise_matches['buildingID'].tolist()
        list_scores = precise_matches['lScore_sc'].tolist()
        nodes_gdf.set_value(row[0], 'loc_land', list_local)
        nodes_gdf.set_value(row[0], 'loc_scor', list_scores)
    
    return nodes_gdf


def distant_landmarks(nodes_gdf, landmarks_gdf, sight_lines, smaller_area = False, all_nodes = None):
    
    # keeping relevant landmarks and the sight lines that point at them
    
    if smaller_area == True:
        nodes_gdf['oldID'] = nodes_gdf['nodeID']
        nodes_gdf['coordinates'] = nodes_gdf[['x', 'y']].apply(tuple, axis=1)
        all_nodes['coordinates'] = all_nodes[['x', 'y']].apply(tuple, axis=1)
        nodes_merged = pd.merge(nodes_gdf, all_nodes[['nodeID','coordinates']], left_on = 'coordinates', right_on = 'coordinates', 
                                how = 'left')
        nodes_merged['nodeID'] = nodes_merged['nodeID_y'] 
        nodes_merged.drop(['nodeID_x', 'nodeID_y', 'coordinates'], axis = 1, inplace = True)
        node_gdf = nodes_merged
    
    relevant = landmarks_gdf[landmarks_gdf.gScore_sc > 0.5]
    relevant = relevant.round({'gScore_sc':3})
    index_relevant = landmarks_gdf['buildingID'].values.astype(int) 
    sightLines_to_relevant = sight_lines[sight_lines['buildingID'].isin(index_relevant)]

    # per each node, the sight lines to the relevant landmarks are extracted. 
    # The visibility score of each node is the sum of the score of the visible landmarks visible from it, regardless the direction

    nodes_gdf['dist_land'] = None
    nodes_gdf['dist_scor'] = None
    nodes_gdf['anchors'] = None

    index_nodeID = nodes_gdf.columns.get_loc("nodeID")+1
    index_geometry = nodes_gdf.columns.get_loc("geometry")+1
    
    sindex = relevant.sindex
    for row in nodes_gdf.itertuples():
        sight_node = sightLines_to_relevant[sightLines_to_relevant['nodeID'] == row[index_nodeID]] 
        index_relevant_fromNode = list(sight_node['buildingID'].values.astype(int))
        relevant_fromNode = relevant[relevant['buildingID'].isin(index_relevant_fromNode)] 
        relevant_fromNode.sort_values(by = 'gScore_sc', ascending = False, inplace = True)
        if len(relevant_fromNode) > 0:
            nodes_gdf.set_value(row[0], 'dist_land', relevant_fromNode['buildingID'].tolist())  
            nodes_gdf.set_value(row[0], 'dist_scores',  relevant_fromNode['gScore_sc'].tolist())  
        
        # anchors around the node (intended as a destination)
        b = row[index_geometry].buffer(500)
        possible_matches_index = list(sindex.intersection(b.bounds))
        possible_matches = relevant.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(b)]
        if len(precise_matches) > 0:
            precise_matches.sort_values(by = 'gScore_sc', ascending = False, inplace = True)
            anchors = precise_matches.iloc[0:5] # first five rows of dataframe
            nodes_gdf.set_value(row[0], 'anchors',  anchors['buildingID'].tolist())  
        else: continue
#         nodes_gdf.set_value(row[0], 'anchors_scores',  list(anchors['gScore'].values.astype(float)))    
    
    
    if smaller_area == True:
        nodes_gdf['nodeID'] = nodes_gdf['oldID']
        nodes_gdf.drop(['oldID'], axis = 1, inplace = True)
    
    return nodes_gdf 


def advance_visibility_buildings(landmarks_gdf, obstructions_gdf):
    
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


def advance_visibility_nodes(nodes_gdf, edges_gdf, landmarks_gdf, field_view = True):
    
    visibility_polygons = pd.DataFrame(columns = ['nodeID', 'geometry'])
    if field_view == True: visibility_polygons = pd.DataFrame(columns = ['nodeID', 'from_Node', 'geometry'])
    
    degree_field_view = 140
    index_x = nodes_gdf.columns.get_loc("x")+1 
    index_y = nodes_gdf.columns.get_loc("y")+1 
    index_nodeID = nodes_gdf.columns.get_loc("nodeID")+1 
    index_geometry = nodes_gdf.columns.get_loc("geometry")+1 
    
    index_v = edges_gdf.columns.get_loc("v")+1 
    index_u = edges_gdf.columns.get_loc("u")+1 
    
    to_sub = (360-degree_field_view)/2
    sindex = landmarks_gdf.sindex
    coming_from = 0
    counter = 0
    
    for row in nodes_gdf.itertuples():
        sys.stdout.write('\r')
        sys.stdout.write("progress: "+str(int(counter/len(nodes_gdf)*100))+" %")
        sys.stdout.flush()
        sleep(0.25)
        counter += 1
        
        origin = (float("{0:.10f}".format(row[index_x])), float("{0:.10f}".format(row[index_y])))
        
        possible_obstacles_index = list(sindex.intersection(row[index_geometry].buffer(220).bounds))
        possible_obstacles = landmarks_gdf.iloc[possible_obstacles_index]
        
        if field_view == True:
                
            segments = edges_gdf[(edges_gdf['v'] == row[index_nodeID]) | (edges_gdf['u'] == row[index_nodeID])]
            zeroCoord = uf.get_coord_angle(origin, distance = 200, angle = 0.5)
            
            for s in segments.itertuples():

                if s[index_v] == row[index_nodeID]:     # take ndode u
                    otherCoord = (nodes_gdf.x[nodes_gdf['nodeID'] == s[index_u]].tolist()[0], 
                           nodes_gdf.y[nodes_gdf['nodeID'] == s[index_u]].tolist()[0])
                    coming_from = s[index_u]

                else: #takes node v
                    otherCoord = (nodes_gdf.x[nodes_gdf['nodeID'] ==  s[index_v]].tolist()[0], 
                           nodes_gdf.y[nodes_gdf['nodeID'] == s[index_v]].tolist()[0])
                    coming_from = s[index_v]

                line0 = ((origin), (zeroCoord))
                lineB = ((origin), (otherCoord))

                if otherCoord[0] > origin[0]: # east

                    diff = uf.ang(line0, lineB)
                    start = diff+to_sub
                    end = start+degree_field_view
                    if end >360: end = end-360

                else: # west
                    diff = 180-uf.ang(line0, lineB)
                    start = 180+diff+to_sub
                    if start > 360: start = start-360

                    end = start+degree_field_view 
                    if end > 360: end = end-360

                first_lineCoord = uf.get_coord_angle(origin, distance = 200, angle = start)
                last_lineCoord = uf.get_coord_angle(origin, distance = 200, angle = end)

                first_line = LineString([row[index_geometry], Point(first_lineCoord)])
                last_line = LineString([row[index_geometry], Point(last_lineCoord)])


                i = start
                list_lines = []

                while (((end > start) & ((i <= end) & (i >= start))) |
                       ((end < start) & (((i >= start) & (i > end)) | ((i < start) & (i <= end))))):

                    coords = uf.get_coord_angle(origin, distance = 200, angle = i)
                    line = LineString([row[index_geometry], Point(coords)])

                    obstacles = possible_obstacles[possible_obstacles.crosses(line)]
                    ob = cascaded_union(obstacles.geometry)

                    if len(obstacles > 0):
                        t = line.intersection(ob)

                        try:
                            intersection = t[0].coords[0]
                        except:
                            intersection = t.coords[0]

                        lineNew = LineString([row[index_geometry], Point(intersection)])

                    else: lineNew = line

                    list_lines.append(lineNew)
                    i = i+1
                    if i > 360: i = i - 360


                list_points = [Point(origin)]
                for i in list_lines: list_points.append(Point(i.coords[1]))
                list_points.append(Point(origin))
                poly = Polygon([[p.x, p.y] for p in list_points])
                
                visibility_polygons.loc[-1] = [row[index_nodeID], coming_from, poly] 
                visibility_polygons.index = visibility_polygons.index + 1
            
        else:
            start = 0.5
            i = start
            list_lines = []
            
            while(i <= 360):

                coords = uf.get_coord_angle(origin, distance = 200, angle = i)
                line = LineString([row[index_geometry], Point(coords)])

                obstacles = possible_obstacles[possible_obstacles.crosses(line)]
                ob = cascaded_union(obstacles.geometry)

                if len(obstacles > 0):
                    t = line.intersection(ob)

                    try:
                        intersection = t[0].coords[0]
                    except:
                        intersection = t.coords[0]

                    lineNew = LineString([row[index_geometry], Point(intersection)])

                else: lineNew = line

                list_lines.append(lineNew)
                i = i+1

            list_points = [Point(origin)]
            for i in list_lines: list_points.append(Point(i.coords[1]))
            list_points.append(Point(origin))
            poly = Polygon([[p.x, p.y] for p in list_points])
            
            visibility_polygons.loc[-1] = [row[index_nodeID], poly] 
            visibility_polygons.index = visibility_polygons.index + 1
        
    
    visibility_polygons_gdf = gpd.GeoDataFrame(visibility_polygons.loc[:, visibility_polygons.columns != 'geometry'], 
                                               crs = nodes_gdf.crs, geometry = visibility_polygons['geometry'])
    return visibility_polygons_gdf
    
    

def visibility_matrix(landmarks_gdf, visibility_polygons_gdf, nodes_gdf):
    
    columns = ['buildingID'] + nodes_gdf['nodeID'].astype(str).tolist()
    visibility_matrix = pd.DataFrame(columns = columns)
    visibility_matrix['buildingID'] = landmarks_gdf['buildingID']

    index_buildingID = visibility_matrix.columns.get_loc("buildingID")+1  
    index_polygon = visibility_polygons_gdf.columns.get_loc("geometry")+1  
    index_nodeID = visibility_polygons_gdf.columns.get_loc("nodeID")+1  
    sindex = visibility_polygons_gdf.sindex

    for row in visibility_matrix.itertuples(): 

        l_polygon = landmarks_gdf['geometry'][landmarks_gdf['buildingID'] == row[index_buildingID]].tolist()[0]
        b = l_polygon.buffer(200)
        possible_matches_index = list(sindex.intersection(b.bounds))
        possible_matches = visibility_polygons_gdf.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(b)]

        for p in precise_matches.itertuples():
            v_polygon = p[index_polygon]
            column = str(p[index_nodeID])

            if ((l_polygon.intersects(v_polygon)) | (l_polygon.touches(v_polygon))):
                visibility_matrix.set_value(row[0], column, True)

    visibility_matrix.fillna(value = False, inplace=True)
    visibility_matrix = visibility_matrix.astype(int)
    
    return(visibility_matrix)
    
