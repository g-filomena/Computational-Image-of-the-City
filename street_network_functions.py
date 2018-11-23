import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import osmnx as ox, networkx as nx, matplotlib.cm as cm, pandas as pd, numpy as np, geopandas as gpd
import functools
import community
import math
from math import sqrt
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.sparse import linalg
import pysal as ps

from shapely.geometry import Point, LineString, Polygon, MultiPolygon, mapping, MultiLineString
from shapely.ops import cascaded_union, linemerge, nearest_points
pd.set_option('precision', 10)

import utilities as uf

"""
This set of functions is designed for extracting the computational image of the city,
Nodes, paths and districts are extracted with street network analysis, employing the primal and the dual graph representations.
Landmarks are extracted via a salience assessment process.

While the use of the terms "nodes" and "edges" can be cause confusion between the graph component and the Lynch components, nodes and edges are here used instead of vertexes and links to be consistent with NetworkX definitions.

"""

	
## Graph preparation functions ###############
	
def get_fromOSM(type_download, place, network_type, epsg, distance = 7000): 
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
    
    if type_download == 'shapefilePolygon':
        file = gpd.read_file(place)
        polygon = file.geometry.loc[0]
        G = ox.graph_from_polygon(polygon, network_type = network_type, simplify = True)

    elif type_download == 'OSMpolygon':
        query = ox.osm_polygon_download(place, limit=1, polygon_geojson=1)
        OSMplace = query[0]['display_name']
        G = ox.graph_from_place(OSMplace, network_type = network_type, simplify = True)
        
    elif type_download == 'distance_from_address':
        G = ox.graph_from_address(place, network_type = network_type, distance = distance, simplify = True)

    else: # (type_download == 'OSMplace')
        G = ox.graph_from_place(place, network_type = network_type, simplify = True)
    
    if project: G = ox.project_graph(G)
    
    for i, item in enumerate(G.edges()):
        if isinstance(G[item[0]][item[1]][0]['osmid'], (list,)):
            G[item[0]][item[1]][0]['osmid'] = G[item[0]][item[1]][0]['osmid'][0]
            
    # using OSMNx to download data from OpenStreetMap      
    nodes = ox.graph_to_gdfs(G, nodes=True, edges=False, node_geometry=True, fill_edge_geometry=False)
    nodes_gdf = nodes.drop(nodes[['highway', 'ref']], axis=1)
    edges = ox.graph_to_gdfs(G, nodes=False, edges=True, node_geometry=False, fill_edge_geometry=True)
    edges_gdf = edges[['geometry', 'length', 'osmid', 'u','v', 'highway','key', 'oneway', 'maxspeed','name']]
    
    # getting rid of OSMid and preparing geodataframes
    edges_gdf = edges_gdf.rename(columns = {'u':'old_u', 'v':'old_v'})
    nodes_gdf = nodes_gdf.reset_index(drop=True)
    nodes_gdf['old_nodeID'] = nodes_gdf.osmid.astype('int64')
    nodes_gdf['nodeID'] = nodes_gdf.index.values.astype(int)
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['old_nodeID', 'nodeID']], how='left', left_on="old_u", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'u'})
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['old_nodeID', 'nodeID']], how='left', left_on="old_v", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'v'})
    
    # resetting index                          
    edges_gdf = edges_gdf.reset_index(drop=True)
    edges_gdf['streetID'] = edges_gdf.index.values.astype(int)
                                            
    # columns to keep (u and v represent "from" and "to" node)
    nodes_gdf = nodes_gdf[['nodeID','x','y','geometry']]
    edges_gdf = edges_gdf[['streetID','u','v','key','geometry', 'length', 'highway','oneway', 'maxspeed','name']]
    
    # resolving lists 
    edges_gdf['highway'] = [x[0] if type(x) == list else x for x in edges_gdf['highway']]
    edges_gdf['name'] = [x[0] if type(x) == list else x for x in edges_gdf['name']]
    
    # finalising geodataframes
    nodes_gdf, edges_gdf = nodes_gdf.to_crs(epsg = epsg), edges_gdf.to_crs(epsg = epsg)
    nodes_gdf['x'], nodes_gdf['y'] = list(zip(*[(r.coords[0][0], r.coords[0][1]) for r in nodes_gdf.geometry]))
    return(nodes_gdf, edges_gdf)

def get_fromSHP(directory, epsg, crs, area = None, roadType_field = None, direction_field = None, speed_field = None, name_field = None):
    
    """
    The function loads a vector lines shapefile from a specified directory, along with the epsg coordinate code.
    It creates two geopandas dataframe, one for street junctions (vertexes) and one for street segments (links).
    
    The file (e.g. street network shapefile) is supposed to be already simplified (e.g. not to have pseudo-nodes).
    The geopandas dataframes are built assuming a planar undirected graph. 
     
    Parameters
    ----------
    directory: string
    epsg: int
    roadType_field: indicates the column name where the street type is stored
    direction_field: indicates the column name where information about direction (one-way, two-way) is stored
    speed_field: string, indicates the column name where the speed limit is stored
    name_field: string, indicates the column name where the street name is stored
    area: int
    
    Returns
    -------
    GeoDataFrames
    """
    
    # try reading street network from directory
    streets_gdf = gpd.read_file(directory).to_crs(epsg=epsg)
        
    # using a buffer to clip the area of study
    if (area != None):
        cn = streets_gdf.geometry.unary_union.centroid
        buffer = cn.buffer(area) 
        streets_gdf = streets_gdf[streets_gdf.geometry.within(buffer)]
        
    columns = [roadType_field, direction_field, speed_field, name_field]
    new_columns = ['highway','oneway', 'maxspeed','name']
    streets['from'] = "NaN"
    streets['to'] = "NaN"
    
    # creating the dataframes
    for n, i in enumerate(columns): 
        if (i is not None): streets_gdf[new_columns[n]] = streets_gdf[i]
     
    standard_columns = ['geometry', 'from', 'to']
    streets_gdf = streets_gdf[standard_columns + [new_columns[n] for n, i in enumerate(columns) if i is not None]]
    index_geo = streets_gdf.columns.get_loc("geometry")+1
    
    for row in streets_gdf.itertuples():
        new_line = []
        
        # removing Z coordinates (assuming a planar graph)
        line_coords = list(row[index_geo].coords)
        from_node = coord[0][0:2]
        to_node = coord[-1][0:2]
                
        for i in range(0, len(line_coord)):
            point = line_coords[i][0:2]
            new_line.append(point)

        geo_line = LineString([coor for coor in new_line])
        streets_gdf.set_value(row[0],'geometry', geo_line)
        streets_gdf.set_value(row[0], 'from', from_node)
        streets_gdf.set_value(row[0], 'to', to_node)
        
    streets_gdf = streets_gdf.loc[streets_gdf['from'] != streets_gdf['to']]
    unique_nodes_tmp = list(streets_gdf['to'].unique()) + list(streets_gdf['from'].unique())
    unique_nodes = list(set(unique_nodes_tmp))
    
    # assigning indexes
    streets_gdf.reset_index(inplace=True, drop=True)
    streets_gdf['streetID'] = streets_gdf.index.values.astype(int) 
    
    #preparing nodes geodataframe
    nodes_data = pd.DataFrame.from_records(unique_nodes, columns=['x', 'y']).astype('float')
    geometry = [Point(xy) for xy in zip(nodes_data.x, nodes_data.y)]
    nodes_gdf = gpd.GeoDataFrame(nodes_data, crs=crs, geometry=geometry)
    nodes_gdf.reset_index(drop=True, inplace = True)
    nodes_gdf['nodeID'] = nodes_gdf.index.values.astype(int)
    nodes_gdf['coordinates'] = list(zip(nodes_gdf.x, nodes_gdf.y))

    edges_tmp = pd.merge(streets_gdf, nodes_gdf[['nodeID','coordinates']], how='left', left_on="from", right_on="coordinates")
    edges_tmp.drop(edges_tmp[['coordinates']], axis = 1, inplace = True)
    edges_tmp.rename(columns = {'nodeID':'u'}, inplace = True)
    
    edges_gdf = pd.merge(edges_tmp, nodes_gdf[['nodeID','coordinates']], how='left', left_on="to", right_on="coordinates")
    edges_gdf = edges_gdf.drop(edges_gdf[['coordinates', 'from', 'to']], axis = 1)
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'v'})
    edges_gdf['length'] = gpd.GeoSeries(edges_gdf['geometry'].length) # computing length
    nodes_gdf.drop(['coordinates'], axis = 1, inplace = True)
        
    return(nodes_gdf, edges_gdf)
	

def reset_index_gdf(nodes_gdf, edges_gdf):
    """
    The function simply reset the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    
    edges_gdf = edges_gdf.rename(columns = {'u':'old_u', 'v':'old_v'})
    nodes_gdf = nodes_gdf.reset_index(drop = False)
    nodes_gdf['old_nodeID'] = nodes_gdf['index'].astype('int64')
    nodes_gdf['nodeID'] = nodes_gdf.index.values.astype('int64')
    
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['old_nodeID', 'nodeID']], how='left', left_on="old_u", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'u'})
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['old_nodeID', 'nodeID']], how='left', left_on="old_v", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'v'})

    edges_gdf.drop(['old_u', 'old_nodeID_x', 'old_nodeID_y', 'old_v'], axis = 1, inplace = True)
    nodes_gdf.drop(['old_nodeID', 'index'], axis = 1, inplace = True)
    edges_gdf = edges_gdf.reset_index(drop=True)
    edges_gdf['streetID'] = edges_gdf.index.values.astype(int)
    
    return(nodes_gdf, edges_gdf)

## Cleaning functions ###############

def double_nodes(nodes_gdf, edges_gdf):
    """
    The function simply resets the indexes of the two dataframes without losing the relative links
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    # the index of nodes_gdf has to be nodeID
    if nodes_gdf.index.values != nodes_gdf.nodeID.values: nodes_gdf.index =  nodes_gdf.nodeID
    nodes_gdf, edges_gdf =  nodes_gdf.copy(), edges_gdf.copy()
    
    print('Eliminating duplicate geometries - nodes..')
    sindex = nodes_gdf.sindex
    
    # detecting duplicate geometries
    G = nodes_gdf["geometry"].apply(lambda geom: geom.wkb)
    new_nodes = nodes_gdf.loc[G.drop_duplicates().index]
    to_edit = list(set(nodes_gdf.index.values.tolist()) - set((new_nodes.index.values.tolist())))
    if len(to_edit) == 0: return(nodes_gdf, edges_gdf) 
    else:
        
        # readjusting edges' nodes too, accordingly
        sindex = new_nodes.sindex
        for node in to_edit:

            geo = nodes_gdf.loc[node].geometry
            possible_matches_index = list(sindex.intersection(geo.bounds))
            possible_matches = new_nodes.iloc[possible_matches_index]
            precise_matches  = possible_matches[possible_matches.intersects(geo)]
            index = precise_matches.iloc[0].name
            
            # assigning the unique index to edges
            edges_gdf.loc[edges_gdf.u == node,'u'] = index
            edges_gdf.loc[edges_gdf.v == node,'v'] = index
        
    return(new_nodes, edges_gdf)
    

def fix_dead_ends(nodes_gdf, edges_gdf):
    """
    The function removes dead-ends. Eliminate nodes shared by only one segment, and the relative segment.
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    nodes_gdf =  nodes_gdf.copy()
    edges_gdf = edges_gdf.copy()
    
    print('Removing dead ends..')
    dd_u = dict(edges_gdf['u'].value_counts())
    dd_v = dict(edges_gdf['v'].value_counts())
    dd = {k: dd_u.get(k, 0) + dd_v.get(k, 0) for k in set(dd_u) | set(dd_v)}
    to_delete = {k: v for k, v in dd.items() if v == 1}
    if len(to_delete) == 0: return(nodes_gdf, edges_gdf)
    
    # removing edges and nodes
    to_delete_list = list(to_delete.keys())
    nodes_gdf.drop(to_delete_list, axis = 0 , inplace = True)
    edges_gdf = edges_gdf[~edges_gdf['u'].isin(to_delete_list)]
    edges_gdf = edges_gdf[~edges_gdf['v'].isin(to_delete_list)]

    return(nodes_gdf, edges_gdf)

def nodes_simplified(edges_gdf):
    """
    The function checks the presence of pseudo-junctions, by using the edges_gdf geodataframe.
     
    Parameters
    ----------
    edges_gdf: GeoDataFrame, street segments
   
    Returns
    -------
    boolean
    """
    
    simplified = True
    dd_u = dict(edges_gdf['u'].value_counts())
    dd_v = dict(edges_gdf['v'].value_counts())
    dd = {k: dd_u.get(k, 0) + dd_v.get(k, 0) for k in set(dd_u) | set(dd_v)}
    to_edit = {k: v for k, v in dd.items() if v == 2}
    if len(to_edit) == 0: return(simplified)
    simplified = False
            
    return(simplified)

def edges_simplified(edges_gdf):
    """
    The function checks the presence of possible duplicate geometries in the edges_gdf geodataframe.
     
    Parameters
    ----------
    edges_gdf: GeoDataFrame, street segments
   
    Returns
    -------
    boolean
    """
    
    simplified = True 
    edges_gdf['code'] = None
    edges_gdf['code'][edges_gdf['v'] >= edges_gdf['u']] = edges_gdf.u.astype(str)+"-"+edges_gdf.v.astype(str)
    edges_gdf['code'][edges_gdf['v'] < edges_gdf['u']] = edges_gdf.v.astype(str)+"-"+edges_gdf.u.astype(str)
    dd = dict(edges_gdf['code'].value_counts())
    dd = {k: v for k, v in dd.items() if v > 1}
    if len(dd) > 0: 
        print("Potential duplicate edges: ", len(dd))
        simplified = False

    return(simplified)

def simplify_graph(nodes_gdf, edges_gdf):
    """
    The function identify pseudo-nodes, namely nodes that represent intersection between only 2 segments.
    The segments are merged and the node is removed from the nodes_gdf geodataframe.
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    
    print('Removing pseudo-nodes..')
    nodes_gdf, edges_gdf = nodes_gdf.copy(), edges_gdf.copy()
    edges_gdf = edges_gdf[edges_gdf.highway != 'primary_link']
    
    # keeping only one item per node and counting its "appearances"
    dd_u, dd_v = dict(edges_gdf['u'].value_counts()), dict(edges_gdf['v'].value_counts())
    dd = {k: dd_u.get(k, 0) + dd_v.get(k, 0) for k in set(dd_u) | set(dd_v)}
    
    # editing the ones which onlu connect two edges
    to_edit = {k: v for k, v in dd.items() if v == 2}
    if len(to_edit) == 0: return(nodes_gdf, edges_gdf)
    to_edit_list = list(to_edit.keys())
    for nodeID in to_edit_list:
        tmp = edges_gdf[(edges_gdf['u'] == nodeID) | (edges_gdf['v'] == nodeID)].copy()           
        index_first = tmp.iloc[0].name # first segment index
        index_second = tmp.iloc[1].name # second segment index
        
        # Identifying the relationship between the two segments.
        # New node_u and node_v are assigned accordingly. A list of ordered coordinates is obtained for 
        # merging the geometries. 4 conditions:
        if (tmp.iloc[0]['u'] == tmp.iloc[1]['u']):  
            edges_gdf.set_value(index_first,'u', edges_gdf.loc[index_first]['v'])
            edges_gdf.set_value(index_first,'v', edges_gdf.loc[index_second]['v'])
            line_coordsA = list(tmp.iloc[0]['geometry'].coords)
            line_coordsB = list(tmp.iloc[1]['geometry'].coords)    
            line_coordsA.reverse()
        
        elif (tmp.iloc[0]['u'] == tmp.iloc[1]['v']): 
            edges_gdf.set_value(index_first,'u', edges_gdf.loc[index_second]['u'])
            line_coordsB = list(tmp.iloc[0]['geometry'].coords)
            line_coordsA = list(tmp.iloc[1]['geometry'].coords)                
        
        elif (tmp.iloc[0]['v'] == tmp.iloc[1]['u']): 
            edges_gdf.set_value(index_first,'v', edges_gdf.loc[index_second]['v'])
            line_coordsA = list(tmp.iloc[0]['geometry'].coords)
            line_coordsB = list(tmp.iloc[1]['geometry'].coords)  
        
        else: # (tmp.iloc[0]['v'] == tmp.iloc[1]['v']) 
            edges_gdf.set_value(index_first,'v', edges_gdf.loc[index_second]['u'])
            line_coordsA = list(tmp.iloc[0]['geometry'].coords)
            line_coordsB = list(tmp.iloc[1]['geometry'].coords)    
            line_coordsB.reverse()
            
        # checking that none edges with node_u == node_v have been created, if yes: drop them
        if edges_gdf.loc[index_first].u == edges_gdf.loc[index_first].v: 
            edges_gdf.drop([index_first, index_second], axis = 0, inplace = True)
            nodes_gdf.drop(nodeID, axis = 0, inplace = True)
            continue
        
        # obtaining coordinates-list in coherent order and merging
        new_line = line_coordsA + line_coordsB
        merged_line = LineString([coor for coor in new_line]) 
        edges_gdf.set_value(index_first, 'geometry', merged_line)
        
        # dropping the second segment, as the new geometry was assigned to the first edge
        edges_gdf.drop(index_second, axis = 0, inplace = True)
        nodes_gdf.drop(nodeID, axis = 0, inplace = True)
    
    return(nodes_gdf, edges_gdf)

def clean_network(nodes_gdf, edges_gdf, dead_ends = False):
    """
    It calls a series of functions (see above) to clean and remove dubplicate geometries or possible parallel short edges.
    It moreover removes pseudo-nodes and, optionally, dead ends.
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
    dead_ends: boolean
   
    Returns
    -------
    GeoDataFrames
    """
    
    nodes_gdf, edges_gdf =  nodes_gdf.copy(), edges_gdf.copy()
    nodes_gdf, edges_gdf = double_nodes(nodes_gdf, edges_gdf)

    index_key = edges_gdf.columns.get_loc("key")+1
    index_type = edges_gdf.columns.get_loc("highway")+1
    index_u, index_v = edges_gdf.columns.get_loc("u")+1, edges_gdf.columns.get_loc("v")+1
    index_geo = edges_gdf.columns.get_loc("geometry")+1
    edges_gdf = edges_gdf[edges_gdf['u'] != edges_gdf['v']] #eliminate node-lines
    edges_gdf.sort_index(inplace = True)  
    
    print('Cleaning and simplyfing network:')
    edges_gdf['code'], edges_gdf['coords'] = None, None
    cycle = 0
    
    while ((not edges_simplified(edges_gdf)) | (not nodes_simplified(edges_gdf))):

        processed = []
        edges_gdf['length'] = edges_gdf['geometry'].length # recomputing length, to account for small changes
        print('Cycle nr. ', cycle, '=============')
        cycle += 1
        
        # Assigning codes based on the edge's nodes. 
        # The string is formulated putting the node with lower ID first, regardless it being 'u' or 'v'
        edges_gdf['code'][edges_gdf['v'] >= edges_gdf['u']] = edges_gdf.u.astype(str)+"-"+edges_gdf.v.astype(str)
        edges_gdf['code'][edges_gdf['v'] < edges_gdf['u']] = edges_gdf.v.astype(str)+"-"+edges_gdf.u.astype(str)
        
        # Reordering coordinates to allow for comparison between edges
        edges_gdf['coords'] = [list(c.coords) for c in edges_gdf.geometry]
        edges_gdf['coords'][(edges_gdf.u.astype(str)+"-"+edges_gdf.v.astype(str)) != edges_gdf.code] = [
            list(x.coords)[::-1] for x in edges_gdf.geometry]
        
        # dropping duplicate-geometries Edges
        print('Eliminating duplicate geometries - edges..')
        G = edges_gdf['geometry'].apply(lambda geom: geom.wkb)
        edges_gdf = edges_gdf.loc[G.drop_duplicates().index]
        
        # dropping edges with same geometry but with coords in different orders (depending on their directions)    
        edges_gdf['tmp'] = edges_gdf['coords'].apply(tuple, 1)
        edges_gdf.drop_duplicates(['tmp'], keep = 'first', inplace = True)
        
        print('Checking lines with same u-v combination of nodes..')
        dd = dict(edges_gdf['code'].value_counts())
        dd = {k: v for k, v in dd.items() if v > 1} # keeping u-v combinations that appear more than once
        
        # iterate through possible double edges for each specific combination of possible duplicates
        for key, value in dd.items():
            tmp = edges_gdf[edges_gdf.code == key].copy()
            
            # sorting the temporary gdf by length, the shortest is then used as a term of comparison
            tmp.sort_values(['length'], ascending = True, inplace = True)
            u, v, geo_line, index_line = tmp.iloc[0]['u'], tmp.iloc[0]['v'], tmp.iloc[0]['geometry'], tmp.iloc[0].name 
            
            # iterate through all the other edges with same u-v nodes                                
            for rowC in tmp.itertuples():
                if rowC[0] == index_line: continue
                uC, vC, geo_lineC, index_lineC = rowC[index_u], rowC[index_v], rowC[index_geo], rowC[0] 
                
                # if this edge is 20% longer than the edge identified in the outer loop, delete it
                if geo_lineC.length > (geo_line.length * 1.20): 
                    edges_gdf.drop(index_lineC, axis = 0, inplace = True) 
                    continue
                
                # else draw a center-line, replace the geometry of the outer-loop segment with the CL, drop the segment of the inner-loop
                else:
                    cl = center_line(u, v, uC, vC, geo_line, geo_lineC)
                    edges_gdf.set_value(index_line,'geometry', cl)
                    edges_gdf.drop(index_lineC, axis = 0, inplace = True)

        # dead-ends and                                  
        if dead_ends: nodes_gdf, edges_gdf = fix_dead_ends(nodes_gdf, edges_gdf)
        nodes_gdf, edges_gdf = simplify_graph(nodes_gdf, edges_gdf)  
    
    # keeps nodes which are actually used by the edges in the geodataframe
    to_keep = list(set(list(edges_gdf['u'].unique()) + list(streets_gdf['v'].unique())))
    nodes_gdf = nodes_gdf[nodes_gdf['nodeID'].isin(to_keep)]
    
    edges_gdf.drop(['code', 'coords', 'tmp'], axis = 1, inplace = True)
    nodes_gdf['nodeID'] = nodes_gdf.nodeID.astype(int)
    print("Done all =========")  
    
    return(nodes_gdf, edges_gdf)

def correct_edges(nodes_gdf, edges_gdf):
    """
    The function checks that the edges coordinates are consistent with their relative u and v nodes'coordinates.
    It can be necessary to run the function after having cleaned the network
     
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments
   
    Returns
    -------
    GeoDataFrames
    """
    
    print("Correcting edges coordinates..")
    
    index_u, index_v  = edges_gdf.columns.get_loc("u")+1, edges_gdf.columns.get_loc("v")+1 
    index_geo = edges_gdf.columns.get_loc("geometry")+1 
    
    for row in edges_gdf.itertuples():
        u = nodes_gdf.loc[row[index_u]]['nodeID']
        v = nodes_gdf.loc[row[index_v]]['nodeID']
        line_coords = list(row[index_geo].coords)
        line_coords[0] = (nodes_gdf.loc[u]['x'], nodes_gdf.loc[u]['y'])
        line_coords[-1] = (nodes_gdf.loc[v]['x'], nodes_gdf.loc[v]['y'])
        geo_line = (LineString([coor for coor in line_coords]))
        edges_gdf.set_value(row[0],"geometry", geo_line)

    print("Done")
                                            
    return(nodes_gdf, edges_gdf)


## Obtaining graph ###############

def graph_fromGDF(nodes_gdf, edges_gdf, nodes_attributes, edges_costs):
    """
    It creates from two geopandas dataframes (street junctions and street segments) a NetworkX graph, passing by a OSMnx function.
    In the lists 'nodes_attributes' and 'edges_costs' please specify attributes that you want to preserve and store in the graph
    representation.
    
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments  
    nodes_attributes, edges_costs: lists
    
    Returns
    -------
    GeoDataFrames
    """
    
    # making use of OSMNx
    nodes_gdf.gdf_name, edges_gdf.gdf_name  = 'Nodes_list', 'Edges_list' 
    edges_gdf['key'] = 0
    G = ox.gdfs_to_graph(nodes_gdf, edges_gdf)
    
    # remapping the network in NetworkX
    # reassigning coordinates
    t = G.nodes()
    pos = {}
    for i, item in enumerate(t): pos[i] = (t[item]['x'],t[item]['y'], t[item]['nodeID'])
    Ng = nx.Graph() #Empty graph
    Ng = Ng.to_undirected()
    Ng.add_nodes_from(pos.keys()) #Add nodes preserving coordinates

    for i, item in enumerate(Ng.nodes()):
        Ng.node[item]['x'] = pos[item][0]
        Ng.node[item]['y'] = pos[item][1]
        Ng.node[item]['nodeID'] = pos[item][2]
        for attribute in nodes_attributes: Ng.node[item][attribute] = nodes[attribute][nodes['nodeID'] == pos[item][2]].tolist()[0]
                                                                                       
    for i, item in enumerate(G.edges()):
        Ng.add_edge(item[0], item[1])
        Ng[item[0]][item[1]]['streetID'] = G[item[0]][item[1]][0]['streetID']
        for cost in edges_costs: Ng[item[0]][item[1]][cost] = G[item[0]][item[1]][0][cost]
        
    return(Ng)

## Building geo-dataframes for dual graph representation ###############

def dual_gdf(nodes_gdf, edges_gdf, crs):
    """
    It creates two dataframes that are supposed to generate the dual graph of a street network. The nodes_dual gdf contains edges 
    centroids, the edges_dual gdf, instead, contains links between the street segment centroids. Those dual edges link real street segment 
    that share a junction. The centroids are stored with the original edge streetID, while the dual edges are associated with several
    attributes computed on the original street segments (distance between centroids, deflection angle).
    
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames, nodes and street segments  
    crs: dictionary
    
    Returns
    -------
    GeoDataFrames
    """
    if edges_gdf.index.values != edges_gdf.streetID.values: edges_gdf.index =  edges_gdf.streetID
    
    # computing centroids                                       
    centroids_gdf = edges_gdf.copy()
    centroids_gdf['centroid'] = centroids_gdf['geometry'].centroid
    centroids_gdf['intersecting'] = None
    
    index_u = centroids_gdf.columns.get_loc("u")+1
    index_v = centroids_gdf.columns.get_loc("v")+1
    index_streetID = centroids_gdf.columns.get_loc("streetID")+1
         
    # find_intersecting segments and storing them in the centroids gdf
    processed = []
    for c in centroids_gdf.itertuples():
        intersections = []
        from_node = c[index_u]
        to_node = c[index_v]
    
        possible_intersections = centroids_gdf.loc[(centroids_gdf['u'] == from_node) |
                        (centroids_gdf['u'] == to_node) |
                        (centroids_gdf['v'] == to_node) |
                        (centroids_gdf['v'] == from_node)]

        for p in possible_intersections.itertuples():
            if ((c[0]==p[0]) | ((c[0], p[0]) in processed) | ((p[0], c[0]) in processed)): continue
        
            else:
                intersections.append(p[index_streetID])  # appending streetID
                processed.append((p[0],c[0]))
    
        centroids_gdf.set_value(c[0],'intersecting', intersections)
        
    # creating vertexes representing street segments (centroids)
    centroids_data = centroids_gdf[['streetID', 'intersecting', 'length']]
    geometry = centroids_gdf['centroid']
    nodes_dual = gpd.GeoDataFrame(centroids_data, crs=crs, geometry=geometry)
    nodes_dual['x'] = [x.coords.xy[0][0] for x in centroids_gdf['centroid']]
    nodes_dual['y'] = [y.coords.xy[1][0] for y in centroids_gdf['centroid']]
    
    # creating fictious links between centroids
    edges_dual = pd.DataFrame(columns=['u','v', 'key', 'geometry', 'length'])

    index_length = nodes_dual.columns.get_loc("length")+1
    index_streetID_nd = nodes_dual.columns.get_loc("streetID")+1
    index_intersecting = nodes_dual.columns.get_loc("intersecting")+1
    index_geo = nodes_dual.columns.get_loc("geometry")+1
    
    # connecting nodes which represent street segments thare a linked in the actual street network                                        
    for row in nodes_dual.itertuples():
        
        streetID = row[index_streetID_nd] #streetID of the relative segment
        length = row[index_length]
                                            
        # intersecting segments:  # i is the streetID                                      
        for i in list(row[index_intersecting]):     
            length_i =  nodes_dual['length'][nodes_dual.streetID == i][i]
            distance = (length+length_i)/2
        
            # adding a row with u-v, key fixed as 0, Linestring geometry 
            # from the first centroid to the centroid intersecting segment 
            ls = LineString([row[index_geo], nodes_dual.loc[i]['geometry']])
            edges_dual.loc[-1] = [streetID, i, 0, ls, distance] 
            edges_dual.index = edges_dual.index + 1
            
    edges_dual = edges_dual.sort_index(axis=0)
    geometry = edges_dual['geometry']
    edges_dual = gpd.GeoDataFrame(edges_dual[['u','v', 'key', 'length']], crs=crs, geometry=geometry)

    index_lineA = edges_dual.columns.get_loc("u")+1
    index_lineB = edges_dual.columns.get_loc("v")+1
    
    for row in edges_dual.itertuples():

        # retrieveing original lines from/to
        geo_lineA = edges_gdf[edges_gdf.index == row[index_lineA]].geometry.iloc[0]
        geo_lineB = edges_gdf[edges_gdf.index == row[index_lineB]].geometry.iloc[0]
        
        # computing angles in degrees and radians
        deflection = uf.ang_geoline(geo_lineA, geo_lineB, degree = True)
        deflection_rad = uf.ang_geoline(geo_lineA, geo_lineB, degree = False)
           
        # setting values                                    
        edges_dual.set_value(row[0],'deg', deflection)
        edges_dual.set_value(row[0],'rad', deflection_rad)
        
    return(nodes_dual, edges_dual)

def get_dual_graph(nodes_dual, edges_dual, edge_costs):
    """
    The function generates a NetworkX graph from dual-nodes and -edges geopandas dataframes.
            
    Parameters
    ----------
    nodes_dual, edges_dual: GeoDataFrames
    edge_costs: list of edges attributes/costs

    Returns
    -------
    GeoDataFrames
    """
   
    nodes_dual.gdf_name = 'Dual_list'
    Gr = ox.gdfs_to_graph(nodes_dual, edges_dual)
    n = Gr.nodes()
    pos = {}

    for l, item in enumerate(n): pos[l] = (n[l]['x'], n[l]['y'], n[l]['streetID'])
        
    DG = nx.Graph() #Empty graph
    DG = DG.to_undirected()
    DG.add_nodes_from(pos.keys()) # Add nodes preserving coordinates
    
    for i, item in enumerate(DG.nodes()):
        DG.node[i]['x']=pos[i][0]
        DG.node[i]['y']=pos[i][1]
        DG.node[i]['streetID']=pos[i][2]
        
    for i, item in enumerate(Gr.edges()):
        DG.add_edge(item[0], item[1])
        for cost in edge_costs: DG[item[0]][item[1]][cost] = Gr[item[0]][item[1]][0][cost]
        
    return(DG)

def dual_id_dict(dict_values, G, nodeAttribute):
    """
    It could be used when one deals with a dual graph and wants to reconnect some analysis conducted on this representation to the analysis
    conducted on the primal graph. For instance, it takes the dictionary containing the betweennes-centrality values of the nodes in the
    dual graph, and associates these features to the corresponding edgeID (nodes in dual graph represent real edges).
    
    Parameters
    ----------
    dict_values: dictionary, of nodeID and centrality values (or other computation)
    G: networkx multigraph
    nodeAttribute: string, attribute of the node to link
    
    Returns
    -------
    dictionary
    """
    
    view = dict_values.items()
    ed_list = list(view)
    ed_dict = {}
    for p in ed_list: ed_dict[graph.node[p[0]][nodeAttribute]] = p[1] #Attribute and measure
        
    return(ed_dict)

	
## Natural roads extraction: the function has to be called twice for each segment, to check in both directions ###############
    """
    The concet of natural road regards the cognitive perception/representation of road entities, regardless changes in names or
    interruptions. Rather, different street segments are "mentally merged" according to continuity rules (here based on the deflection
    angle and the egoistic choice process).
    """
                                            
def natural_roads(streetID, naturalID, direction, nodes_gdf, edges_gdf): 
    """
    This function takes a direction "to" or "fr" and two geopandas dataframes, one for roads one for nodes (or junctions).
    Since this function works only for one segment at the time, in one direction, it has to be run in a for loop where iterating through
    all the edges, noth directions. See function below.
    
    Parameters
    ----------
    streetID: int, next road to examine
    naturalID: int, current naturalID 
    direction: string, {'to', 'fr'}
    edges_gdf: geopandas dataframe
    nodes_gdf: geopandas dataframe
    """

    # initialise variable
    angles = {}
    directions_dict = {}
    index_geo = edges_gdf.columns.get_loc("geometry")+1
    index_u, index_v = edges_gdf.columns.get_loc("u")+1, edges_gdf.columns.get_loc("v")+1
    index_nID = edges_gdf.columns.get_loc("naturalID")+1

    to_node, from_node = edges_gdf.loc[streetID]['u'], edges_gdf.loc[streetID]['v']
    geo = edges_gdf.loc[streetID]['geometry']
      
    # continue from the to_node or from the from_node    
    if (direction == "to"): intersecting = edges_gdf[(edges_gdf['u'] == to_node) | (edges_gdf['v'] == to_node)]
    else: intersecting = edges_gdf[(edges_gdf['u'] == from_node) | (edges_gdf['v'] == from_node)]
    if (len(intersecting) == 0): return
    
    # check all possible deflection angles with the intersecting roads identified
    for row_F in intersecting.itertuples():
        if ((streetID == row_F[0]) | (row_F[index_nID] > 0)): continue
        to_node_F = row_F[index_u]
        from_node_F = row_F[index_v]
        geo_F = row_F[index_geo]

        # where to go next, in case?
        if (to_node == to_node_F): towards = "fr"
        elif (to_node == from_node_F): towards = "to"
        elif (from_node == from_node_F): towards = "to"
        else: towards = "fr"

        # measuring deflection angle, adding it to the dictionary, if lower than 45 degrees
        deflection = uf.ang_geoline(geo, geo_F, degree = True)
        if (deflection >= 45): continue
        else:
            angles[row_F[0]] = deflection # dictionary with streetID and angle
            directions_dict[row_F[0]] = towards # dictionary with streetID and direction
    
    # No natural continuations
    if (len(angles) == 0):
        edges_gdf.set_value(row[0], 'naturalID', naturalID)
        return
   
    # selecting the best continuation and continuing in its direction
    else:
        angles_sorted = sorted(angles, key = angles.get) 
                                            
        # taking the streetID of the segment which form the gentlest angle with the segment examined
        matchID = angles_sorted[0] 
        edges_gdf.set_value(row[0], 'naturalID', naturalID)
        natural_roads(matchID, natural_id, directions_dict[matchID], nodes_gdf, edges_gdf)                    
                                                                                    
def run_natural_roads(nodes_gdf, edges_gdf): 
    """
    Run the natural_roads function on an entire geodataframe of street segments.
    The geodataframes are supposed to be cleaned and can be obtained via the functions "get_fromOSM(place)" or "get_fromSHP(directory, 
    epsg)". Please clean the graph before running.
    
    Parameters
    ----------
    nodes_gdf, edges_gdf: GeoDataFrames
    
    Returns
    -------
    GeoDataFrames
    """
    # only before starting the loop - the function is also called within the loop
    
    if (not nodes_simplified(edges_gdf)) | (not edges_simplified(edges_gdf)): raise GraphError('The graph is not simplified!')
    edges_gdf['naturalID'] = 0
    edges_gdf.index = edges_gdf.streetID
    nodes_gdf.index = nodes_gdf.nodeID
    naturalID = 1  
    
    for row in edges_gdf.itertuples():
        if (row[index_nID] > 0): continue # if already assigned to a natural road
        snf.natural_roads(row[0], naturalID, "fr", nodes_gdf, edges_gdf) 
        snf.natural_roads(row[0], naturalID, "to", nodes_gdf, edges_gdf) 
        naturalID = naturalID + 1
                                            
    return(nodes_gdf, streets_gdf)
    
## Centrality functions ###############

def straightness_centrality(G, weight, normalized = True):
    """
    Straightness centrality compares the length of the path between two nodes with the straight line that links them capturing a 
    centrality that refers to ‘being more directly reachable’. (Porta, S., Crucitti, P. & Latora, V., 2006b. The Network Analysis Of Urban
    Streets: A Primal Approach. Environment and Planning B: Planning and Design, 33(5), pp.705–725.)
    
    Parameters
    ----------
    G: networkx multigraph
    weight: string, edges weight
    normalized: boolean
    
    Returns
    -------
    dictionary
    """
    path_length = functools.partial(nx.single_source_dijkstra_path_length, weight = weight)
    nodes = G.nodes()
    straightness_centrality = {}

    # Initialize dictionary containing all the node id and coordinates
    # coord_nodes = get_nodes_coords(Node, Session)
    coord_nodes = nodes_dict(G)

    for n in nodes:
        straightness = 0
        sp = path_length(G,n)

        if len(sp) > 0 and len(G) > 1:
            # start computing the sum of euclidean distances
            for target in sp:
                if n != target and target in coord_nodes:
                    network_dist = sp[target]
                    euclidean_dist = euclidean_distance(*coord_nodes[n]+coord_nodes[target])
                    straightness = straightness + (euclidean_dist/network_dist)

            straightness_centrality[n] = straightness
            if normalized: straightness_centrality[n] = straightness * (1.0/(len(G)-1.0) )
        else:
            straightness_centrality[n]=0.0

    return straightness_centrality

def weight_nodes(nodes_gdf,  points_gdf, G, buffer, name):
    
    """
    Given a nodes and a services/points geodataframes, the function assign an attribute to nodes in the graph G (prevously derived from 
    nodes_gdf) based indeed on the amount of features in poinst_gdf in a buffer around each node. 
    
    Parameters
    ----------
    nodes_gdf: geodataframe
    points_gdf: geodataframe
    G: networkx multigraph
    buffer: float, distance around the node within looking for services or point features
    name: string, attribute name
    
    Returns
    -------
    networkx multidigraph
    """
    nodes_gdf[name] = None
    index_geo = nodes_gdf.columns.get_loc("geometry")+1
    sindex = points_gdf.sindex
    for row in nodes_gdf.itertuples():

        g = row[index_geo]
        fil = g.buffer(buffer)
        possible_matches_index = list(sindex.intersection(fil.bounds))
        possible_matches = points_gdf.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(fil.bounds)]
        weight = len(precise_matches)
        nodes_gdf.set_value(row[0], name, weight)

    for n in G.nodes():
        G.node[n][name] = nodes_gdf[name].loc[n]
    
    return G

def reach_centrality(G, weight, radius, attribute):
    """
    The measure contemplates the assignment of attributes (e.g. number of activities, population, employees in an area) to nodes and
    accounts for opportunities that are reachable along the actual street network as perceived by pedestrians’. The reach centrality of a
    node j, indicates the number of other nodes reachable from i, at the shortest path distance of r, where nodes are rewarded with a
    score (indicated by 'attribute') which indicates their importance. The function is readapted from: Sevtsuk, A. & Mekonnen, M., 2012.
    Urban Network Analysis: A New Toolbox For ArcGIS. Revue internationale de géomatique, 2, pp.287–305.

    Parameters
    ----------
    G: networkx multigraph
    weight: string, edges weight
    radius: float, distance from node within looking for other reachable nodes
    attribute: string, node attribute used to compute reach centralily. It indicates the importance of the node (e.g. number of services in     50mt buffer)
    
    Returns
    -------
    dictionary
    """
    
    path_length = functools.partial(nx.single_source_dijkstra_path_length, weight = weight)

    nodes = G.nodes()
    reach_centrality = {}
    coord_nodes = nodes_dict(G)

    for n in nodes:
        reach = 0
        sp = path_length(G, n)
        sp_radium = dict((k, v) for k, v in sp.items() if v <= radius)
        
        if len(sp_radium) > 0 and len(G) > 1:
            
            for target in sp_radium:
                if n != target and target in coord_nodes:
                    weight_target = G.node[target][attribute]
                    reach = reach + weight_target
            reach_centrality[n] = reach
        else:               
            reach_centrality[n]=0.0

    return reach_centrality

def local_betweenness_centrality(G, weight, radius):
    """
    The function computes betweenness centrality at the local level.
    ## Still to be optimised. 
      
    Parameters
    ----------
    G: networkx multigraph
    weight: string, edges weight
    radius: float, distance from node, within local betweenness is computed  
    
    Returns
    -------
    dictionary
    """
           
    path_length = functools.partial(nx.single_source_dijkstra_path_length, weight = weight)

    nodes = G.nodes()
    cb = {}
    coord_nodes = nodes_dict(G)

    for obj in nodes:
        sp = path_length(G,obj)
        sp_radium = dict((k, v) for k, v in sp.items() if v <= radius)

        to_keep = list(sp_radium.keys())
        G_small = nx.Graph(G.subgraph(to_keep))
        
        be = nx.betweenness_centrality(G_small, k=None, weight = 'length', normalized=False)
        cb[obj] = be[obj]
     
    return cb
	


    
