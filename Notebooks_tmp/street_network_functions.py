import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
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
from shapely.ops import nearest_points
from math import sqrt
import pandas as pd
from shapely.ops import cascaded_union, linemerge
pd.set_option('precision', 10)

"""
This set of functions is designed for extracting the computational image of the city,
Nodes, paths and districts are extracted with street network analysis, employing the primal and the dual graph representations.
Landmarks are extracted via a salience assessment process.

While the use of the terms "nodes" and "edges" can be cause confusion between the graph component and the Lynch components, nodes and edges are here used instead of vertexes and links to be consistent with NetworkX definitions.

"""

# utilities
	
def plot_grad_GDF_manualBreaks(gdf, column, title):
    bins = [0.12, 0.25, 0.50, 0.75, 1.00]
    cl = ps.User_Defined(gdf[column], bins)
    
    f, ax = plt.subplots(1, figsize=(15, 15))
    gdf.assign(cl=cl.yb).plot(ax=ax, column='cl', categorical=True, k=5, cmap='OrRd', linewidth=0.5,
                              legend=True)
    f.suptitle(title)
    plt.axis('equal')
    
    leg = ax.get_legend()
    leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
    leg.get_texts()[0].set_text('0.00 - 0.12')
    leg.get_texts()[1].set_text('0.12 - 0.25')
    leg.get_texts()[2].set_text('0.25 - 0.50')
    leg.get_texts()[3].set_text('0.50 - 0.75')
    leg.get_texts()[4].set_text('0.75 - 1.00')
    
    
    leg.get_frame().set_linewidth(0.0)
    
    ax.set_axis_off()
    plt.show()
    
def plot_grad_GDF(gdf, column, title, black_back = True, legend = False, cmap = 'Greys_r'):
    
    f, ax = plt.subplots(1, figsize=(15, 15))
    gdf.plot(ax = ax, column = column, k = 7, cmap = cmap, linewidth = 1.2, scheme = 'fisher_jenks', legend = legend)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = gdf[column].min(), vmax = gdf[column].max()))
    f.suptitle(title)
    ax.set_axis_off()
    plt.axis('equal')
    if black_back == True: plt.rcParams['figure.facecolor'] = 'black'
    else: plt.rcParams['figure.facecolor'] = 'white'     
    if legend == True:
        leg = ax.get_legend()  
        leg.get_frame().set_linewidth(0.0)
    sm._A = []
    f.colorbar(sm)
    
    plt.show()

def scaling_columnDF(df, i, inverse = False):
    
    """
    it scales a column of a dataframe from 0 to 1
    
    Parameters
    df: pandas dataframe
    i: string (column name)
    ----------
    """
    df[i+'_sc'] = (df[i]-df[i].min())/(df[i].max()-df[i].min())
    if (inverse == True): df[i+'_sc'] = 1-(df[i]-df[i].min())/(df[i].max()-df[i].min())
        
	
def nodes_dict(G):
    """
    it creates a dictionary where keys represent the id of the node, the item is a tuple of coordinates
    
    Parameters
    G: NetworkX graph
    ----------
    """
    nodes_list = G.nodes()
    nodes_dict = {}

    for i, item in enumerate(nodes_list):
        cod = item
        x = nodes_list[item]['x']
        y = nodes_list[item]['y']
        nodes_dict[cod] = (x,y)
    
    return nodes_dict

def dict_to_df(list_dict, list_col):
    """
    take a list of dictionary and merge them in a df,
    
    Parameters
    list_dict: list of dictionaries that will become df columns
    list_col: list of the names that will be used as colums heading
    ----------
    """
    
    df = pd.DataFrame(list_dict).T
    df.columns = ['d{}'.format(i) for i, col in enumerate(df, 1)]
    df.columns = list_col
    
    return(df)

def dual_id_dict(dict_values, graph, nodeAttribute):
    """
    It could be used when one deals with a dual graph and wants to reconnect some analysis conducted on this representation to the analysis
    conducted on the primal graph. It takes for example the nodes betweennes centrality dictionary of the dual, and associate the result 
    with the original edgeID 
    
    Parameters
    ----------
    dict_values: dictionary of nodeID and centrality values (or other computation)
    graph: NetworkX graph
    nodeAttribute: string, attribute of the node
    """
    
    view = dict_values.items()
    ed_list = list(view)
    ed_dict = {}

    for p in ed_list: ed_dict[graph.node[p[0]][nodeAttribute]] = p[1] #Attribute and measure
        
    return(ed_dict)

# math functions for angle computations
# from Abhinav Ramakrishnan answer in https://stackoverflow.com/a/28261304/7375309

def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]

def ang_obsolete(lineA, lineB):
    # Get nicer vector form
    vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
    # Get dot prod
    dot_prod = dot(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    # Get cosine value
    cos_ = dot_prod/magA/magB
    # Get angle in radians and then convert to degrees
    angle = math.acos(dot_prod/magB/magA)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle)%360
    return ang_deg
     #if ang_deg-180<=0:
     #    return 180 - ang_deg
    # else:
     #    return ang_deg      
        
def ang(geolineA, geolineB, degree = False):

    coordsA = list(geolineA.coords)
    coordsB = list(geolineB.coords)   

    x_fA = float("{0:.10f}".format(coordsA[0][0]))
    y_fA = float("{0:.10f}".format(coordsA[0][1]))
    x_tA = float("{0:.10f}".format(coordsA[-1][0]))
    y_tA = float("{0:.10f}".format(coordsA[-1][1]))
    
    x_fB = float("{0:.10f}".format(coordsB[0][0]))
    y_fB = float("{0:.10f}".format(coordsB[0][1]))
    x_tB = float("{0:.10f}".format(coordsB[-1][0]))
    y_tB = float("{0:.10f}".format(coordsB[-1][1]))
    
    if ((x_tA, y_tA) == (x_tB, y_tB)):
        lineA = ((x_fA, y_fA),(x_tA,y_tA))
        lineB = ((x_tB, y_tB),(x_fB, y_fB))

    elif ((x_tA, y_tA) == (x_fB, y_fB)):
        lineA = ((x_fA, y_fA),(x_tA,y_tA))
        lineB = ((x_fB, y_fB),(x_tB, y_tB))

    elif ((x_fA, y_fA) == (x_fB, y_fB)):
        lineA = ((x_tA, y_tA),(x_fA,y_fA))
        lineB = ((x_fB, y_fB),(x_tB, y_tB))

    else: #(from_node == to_node2)
        lineA = ((x_tA, y_tA),(x_fA,y_fA))
        lineB = ((x_tB, y_tB),(x_fB, y_fB))
    
    
    # Get nicer vector form
    vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
    
    try:
        # Get dot prod
        dot_prod = dot(vA, vB)
        # Get magnitudes
        magA = dot(vA, vA)**0.5
        magB = dot(vB, vB)**0.5
        # Get cosine value
        cos_ = dot_prod/magA/magB
        # Get angle in radians and then convert to degrees
        angle_rad = math.acos(dot_prod/magB/magA)
        # Basically doing angle <- angle mod 360
        angle_deg = math.degrees(angle_rad)%360
        
    except:
        angle_deg = 0
        angle_rad = 0
        
    if degree == True: return angle_deg
    else: return angle_rad

def isparallel(geolineA, geolineB):
    
    angle = ang(geolineA, geolineB, degree = True)
    if ((angle <= 20) | (angle >= 160)): return True
    else: return False

def euclidean_distance(xs, ys, xt, yt):
    """ xs stands for x source and xt for x target """
    return sqrt((xs - xt)**2 + (ys - yt)**2)

	
# preparation functions
	
def get_fromOSM(type_download, place, network_type, epsg, distance = 7000, project = False): 
    """
    
    The function downloads and creates a simplified OSMNx graph for a selected area.
    Afterwards, geopandas dataframes for nodes and edges are created, assignind new nodeID and streeID identifiers.
    osmid are indeed heavy and confusing.
        
    Parameters
    place: string, name of cities or areas in OSM
    ----------
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
    
    if project == True: G = ox.project_graph(G)
    
    for i, item in enumerate(G.edges()):
        if isinstance(G[item[0]][item[1]][0]['osmid'], (list,)):
            G[item[0]][item[1]][0]['osmid'] = G[item[0]][item[1]][0]['osmid'][0]
            
            
    nodes = ox.graph_to_gdfs(G, nodes=True, edges=False, node_geometry=True, fill_edge_geometry=False)
    nodes_gdf = nodes.drop(nodes[['highway', 'ref']], axis=1)
    edges = ox.graph_to_gdfs(G, nodes=False, edges=True, node_geometry=False, fill_edge_geometry=True)
    edges_gdf = edges[['geometry', 'length', 'osmid', 'u','v', 'highway','key', 'oneway', 'maxspeed','name']]
    
    edges_gdf = edges_gdf.rename(columns = {'u':'old_u'})
    edges_gdf = edges_gdf.rename(columns = {'v':'old_v'})
    
    nodes_gdf = nodes_gdf.reset_index(drop=True)
    nodes_gdf['old_nodeID'] = nodes_gdf.osmid.astype('int64')
    nodes_gdf['nodeID'] = nodes_gdf.index.values.astype(int)

    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['old_nodeID', 'nodeID']], how='left', left_on="old_u", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'u'})
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['old_nodeID', 'nodeID']], how='left', left_on="old_v", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'v'})
    
    edges_gdf = edges_gdf.reset_index(drop=True)
    edges_gdf['streetID'] = edges_gdf.index.values.astype(int)
    
    nodes_gdf = nodes_gdf[['nodeID','x','y','geometry']]
    nodes_gdf.gdf_name = 'Nodes_gdf' #for OSMNx
    edges_gdf = edges_gdf[['streetID','u','v','key','geometry', 'length', 'highway','oneway', 'maxspeed','name']]
    edges_gdf['highway'] = [x[0] if type(x) == list else x for x in edges_gdf['highway']]
    edges_gdf['name'] = [x[0] if type(x) == list else x for x in edges_gdf['name']]
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
    roadType_field: string
    direction_field: string
    speed_field: string
    name_field: string
    area: int
    """
    
    #try reading street network

    streets = gpd.read_file(directory)
    streets = streets_gdf.to_crs(epsg=epsg)
    
    if (area != None):
        cn = streets.geometry.unary_union.centroid
        buffer = cn.buffer(area) 
        streets = streets[streets.geometry.within(buffer)]
        
    columns = [roadType_field, direction_field, speed_field, name_field]
    new_columns = ['highway','oneway', 'maxspeed','name']
    streets['from'] = "NaN"
    streets['to'] = "NaN"
    
    for n, i in enumerate(columns): 
        if (i is not None): streets[new_columns[n]] = streets_gdf[i]
     
    standard_columns = ['geometry', 'from', 'to']
    streets = streets_[standard_columns + [new_columns[n] for n, i in enumerate(columns) if i is not None]]
    index_geometry = streets.columns.get_loc("geometry")+1
    
    for row in streets.itertuples():
        line = []
        
        # to remove Z coordinates (assuming a planar graph)
        coord = list(row[index_geometry].coords)
        from_node = coord[0][0:2]
        to_node = coord[-1][0:2]
                
        for i in range(0, len(coord)):
            point = coord[i][0:2]
            line.append(point)

        t = LineString([coor for coor in line])
        streets.set_value(row[0],'geometry', t)
        streets.set_value(row[0], 'from', from_node)
        streets.set_value(row[0], 'to', to_node)
        
    streets = streets.loc[streets['from'] != streets['to']]
    unique_nodes_tmp = list(streets['to'].unique()) + list(streets['from'].unique())
    unique_nodes = list(set(unique_nodes_tmp))
    
#     if (simplify == True):
#         for i in unique_nodes:
#             tmp = streets[(streets['from'] == i) | (streets['to'] == i)].copy()
            
#             if(len(tmp) != 2): continue
#             else: 
#                 index_first = tmp.iloc[0].name
#                 index_second = tmp.iloc[1].name

#                 if (tmp.iloc[0]['from'] == tmp.iloc[1]['from']):
#                     streets.loc[index_first]['from'] = tmp.iloc[1]['to']
#                 elif (tmp.iloc[0]['from'] == tmp.iloc[1]['to']):
#                     streets.loc[index_first]['from'] = tmp.iloc[1]['from']
#                 elif (tmp.iloc[0]['to'] == tmp.iloc[1]['from']):
#                     streets.loc[index_first]['to'] = tmp.iloc[1]['to'] 
#                 else: #(tmp.iloc[0]['to'] == tmp.iloc[1]['to'])
#                     streets_gdf.loc[index_first]['to'] = tmp.iloc[1]['from']

#                 multi_line = MultiLineString([tmp.iloc[0]['geometry'], tmp.iloc[1]['geometry']])
#                 merged_line = linemerge(multi_line)
#                 streets_gdf.loc[index_first]['geometry'] = merged_line
#                 streets_gdf.drop([index_second], axis=0, inplace = True)
       
#     # getting unique nodes again
#     streets_gdf = streets_gdf.loc[streets_gdf['from'] != streets_gdf['to']]
#     unique_nodes_tmp = list(streets_gdf['to'].unique()) + list(streets_gdf['from'].unique())
#     unique_nodes = list(set(unique_nodes_tmp))
    
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
    edges_gdf['key'] = 0 #for OSMNx
    edges_gdf['length'] = gpd.GeoSeries(edges_gdf['geometry'].length)
    nodes_gdf.drop(['coordinates'], axis = 1, inplace = True)
    nodes_gdf.gdf_name = 'Nodes_gdf' #for OSMNx
        
    return(nodes_gdf, edges_gdf)
	

def reset_index_gdf(nodes_gdf, edges_gdf):
    edges_gdf = edges_gdf.rename(columns = {'u':'old_u'})
    edges_gdf = edges_gdf.rename(columns = {'v':'old_v'})
    
    nodes_gdf = nodes_gdf.reset_index(drop = False)
    nodes_gdf['old_nodeID'] = nodes_gdf['index'].astype('int64')
    nodes_gdf['nodeID'] = nodes_gdf.index.values.astype('int64')
    
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['old_nodeID', 'nodeID']], how='left', left_on="old_v", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'v'})

    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['old_nodeID', 'nodeID']], how='left', left_on="old_u", right_on="old_nodeID")
    edges_gdf = edges_gdf.rename(columns = {'nodeID':'u'})

    edges_gdf.drop(['old_u', 'old_nodeID_x', 'old_nodeID_y', 'old_v'], axis = 1, inplace = True)
    nodes_gdf.drop(['old_nodeID', 'index'], axis = 1, inplace = True)
    
    edges_gdf = edges_gdf.reset_index(drop=True)
    edges_gdf['streetID'] = edges_gdf.index.values.astype(int)
    
    return(nodes_gdf, edges_gdf)


def double_nodes(nodes_gdf, edges_gdf):
    
    nodes_gdf =  nodes_gdf.copy()
    edges_gdf = edges_gdf.copy()
    
    print('Eliminating duplicate geometries - nodes..')
    sindex = nodes_gdf.sindex
    index_geometry = nodes_gdf.columns.get_loc("geometry")+1
    
    # detect duplicate geometries
    G = nodes_gdf["geometry"].apply(lambda geom: geom.wkb)
    new_nodes = nodes_gdf.loc[G.drop_duplicates().index]
    
    to_edit = list(set(nodes_gdf.index.values.tolist()) - set((new_nodes.index.values.tolist())))
    if len(to_edit) == 0: return(nodes_gdf, edges_gdf)
    else:
        sindex = new_nodes.sindex
        for node in to_edit:
            # find new index
            geo = nodes_gdf.loc[node].geometry
            possible_matches_index = list(sindex.intersection(geo.bounds))
            possible_matches = new_nodes.iloc[possible_matches_index]
            precise_matches  = possible_matches[possible_matches.intersects(geo)]
            index = precise_matches.iloc[0].name
            edges_gdf.loc[edges_gdf.u == node,'u'] = index
            edges_gdf.loc[edges_gdf.v == node,'v'] = index
        return(new_nodes, edges_gdf)
    

def fix_dead_ends(nodes_gdf, edges_gdf):
    nodes_gdf =  nodes_gdf.copy()
    edges_gdf = edges_gdf.copy()
    
    print('Removing dead ends..')
    dd_u = dict(edges_gdf['u'].value_counts())
    dd_v = dict(edges_gdf['v'].value_counts())
    dd = {k: dd_u.get(k, 0) + dd_v.get(k, 0) for k in set(dd_u) | set(dd_v)}
    to_delete = {k: v for k, v in dd.items() if v == 1}
    
    if len(to_delete) == 0: return(nodes_gdf, edges_gdf)
    to_delete_list = list(to_delete.keys())
    nodes_gdf.drop(to_delete_list, axis = 0 , inplace = True)
    edges_gdf = edges_gdf[~edges_gdf['u'].isin(to_delete_list)]
    edges_gdf = edges_gdf[~edges_gdf['v'].isin(to_delete_list)]

    return(nodes_gdf, edges_gdf)

def nodes_simplified(edges_gdf):
    
    simplified = True
    dd_u = dict(edges_gdf['u'].value_counts())
    dd_v = dict(edges_gdf['v'].value_counts())
    dd = {k: dd_u.get(k, 0) + dd_v.get(k, 0) for k in set(dd_u) | set(dd_v)}
    to_edit = {k: v for k, v in dd.items() if v == 2}
    if len(to_edit) == 0: return(simplified)
    simplified = False
            
    return(simplified)

def edges_simplified(edges_gdf):
    
#     index_u = edges_gdf.columns.get_loc("u")+1
#     index_v = edges_gdf.columns.get_loc("v")+1
       
    simplified = True 
    edges_gdf['code'] = 'NA'
    edges_gdf['code'][edges_gdf['v'] >= edges_gdf['u']] = edges_gdf.u.astype(str)+"-"+edges_gdf.v.astype(str)
    edges_gdf['code'][edges_gdf['v'] < edges_gdf['u']] = edges_gdf.v.astype(str)+"-"+edges_gdf.u.astype(str)
    dd = dict(edges_gdf['code'].value_counts())
    dd = {k: v for k, v in dd.items() if v > 1}
    if len(dd) > 0: 
        print("Potential duplicate edges: ", len(dd))
        simplified = False
    
    
#     for row in edges_gdf.itertuples():

#         #                     if row[index_key] != rowC[index_key]: continue
#         #         #       if (row[index_type] != row[index_type]): continue

#         u = row[index_u]
#         v = row[index_v]
#         tmp = edges_gdf[(((edges_gdf.u == u) & (edges_gdf.v == v)) | ((edges_gdf.u == v) & (edges_gdf.v == u)))].copy()

#         for rowC in tmp.itertuples():
#             if (((row[0], rowC[0]) in to_ignore) | ((rowC[0], row[0]) in to_ignore)): tmp.drop(rowC[0], axis = 0, inplace = True)
        
#         if (len(tmp) == 1): continue # no double edges

#         else: 
#             simplified = False
#             break
    
    return(simplified)

def simplify_graph(nodes_gdf, edges_gdf):
    
    print('Removing pseudo-nodes..')
    nodes_gdf = nodes_gdf.copy()
    edges_gdf = edges_gdf.copy()
    edges_gdf = edges_gdf[edges_gdf.highway != 'primary_link']
    
    dd_u = dict(edges_gdf['u'].value_counts())
    dd_v = dict(edges_gdf['v'].value_counts())
    dd = {k: dd_u.get(k, 0) + dd_v.get(k, 0) for k in set(dd_u) | set(dd_v)}
    to_edit = {k: v for k, v in dd.items() if v == 2}
    if len(to_edit) == 0: return(nodes_gdf, edges_gdf)
    to_edit_list = list(to_edit.keys())
    
    for nodeID in to_edit_list:
        
        tmp = edges_gdf[(edges_gdf['u'] == nodeID) | (edges_gdf['v'] == nodeID)].copy()           
        index_first = tmp.iloc[0].name
        index_second = tmp.iloc[1].name
        
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

        new_line = line_coordsA + line_coordsB
        if edges_gdf.loc[index_first].u == edges_gdf.loc[index_first].v: 
            edges_gdf.drop([index_first, index_second], axis = 0, inplace = True)
            nodes_gdf.drop(nodeID, axis = 0, inplace = True)
            continue
            
        merged_line = LineString([coor for coor in new_line]) 
        edges_gdf.set_value(index_first, 'geometry', merged_line)
        
        edges_gdf.drop(index_second, axis = 0, inplace = True)
        nodes_gdf.drop(nodeID, axis = 0, inplace = True)
    return(nodes_gdf, edges_gdf)

def simplify_network(nodes_gdf, edges_gdf, dead_ends = False):
    
    nodes_gdf =  nodes_gdf.copy()
    edges_gdf = edges_gdf.copy()
    nodes_gdf, edges_gdf = double_nodes(nodes_gdf, edges_gdf)

    index_key = edges_gdf.columns.get_loc("key")+1
    index_type = edges_gdf.columns.get_loc("highway")+1
    index_u = edges_gdf.columns.get_loc("u")+1
    index_v = edges_gdf.columns.get_loc("v")+1
    index_geometry = edges_gdf.columns.get_loc("geometry")+1

    edges_gdf.sort_index(inplace = True)
    edges_gdf = edges_gdf[edges_gdf['u'] != edges_gdf['v']] #eliminate node-lines
   
    print('Cleaning and simplyfing network:')
    edges_gdf['code'] = 'NA'
    edges_gdf['coords'] = 'NA'
    cycle = 0
    
    while ((edges_simplified(edges_gdf) is False) | (nodes_simplified(edges_gdf) is False)):

        processed = []
        edges_gdf['length'] = edges_gdf['geometry'].length
        print('Cycle nr. ', cycle, '=============')
        cycle += 1
        
        edges_gdf['code'][edges_gdf['v'] >= edges_gdf['u']] = edges_gdf.u.astype(str)+"-"+edges_gdf.v.astype(str)
        edges_gdf['code'][edges_gdf['v'] < edges_gdf['u']] = edges_gdf.v.astype(str)+"-"+edges_gdf.u.astype(str)
        edges_gdf['coords'] = [list(x.coords) for x in edges_gdf.geometry]
        edges_gdf['coords'][(edges_gdf.u.astype(str)+"-"+edges_gdf.v.astype(str)) != edges_gdf.code] = [list(x.coords)[::-1] for x in edges_gdf.geometry]
        
        print('Eliminating duplicate geometries - edges..')
        G = edges_gdf['geometry'].apply(lambda geom: geom.wkb)
        edges_gdf = edges_gdf.loc[G.drop_duplicates().index]
                  
        edges_gdf['tmp'] = edges_gdf['coords'].apply(tuple, 1)
        edges_gdf.drop_duplicates(['tmp'], keep = 'first', inplace = True)
        
        print('Checking lines with same nodes..')
        dd = dict(edges_gdf['code'].value_counts())
        dd = {k: v for k, v in dd.items() if v > 1}
        for key, value in dd.items():
            tmp = edges_gdf[edges_gdf.code == key].copy()
            tmp.sort_values(['length'], ascending = True, inplace = True)
            u, v, geo, index_line = tmp.iloc[0]['u'], tmp.iloc[0]['v'], tmp.iloc[0]['geometry'], tmp.iloc[0].name, 

            for rowC in tmp.itertuples():
                if rowC[0] == index_line: continue
                uC, vC, geoC, index_lineC = rowC[index_u], rowC[index_v], rowC[index_geometry], rowC[0] 
                if geoC.length > (geo.length * 1.20):
                    edges_gdf.drop(index_lineC, axis = 0, inplace = True) 
                    continue
                else:
                    cl = center_line(u, v, uC, vC, geo, geoC)
                    edges_gdf.set_value(index_line,'geometry', cl)
                    edges_gdf.drop(index_lineC, axis = 0, inplace = True)

        if dead_ends is True: nodes_gdf, edges_gdf = fix_dead_ends(nodes_gdf, edges_gdf)
        nodes_gdf, edges_gdf = simplify_graph(nodes_gdf, edges_gdf)
        
    
    edges_gdf.drop(['code', 'coords', 'tmp'], axis = 1, inplace = True)
    nodes_gdf['nodeID'] = nodes_gdf.nodeID.astype(int)
    print("Done all =========")  
    return(nodes_gdf, edges_gdf)

# def is_t_junciont(nodes_gdf, edges_gdf, nodeID):
#     tmp = edges_gdf[edges_gdf.u == nodeID | edges_gdf.v == nodeID].copy()
#     if len(tmp) != 3: return(False)
    
#     if ((isparallel(tmp.iloc[0].geometry, tmp.iloc[1].geometry) == True) 
#         & (isperpendicular(tmp.iloc[0].geometry, tmp.iloc[2].geometry) == True))
         
#     elif ((isparallel(tmp.iloc[0].geometry, tmp.iloc[2].geometry) == True) & 
#          (isperpendicular(tmp.iloc[0].geometry, tmp.iloc[1].geometry) == True)):
#     elif ((isparallel(tmp.iloc[1].geometry, tmp.iloc[2].geometry) == True) &
#          (isperpendicular(tmp.iloc[0].geometry, tmp.iloc[1].geometry) == True)):

def simplify_junctions(nodes_gdf, edges_gdf, radius = 10):   
    
    buffered_nodes = nodes_gdf.buffer(radius).unary_union
    if isinstance(buffered_nodes, Polygon): buffered_nodes = [buffered_nodes]
        
    buffered_nodes_geoS = gpd.GeoSeries(list(buffered_nodes))
    buffered_nodes_df =  pd.concat([buffered_nodes_geoS.rename('geometry'), pd.Series(buffered_nodes_geoS.index).rename('code')], axis=1)

    buffered_nodes_gdf = gpd.GeoDataFrame(buffered_nodes_df, geometry = buffered_nodes_df.geometry)
    buffered_nodes_gdf['area']= buffered_nodes_gdf['geometry'].area
    buffered_nodes_gdf['centroid'] = buffered_nodes_gdf.geometry.centroid
    
    junctions_gdf = buffered_nodes_gdf[buffered_nodes_gdf["area"] > (radius*radius*3.14159)]
    junctions_gdf['x'], junctions_gdf['y'] = (junctions_gdf.geometry.centroid.x, junctions_gdf.geometry.centroid.y)
    junctions_gdf.index += nodes_gdf.index.max()+1
    junctions_gdf['code'] = junctions_gdf.index
    
    nodes_gdf['cluster'] = 'NA'
    index_geometry_bn = junctions_gdf.columns.get_loc("geometry")+1
    index_code = junctions_gdf.columns.get_loc("code")+1
    index_geometry = nodes_gdf.columns.get_loc("geometry")+1
    
    sindex = junctions_gdf.sindex    
    # set cluster column values
    for row in nodes_gdf.itertuples(): 
        geo = row[index_geometry]
        possible_matches_index = list(sindex.intersection(geo.buffer(5).bounds)) 
        if len(possible_matches_index) == 0: continue
        possible_matches = junctions_gdf.iloc[possible_matches_index]
        for rowC in possible_matches.itertuples():
            if geo.within(rowC[index_geometry_bn]): 
                nodes_gdf.set_value(row[0],'cluster', rowC[index_code])

    return(nodes_gdf, junctions_gdf)


def find_next_cluster(nodes_gdf, edges_gdf, index_line, search_direction):
    
    if search_direction == 'v': 
        v = edges_gdf.loc[index_line]['v']
        possible_matches = edges_gdf[edges_gdf.u == v].copy()
    else: 
        u = edges_gdf.loc[index_line]['u']
        possible_matches = edges_gdf[edges_gdf.v == u].copy()
    
    line = edges_gdf.loc[index_line].geometry
    name = edges_gdf.loc[index_line]['name']
    
    index_geometry = edges_gdf.columns.get_loc("geometry")+1
    index_name = edges_gdf.columns.get_loc("name")+1
    index_u = edges_gdf.columns.get_loc("u")+1
    index_v = edges_gdf.columns.get_loc("v")+1
    
    nodes_encountered = []
    lines_traversed = []
    line_coords = list(line.coords)
    if search_direction == 'u': line_coords.reverse()
        
    print('looing', edges_gdf.loc[index_line].streetID)
    found = False
    while (found is False):
        print(possible_matches.streetID.values.tolist(), 'indexes')
        if len(possible_matches) == 0: 
            found = True
            return(None, None, None, None, None)
        
        for p in possible_matches.itertuples():
            angle = ang(line, p[index_geometry])
            if isparallel(line, p[index_geometry]) is False:
                possible_matches.drop(p[0], axis = 0, inplace = True)
                continue
            else:
                uCP, vCP = p[index_u], p[index_v]
                
                if search_direction == 'v': cluster = nodes_gdf.loc[vCP].cluster
                else:  cluster = nodes_gdf.loc[uCP].cluster
                                    
                if cluster == 'NA':
                    lines_traversed.append(p[0])
                    if search_direction == 'v': 
                        possible_matches = edges_gdf[edges_gdf.u == vCP].copy()
                        nodes_encountered.append(uCP) 
                        line_coords = line_coords + list(p[index_geometry].coords)
                        
                    else:
                        possible_matches = edges_gdf[edges_gdf.v == uCP].copy()
                        nodes_encountered.append(vCP)
                        tmp = list(p[index_geometry].coords)
                        tmp.reverse()
                        line_coords = line_coords + tmp
                    break
                
                else:
                    found = True
                    lines_traversed.append(p[0])
                    if search_direction == 'v':
                        nodes_encountered.append(uCP)
                        last_node = vCP
                        line_coords = line_coords + list(p[index_geometry].coords)
                    else: 
                        nodes_encountered.append(vCP)
                        last_node = uCP
                        tmp = list(p[index_geometry].coords)
                        tmp.reverse()
                        line_coords = line_coords + tmp
                    break
    
    merged_line = LineString([coor for coor in line_coords])   
    return(cluster, merged_line, lines_traversed, nodes_encountered, last_node)
    
def center_line(u, v, uC, vC, line_geo, line_geoC): 
    
    line_coordsA = list(line_geo.coords)
    line_coordsB = list(line_geoC.coords)
        
    if ((u == vC) & (v == uC)): line_coordsB.reverse()  
    if line_coordsA == line_coordsB: center_line = LineString([coor for coor in line_coordsA]) 
    else:
        while len(line_coordsA) > len(line_coordsB):
            index = int(len(line_coordsA)/2)
            del line_coordsA[index]
        while len(line_coordsB) > len(line_coordsA):
            index = int(len(line_coordsB)/2)
            del line_coordsB[index]          
    
        new_line = line_coordsA
        for n, i in enumerate(line_coordsA):
            link = LineString([coor for coor in [line_coordsA[n], line_coordsB[n]]])
            np = link.centroid.coords[0]           
            new_line[n] = np
            
        new_line[0] = line_coordsA[0]
        new_line[-1] = line_coordsA[-1]
        center_line = LineString([coor for coor in new_line])

    return(center_line)
    
def center_line_cluster(line_geo, line_geoC, nodes_gdf, junctions_gdf, junction_from, junction_to, one_cluster = False): 
    
    coord_from = (junctions_gdf.loc[junction_from]['x'], junctions_gdf.loc[junction_from]['y'])
    if one_cluster == True: coord_from = (nodes_gdf.loc[node_from]['x'], nodes_gdf.loc[node_from]['y'])
    coord_to =  (junctions_gdf.loc[junction_to]['x'], junctions_gdf.loc[junction_to]['y'])
    line_coordsA = list(line_geo.coords)
    line_coordsB = list(line_geoC.coords)
        
    # no need to reverse lines, as they should arrive already in the same order      
    #different number of vertexes, connect the line
    while len(line_coordsA) > len(line_coordsB):
        index = int(len(line_coordsA)/2)
        del line_coordsA[index]
    while len(line_coordsB) > len(line_coordsA):
        index = int(len(line_coordsB)/2)
        del line_coordsB[index]      
    
    new_line = line_coordsA
    for n, i in enumerate(line_coordsA):
        link = LineString([coor for coor in [line_coordsA[n], line_coordsB[n]]])
        np = link.centroid.coords[0]           
        new_line[n] = np
        
    new_line[0] = coord_from
    new_line[-1] = coord_to
    center_line = LineString([coor for coor in new_line])           
        
    return(center_line)

def split_line_interpolation(node_index, line_geo, nodes_gdf, edges_gdf):
    
    point = nodes_gdf.loc[node_index].geometry
    old_list = list(line_geo.coords)
    
    starting_point = Point(line_geo.coords[0])
    np = nearest_points(point, line_geo)[1]
    distance_start = np.distance(starting_point)

    new_line_A = []
    new_line_B = []

    if(len(old_list)) == 2:
        new_line_A = [old_list[0],  np.coords[0]]
        new_line_B = [np.coords[0], old_list[-1]]
        lineA = LineString([coor for coor in new_line_A])
        lineB = LineString([coor for coor in new_line_B])

    else:
        new_line_A.append(old_list[0])
        new_line_B.append(np.coords[0])

        for n, i in enumerate(old_list):
            if (n == 0) | (n == len(old_list)-1): continue
            if (Point(i).distance(starting_point)) < distance_start:
                new_line_A.append(i)
            else:
                new_line_B.append(i)

        new_line_A.append(np.coords[0])
        new_line_B.append(old_list[-1])
        lineA = LineString([coor for coor in new_line_A])
        lineB = LineString([coor for coor in new_line_B])
    
    return((lineA, lineB), np)

def interpolate(u, center_line, last_node, nodes_to_interpolate, nodes_gdf, edges_gdf, index_line):
        
    line = center_line    
    new_index = index_line
    for counter, node in enumerate(nodes_to_interpolate):

        result, np = split_line_interpolation(node, line, nodes_gdf, edges_gdf)
              
        #first part of the segment, adjusting node coordinates
        nodes_gdf.set_value(node, 'x', np.coords[0][0])
        nodes_gdf.set_value(node, 'y', np.coords[0][1])
        nodes_gdf.set_value(node, 'geometry', np)
        
        if counter == 0: edges_gdf.set_value(new_index, 'u', u)
            
        edges_gdf.set_value(new_index, 'geometry', result[0])
        edges_gdf.set_value(new_index, 'v', node)
                
        # second part of the segment
        new_index = max(edges_gdf.index)+1
        
        edges_gdf.loc[new_index] = edges_gdf.loc[index_line]
        edges_gdf.set_value(new_index, 'geometry', result[1])
        edges_gdf.set_value(new_index, 'u', node)
        edges_gdf.set_value(new_index, 'v', last_node) 
        line = result[1]
            
    return(nodes_gdf, edges_gdf)

def interpolate_multi(u, center_line, last_node, list_nodes, nodes_gdf, edges_gdf, index_line):
    
    line = center_line  
    new_index = index_line
    distances = {}
    for node in list_nodes:
        distance = nodes_gdf.loc[node]['geometry'].distance(Point(center_line.coords[0]))
        distances[node] = distance

    distance_sorted = sorted(distances.items(), key=lambda kv: kv[1])               
    for counter, node in enumerate(distance_sorted):
        
        node = distance_sorted[counter][0]
        result, np = split_line_interpolation(node, line, nodes_gdf, edges_gdf)
        
        #first part of the segment, adjusting node coordinates
        nodes_gdf.set_value(node, 'x', np.coords[0][0])
        nodes_gdf.set_value(node, 'y', np.coords[0][1])
        nodes_gdf.set_value(node, 'geometry', np)

        if counter == 0: edges_gdf.set_value(new_index, 'u', u)
            
        edges_gdf.set_value(new_index, 'geometry', result[0])
        edges_gdf.set_value(new_index, 'v', node)
         
        # second part of the segment
        new_index = max(edges_gdf.index)+1
        
        edges_gdf.loc[new_index] = edges_gdf.loc[index_line]
        edges_gdf.set_value(new_index, 'geometry', result[1])
        edges_gdf.set_value(new_index, 'u', node)
        edges_gdf.set_value(new_index, 'v', last_node) 
        line = result[1]
            
    return(nodes_gdf, edges_gdf)

def simplify_dual_lines(nodes_gdf, edges_gdf, junctions_gdf):
    
    edges_gdf = edges_gdf.copy()
    nodes_gdf = nodes_gdf.copy()
    list_cluster = junctions_gdf.index.values.tolist()  
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['cluster', 'nodeID']], how= 'left', left_on= "u", right_on = "nodeID")
    edges_gdf = edges_gdf.rename(columns = {'cluster':'cluster_u'})
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['cluster', 'nodeID']], how= 'left', left_on= "v", right_on = "nodeID")
    edges_gdf = edges_gdf.rename(columns = {'cluster':'cluster_v'})  
    
    index_geometry = edges_gdf.columns.get_loc("geometry")+1
    index_cluster_u = edges_gdf.columns.get_loc("cluster_u")+1
    index_cluster_v = edges_gdf.columns.get_loc("cluster_v")+1
    index_u = edges_gdf.columns.get_loc("u")+1
    index_v = edges_gdf.columns.get_loc("v")+1
    index_name = edges_gdf.columns.get_loc("name")+1
    
    index_cluster = nodes_gdf.columns.get_loc("cluster")+1
    
    edges_gdf['cluster_uR'] = None
    edges_gdf['cluster_vR'] = None
    
    index_cluster_uR = edges_gdf.columns.get_loc("cluster_uR")+1
    index_cluster_vR = edges_gdf.columns.get_loc("cluster_vR")+1
    old_edges_gdf = edges_gdf.copy()

    for row in edges_gdf.itertuples():

        if row[index_cluster_u] == row[index_cluster_v]: continue
        if (row[index_cluster_u] == 'NA') & (row[index_cluster_v] == 'NA'): continue
        if (row[index_cluster_u] != 'NA') & (row[index_cluster_v] != 'NA'): continue
        elif row[index_cluster_u] == 'NA':
            result = find_next_cluster(nodes_gdf, edges_gdf, row[0], 'u')
            destination = result[0]
            edges_gdf.set_value(row[0], 'cluster_uR', destination)
        else:
            result = find_next_cluster(nodes_gdf, edges_gdf, row[0], 'v')
            destination = result[0]
            edges_gdf.set_value(row[0], 'cluster_vR', destination)
    
    ################################ FROM NODES TO CLUSTERED JUNCTIONS
    
    processed = []
    print('Simplifying dual lines: First part')
    junctions_gdf['keep'] = False
    
    for row in nodes_gdf.itertuples():
        if row[index_cluster] != 'NA': continue
        else:
            tmp = edges_gdf[((edges_gdf.u == row[0]) | (edges_gdf.v == row[0]))].copy()

            for rowC in tmp.itertuples():
                if rowC[0] in processed: continue 
                if rowC[index_u] == row[0]:
                    destination = rowC[index_cluster_v]
                    if destination == 'NA': destination = rowC[index_cluster_vR]
                elif rowC[index_v] == row[0]:
                    destination = rowC[index_cluster_u]
                    if destination == 'NA': destination = rowC[index_cluster_uR]
                if destination is None: continue


                group = tmp[((tmp.cluster_u == destination) | (tmp.cluster_uR == destination) 
                         | (tmp.cluster_v == destination) | (tmp.cluster_vR == destination))].copy()
                # orientate everything from "u" to "v"
                
                group['direction'] = 'v'
                for g in group.itertuples():
                    if g[index_v] == row[0]:
                        line_geometry = list(g[index_geometry].coords)
                        line_geometry.reverse() 
                        new_line = LineString([coor for coor in line_geometry])
                        old_u = g[index_u]
                        old_cluster_u = g[index_cluster_u]
                        old_cluster_uR = g[index_cluster_uR]

                        group.set_value(g[0],'geometry', new_line)
                        group.set_value(g[0],'u', g[index_v])
                        group.set_value(g[0],'v', old_u)
                        group.set_value(g[0],'cluster_u', g[index_cluster_v])
                        group.set_value(g[0],'cluster_v', old_cluster_u)
                        group.set_value(g[0],'cluster_uR', g[index_cluster_vR])
                        group.set_value(g[0],'cluster_vR', old_cluster_uR)
                        group.set_value(g[0], 'direction', 'u') # indicates original direction
                
                group = group[(group.cluster_v == destination) | (group.cluster_vR == destination)].copy()
                group = group[~group.index.isin(processed)]
                if len(group) == 1: break
                
                if len(group) == 2:
                    c_v, c_vC, = group.iloc[0]['cluster_v'], group.iloc[1]['cluster_v']
                    u, uC =  group.iloc[0]['u'], group.iloc[1]['u']
                    v, vC = group.iloc[0]['v'], group.iloc[1]['v']
                    geo, geoC = group.iloc[0]['geometry'], group.iloc[1]['geometry']
                    dr, drC = group.iloc[0]['direction'], group.iloc[1]['direction']
                    index_line, index_lineC  = group.iloc[0].name, group.iloc[1].name

                    if (c_v == c_vC) & (c_v != 'NA'):
                        destination = c_v
                        cl = center_line(u, v, uC, vC, geo, geoC)
                        edges_gdf.drop(index_lineC, axis = 0, inplace = True)
                        if dr == 'u':
                            line_geometry = list(cl.coords)
                            line_geometry.reverse() 
                            cl = LineString([coor for coor in line_geometry])
                        
                        edges_gdf.set_value(index_line, 'geometry', cl)
                        processed = processed + [index_line, index_lineC]
                        junctions_gdf.set_value(destination, 'keep', True)
                        break # next group

                    ######################################################## 
                    ## SUB-OPTION 2: only one reaches another cluster:

                    elif (c_v != 'NA') | (c_vC != 'NA'):

                        if c_v != 'NA': 
                            destination = c_v
                            found, geoC, lines_t, nodes_en, vC = find_next_cluster(nodes_gdf,edges_gdf, index_lineC, drC)
                            last_node = vC
                        else: 
                            destination = c_vC
                            found, geo, lines_t, nodes_en, v = find_next_cluster(nodes_gdf,edges_gdf, index_line, dr)
                            last_node = v
                        
                        cl =  center_line_cluster(geo, geoC, nodes_gdf, junctions_gdf, u, destination, one_cluster = True)
                        nodes_gdf, edges_gdf = interpolate(u, cl, last_node, nodes_en, nodes_gdf, edges_gdf, index_line)                                       
                        processed = processed + [index_line, index_lineC] + lines_t
                        lines_t.append(index_lineC)
                        edges_gdf.drop(lines_t, axis = 0, inplace = True, errors = 'ignore')
                        junctions_gdf.set_value(destination, 'keep', True)
                        break # next group                          

                    ####################################################### 
                    # SUB-OPTION 3: none reaches a cluster directly; comparing the first reached cluster
                    else:                    
                        dest, geo, lines_t, nodes_en, v = find_next_cluster(nodes_gdf, edges_gdf, index_line, dr)
                        destC, geoC, lines_tC, nodes_enC, vC = find_next_cluster(nodes_gdf, edges_gdf, index_lineC, drC)    

                        # the center line is built in relation to the variable cluster as 'u', or from_node --> to_node
                        cl =  center_line_ocluster(geo, geoC, nodes_gdf, junctions_gdf, u, dest, one_cluster = True)

                        # last node does not matter, as it will be reassigned to the relative cluster
                        list_nodes = nodes_en + nodes_enC
                        nodes_gdf, edges_gdf = interpolate_multi(u, cl, v, list_nodes, nodes_gdf, edges_gdf, index_line)  

                        lines_tC.append(index_lineC)
                        edges_gdf.drop(lines_tC, axis = 0, inplace = True, errors = 'ignore')
                        processed = processed + [index_line] + lines_tC
                        junctions_gdf.set_value(destination, 'keep', True)
                        break

    ################################ FROM CLUSTERED JUNCTIONS TO CLUSTERED JUNCTIONS
    
    print('Simplifying dual lines: Second part')
    processed = []
    
    for cluster in list_cluster:
        edges_tmp = edges_gdf[((edges_gdf.cluster_u == cluster) | (edges_gdf.cluster_v == cluster))].copy()
        edges_tmp = edges_tmp[edges_tmp.cluster_u != edges_tmp.cluster_v]
        if len(edges_tmp) == 1: continue
        
        for row in edges_tmp.itertuples():
            if row[0] in processed: continue                
            group = edges_tmp.copy()
            
            # eliminate unparallel lines
            for rowC in group.itertuples():
                if row[0] == rowC[0]: continue
                elif rowC[0] in processed: 
                    group.drop(rowC[0], axis = 0, inplace = True)
                    continue
                elif ((isparallel(row[index_geometry], rowC[index_geometry]) is True) |
                      (row[index_name] == rowC[index_name])): continue
                else: group.drop(rowC[0], axis = 0, inplace = True)

            # does the line considered in the loop reach a cluster? if not straight away, at some point?
            
            group['direction'] = 'v'
            # orientate everything from "u" to "v"
            for rowC in group.itertuples():
                if rowC[index_cluster_v] == cluster:
                    line_geometry = list(rowC[index_geometry].coords)
                    line_geometry.reverse() 
                    new_line = LineString([coor for coor in line_geometry])
                    old_u = rowC[index_u]
                    old_cluster_u = rowC[index_cluster_u]
                    old_cluster_uR = rowC[index_cluster_uR]

                    group.set_value(rowC[0],'geometry', new_line)
                    group.set_value(rowC[0],'u', rowC[index_v])
                    group.set_value(rowC[0],'v', old_u)
                    group.set_value(rowC[0],'cluster_u', rowC[index_cluster_v])
                    group.set_value(rowC[0],'cluster_v', old_cluster_u)
                    group.set_value(rowC[0],'cluster_uR', rowC[index_cluster_vR])
                    group.set_value(rowC[0],'cluster_vR', old_cluster_uR)
                    group.set_value(rowC[0], 'direction', 'u') # indicates original direction

            if row[index_cluster_v] != 'NA': group_destination = row[index_cluster_v]
            else: group_destination = row[index_cluster_vR]                
            if group_destination == None: continue

            for rowC in group.itertuples():
                if rowC[index_cluster_v] != 'NA': secondary_destination = rowC[index_cluster_v]
                else: secondary_destination = rowC[index_cluster_vR]
                if (secondary_destination != group_destination): 
                    group.drop(rowC[0], axis = 0, inplace = True)
            
            ######################################################## OPTION 1
            
            if len(group) == 1: 
                break # no parallel streets to row[0] 
                            
            ######################################################## OPTION 2
             
            elif len(group) == 2:
                c_u, c_uC = group.iloc[0]['cluster_u'], group.iloc[1]['cluster_u']
                c_v, c_vC, = group.iloc[0]['cluster_v'], group.iloc[1]['cluster_v']
                u, uC =  group.iloc[0]['u'], group.iloc[1]['u']
                v, vC = group.iloc[0]['v'], group.iloc[1]['v']
                dr, drC = group.iloc[0]['direction'], group.iloc[1]['direction']
                geo, geoC = group.iloc[0]['geometry'], group.iloc[1]['geometry']
                index_line, index_lineC  = group.iloc[0].name, group.iloc[1].name
                
                ######################################################## 
                ## SUB-OPTION 1: they all reach another cluster:
                    
                if (c_v == c_vC) & (c_v != 'NA'):
                    destination = c_v
                    cl = center_line_cluster(geo, geoC, nodes_gdf, junctions_gdf, cluster, destination)
                    edges_gdf.drop(index_lineC, axis = 0, inplace = True)
                    if dr == 'u':
                        line_geometry = list(cl.coords)
                        line_geometry.reverse() 
                        cl = LineString([coor for coor in line_geometry])
                    
                    edges_gdf.set_value(index_line, 'geometry', cl)
                    processed = processed + [index_line, index_lineC]
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break # next group

                ######################################################## 
                ## SUB-OPTION 2: only one reaches another cluster:
                    
                elif (c_v != 'NA') | (c_vC != 'NA'):
                    
                    if c_v != 'NA': 
                        destination = c_v
                        found, geoC, lines_t, nodes_en, vC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)
                        last_node = vC
                    else: 
                        destination = c_vC
                        found, geo, lines_t, nodes_en, v = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                        last_node = v
                    
                    cl = center_line_cluster(geo, geoC, nodes_gdf, junctions_gdf, cluster, destination)
                    nodes_gdf, edges_gdf = interpolate(u, cl, last_node, nodes_en, nodes_gdf, edges_gdf, index_line)                                       
                    processed = processed + [index_line, index_lineC] + lines_t
                    lines_t.append(index_lineC)
                    edges_gdf.drop(lines_t, axis = 0, inplace = True, errors = 'ignore')
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break # next group                          
                    
                ####################################################### 
                # SUB-OPTION 3: none reaches a cluster directly; comparing the first reached cluster
                else: 
                    dest, geo, lines_t, nodes_en, v = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                    destC, geoC, lines_tC, nodes_enC, vC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)    
                        
                    # the center line is built in relation to the variable cluster as 'u', or from_node --> to_node
                    cl = center_line_cluster(geo, geoC, nodes_gdf, junctions_gdf, cluster, dest)
                        
                    # last node does not matter, as it will be reassigned to the relative cluster
                    list_nodes = nodes_en + nodes_enC
                
                    nodes_gdf, edges_gdf = interpolate_multi(u, cl, v, list_nodes, nodes_gdf, edges_gdf, index_line)  
                    
                    lines_tC.append(index_lineC)
                    edges_gdf.drop(lines_tC, axis = 0, inplace = True, errors = 'ignore')
                    processed = processed + [index_line] + lines_tC
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break
            
            ####################################################### OPTION 3
            
            elif len(group) == 3:

                c_u, c_uC, c_uCC = group.iloc[0]['cluster_u'], group.iloc[1]['cluster_u'], group.iloc[2]['cluster_u']
                c_v, c_vC, c_vCC = group.iloc[0]['cluster_v'], group.iloc[1]['cluster_v'], group.iloc[2]['cluster_v']
                u, uC, uCC =  group.iloc[0]['u'], group.iloc[1]['u'], group.iloc[2]['u']
                v, vC, vCC = group.iloc[0]['v'], group.iloc[1]['v'], group.iloc[2]['v']
                dr, drC, drCC = group.iloc[0]['direction'], group.iloc[1]['direction'], group.iloc[2]['direction']
                geo, geoC, geoCC = group.iloc[0]['geometry'], group.iloc[1]['geometry'], group.iloc[2]['geometry']
                index_line, index_lineC, index_lineCC  = group.iloc[0].name, group.iloc[1].name, group.iloc[2].name            
               
                ######################################################## 
                ## SUB-OPTION 1: they all reach another cluster (the same)
                    
                if ((c_v == c_vC) & (c_v == c_vCC) & (c_v != 'NA')):
                    max_dist = 0                        
                    to_delete = []
                    for g in group.itertuples():
                        for gC in group.itertuples():
                            if g[0] == gC[0]: continue
                            distance = (g[index_geometry].centroid).distance(gC[index_geometry].centroid)
                            if distance > max_dist: 
                                max_dist = distance
                                to_delete = [g[0],gC[0]] #the other is the center line
                        
                    group = group[~group.index.isin(to_delete)]
                    edges_gdf.drop(to_delete, axis = 0, inplace = True, errors = 'ignore')
                    # no need to change geometry here

                    processed = processed + [group.iloc[0].name] + to_delete
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break # next group
                        
                ########################################################  
                ## SUB-OPTION 2: two reach another cluster:   
                elif (((c_v == c_vC) & (c_v != 'NA'))| ((c_v == c_vCC) & (c_v != 'NA')) | ((c_vC == c_vCC) & (c_vC != 'NA'))):
                    # the one that doesn't reach the cluster is shorter than the other 2 segments. 
                    # let's keep the one that is the central one (closest to the shortest)
                    
                    group.sort_values(['length'], ascending = False, inplace = True)
                    distance_a = (group.iloc[0]['geometry']).distance(group.iloc[2]['geometry'])
                    distance_b = (group.iloc[1]['geometry']).distance(group.iloc[2]['geometry'])
                    
                    if distance_a > distance_b: 
                            processed.append(group.iloc[0].name)
                            edges_gdf.drop([group.iloc[0].name], axis = 0, inplace = True)
                            group.drop([group.iloc[0].name], axis = 0, inplace = True)                           
                    else: 
                            processed.append(group.iloc[1].name)
                            edges_gdf.drop([group.iloc[1].name], axis = 0, inplace = True)
                            group.drop([group.iloc[1].name], axis = 0, inplace = True)                           
                            
                    c_u, c_uC = group.iloc[0]['cluster_u'], group.iloc[1]['cluster_u']
                    c_v, c_vC = group.iloc[0]['cluster_v'], group.iloc[1]['cluster_v']
                    u, uC =  group.iloc[0]['u'], group.iloc[1]['u']
                    v, vC = group.iloc[0]['v'], group.iloc[1]['v']
                    dr, drC = group.iloc[0]['direction'], group.iloc[1]['direction']
                    geo, geoC = group.iloc[0]['geometry'], group.iloc[1]['geometry']
                    index_line, index_lineC  = group.iloc[0].name, group.iloc[1].name
                    
                    destination, line, lines_t, nodes_en, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)
                    cl = geo                                    
                    # last node does not matter, as it will be reassigned to the relative cluster
                    # the line used as "merged" is the already existing one, that arrives at the cluster
                    nodes_gdf, edges_gdf = interpolate(u, cl, last_node, nodes_en, nodes_gdf, edges_gdf, index_line)
                    
                    lines_t = lines_t + [index_lineCC, index_lineC]
                    edges_gdf.drop(lines_t, axis = 0, inplace = True, errors = 'ignore')
                    processed = processed + [index_line] + lines_t
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break
                
                ########################################################  
                ## SUB-OPTION 4: only one reaches a cluster:

                elif (c_v != 'NA')| (c_vC != 'NA') | (c_vCC != 'NA'):

                    if (c_v != 'NA'):
                        destC, geoC, lines_tC, nodes_enC, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)
                        destCC, geoCC, lines_tCC, nodes_enCC, last_nodeCC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineCC, drCC)
                        lines_t, nodes_en = [], []
                    elif (c_vC != 'NA'):
                        dest, geo, lines_t, nodes_en, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                        destCC, geoCC, lines_tCC, nodes_enCC, last_nodeCC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineCC, drCC)
                        lines_tC, nodes_enC = [], []
                    else:
                        dest, geo, lines_t, nodes_en, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                        destC, geoC, lines_tC, nodes_enC, last_nodeC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)                   
                        lines_tCC, nodes_enCC = [], []
                        
                        # exclude the 2 lines furter away                          
                    
                    max_dist = 0
                    dict_lines = {index_line: geo, index_lineC: geoC, index_lineCC: geoCC}
                    secondary_lines = []
                    
                    for key, value in dict_lines.items():
                        for keyC, valueC in dict_lines.items():
                            if key == keyC: continue
                            distance = (value).distance(valueC)
                            if distance > max_dist: 
                                max_dist = distance
                                secondary_lines = [key, keyC]
                    
                    list_nodes = nodes_en + nodes_enC + nodes_enCC
                    central = [x for x in list(dict_lines.keys()) if x not in secondary_lines][0]

                    if central == index_line: cl = geo
                    elif central == index_lineC: cl = geoC
                    else: cl = geoCC
                    
                    to_drop = secondary_lines + lines_t + lines_tC + lines_tCC
                    nodes_gdf, edges_gdf = interpolate_multi(u, cl, last_node, list_nodes, nodes_gdf, edges_gdf, index_line)     
                    edges_gdf.drop(to_drop, axis = 0, inplace = True, errors = 'ignore')
                    processed = processed + to_drop + [central]
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break
                 
                ########################################################  
                ## SUB-OPTION 4: none reaches a cluster:
                else: 
                    dest, geo, lines_t, nodes_en, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                    destC, geoC, lines_tC, nodes_enC, last_nodeC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)
                    destCC, geoCC, lines_tCC, nodes_enCC, last_nodeCC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineCC, drCC)  
                    # exclude the 2 lines furter away                          
                    
                    max_dist = 0
                    dict_lines = {index_line: geo, index_lineC: geoC, index_lineCC: geoCC}
                    secondary_lines = []
                    
                    for key, value in dict_lines.items():
                        for keyC, valueC in dict_lines.items():
                            if key == keyC: continue
                            distance = (value).distance(valueC)
                            if distance > max_dist: 
                                max_dist = distance
                                secondary_lines = [key, keyC]
                         
                    list_nodes = nodes_en + nodes_enC + nodes_enCC
                    central = [x for x in list(dict_lines.keys()) if x not in secondary_lines][0]

                    if central == index_line: cl = geo
                    elif central == index_lineC: cl = geoC
                    else: cl = geoCC
                    
                    to_drop = secondary_lines + lines_t + lines_tC + lines_tCC
                    nodes_gdf, edges_gdf = interpolate_multi(u, cl, last_node, list_nodes, nodes_gdf, edges_gdf, index_line)     
                    edges_gdf.drop(to_drop, axis = 0, inplace = True, errors = 'ignore')
                    processed = processed + to_drop + [central]

                    junctions_gdf.set_value(cluster, 'keep', True)
                    break
    
    edges_gdf.drop(['nodeID_x', 'nodeID_y','cluster_uR', 'cluster_vR'], axis = 1, inplace = True)
    edges_gdf['streetID'] = edges_gdf.index.values.astype(int)
    nodes_gdf['nodeID'] = nodes_gdf.index.values.astype(int)
    return(nodes_gdf, edges_gdf)

def correct_edges(nodes_gdf, edges_gdf, junctions_gdf):
    
    print("Correcting edges coordinates..")

    edges_gdf = edges_gdf.rename(columns = {'u':'old_u'})
    edges_gdf = edges_gdf.rename(columns = {'v':'old_v'})
    
    edges_gdf['u'] = 0
    edges_gdf['v'] = 0
    index_u = edges_gdf.columns.get_loc("u")+1
    index_v = edges_gdf.columns.get_loc("v")+1 
    index_old_u = edges_gdf.columns.get_loc("old_u")+1
    index_old_v = edges_gdf.columns.get_loc("old_v")+1
    index_geometry = edges_gdf.columns.get_loc("geometry")+1 
    
    index_x = junctions_gdf.columns.get_loc("x")+1
    index_y = junctions_gdf.columns.get_loc("y")+1
    index_centroid = junctions_gdf.columns.get_loc("centroid")+1
    index_check = junctions_gdf.columns.get_loc("keep")+1
    edges_gdf = edges_gdf[edges_gdf.cluster_u != edges_gdf.cluster_v]
    
    for row in edges_gdf.itertuples():

        line_list = list(row[index_geometry].coords)

        u = nodes_gdf.loc[row[index_old_u]]["cluster"]
        v = nodes_gdf.loc[row[index_old_v]]["cluster"]
        old_u = row[index_old_u]
        old_v = row[index_old_v]
               
        if ((u != 'NA') & (v !='NA')):  # change starting and ending node in the list of coordinates for the line
                if (junctions_gdf.loc[u].keep == False) & (junctions_gdf.loc[v].keep == False): 
                    u = old_u
                    v = old_v
                    line_list[0] = (nodes_gdf.loc[old_u]['x'], nodes_gdf.loc[old_u]['y'])
                    line_list[-1] = (nodes_gdf.loc[old_v]['x'], nodes_gdf.loc[old_v]['y'])
                elif junctions_gdf.loc[v].keep == False:
                    v = old_v
                    line_list[0] = (junctions_gdf.loc[u]['x'], junctions_gdf.loc[u]['y'])
                    line_list[-1] = (nodes_gdf.loc[old_v]['x'], nodes_gdf.loc[old_v]['y'])
                elif junctions_gdf.loc[u].keep == False:  
                    u = old_u    
                    line_list[0] = (nodes_gdf.loc[old_u]['x'], nodes_gdf.loc[old_u]['y'])
                    line_list[-1] = (junctions_gdf.loc[v]['x'], junctions_gdf.loc[v]['y'])
                else:
                    line_list[0] = (junctions_gdf.loc[u]['x'], junctions_gdf.loc[u]['y'])
                    line_list[-1] = (junctions_gdf.loc[v]['x'], junctions_gdf.loc[v]['y'])

        elif ((u == 'NA') & (v =='NA')):  # maintain old_u and old_v
                u = old_u
                v = old_v
                line_list[0] = (nodes_gdf.loc[old_u]['x'], nodes_gdf.loc[old_u]['y'])
                line_list[-1] = (nodes_gdf.loc[old_v]['x'], nodes_gdf.loc[old_v]['y'])

        elif ((u == 'NA') & (v != 'NA')) : # maintain old_u
                u = old_u
                line_list[0] = (nodes_gdf.loc[old_u]['x'], nodes_gdf.loc[old_u]['y'])
                
                if junctions_gdf.loc[v].keep == False:
                    v = old_v
                    line_list[-1] = (nodes_gdf.loc[v]['x'], nodes_gdf.loc[v]['y'])
                else:
                    line_list[-1] = (junctions_gdf.loc[v]['x'], junctions_gdf.loc[v]['y'])

        else: #(( u =! 'NA') & (v == 'NA') !: # maintain old_v
                v = old_v
                line_list[-1] = (nodes_gdf.loc[old_v]['x'], nodes_gdf.loc[old_v]['y'])
                if junctions_gdf.loc[u].keep == False:
                    u = old_u
                    line_list[0] = (nodes_gdf.loc[u]['x'], nodes_gdf.loc[u]['y'])
                else:
                    line_list[0] = (junctions_gdf.loc[u]['x'], junctions_gdf.loc[u]['y'])

        line_geo = (LineString([coor for coor in line_list]))
        
        if u == v: 
            edges_gdf.drop(row[0], axis = 0, inplace = True)
            continue
            
        edges_gdf.set_value(row[0],"u", u)
        edges_gdf.set_value(row[0],"v", v)
        edges_gdf.set_value(row[0],"geometry", line_geo)

    edges_gdf.drop(['old_u', 'old_v'], axis = 1, inplace=True)
    edges_gdf['u'] = edges_gdf['u'].astype(int)
    edges_gdf['v'] = edges_gdf['v'].astype(int)
    
#     nodes_gdf = nodes_gdf[nodes_gdf['cluster'] == "NA"]
    nodes_gdf['x'] = nodes_gdf['x'].astype(float)
    nodes_gdf['y'] = nodes_gdf['y'].astype(float)

    for row in junctions_gdf.itertuples():
        if row[index_check] != True: continue
        nodes_gdf.set_value(row[0], 'x', row[index_x])
        nodes_gdf.set_value(row[0], 'y', row[index_y])
        nodes_gdf.set_value(row[0], 'geometry', row[index_centroid])
        nodes_gdf.set_value(row[0], 'nodeID', row[0])
        nodes_gdf.set_value(row[0], 'cluster', 'NA')
    
    nodes_gdf.drop(['cluster'], axis = 1, inplace = True)
    print("Done") 
    return(nodes_gdf, edges_gdf)


def graph_fromGDF(nodes_gdf, edges_gdf, nodes_attributes, edges_costs):
    """
    It creates from due geopandas dataframes (street junctions and street segments) a NetworkX graph, passing by a OSMnx function.
    The length of the segments and their ID is stored in the edges attributes, as well as the nodeID.
    
    Parameters
    ----------
    nodes_gdf: geopandas dataframe
    edges_gdf: geopandas dataframe    
    """
    
    nodes_gdf.gdf_name = 'Nodes_list' #for OSMNx
    edges_gdf.gdf_name = 'Edges_list' #for OSMNx
    edges_gdf['key'] = 0
    G = ox.gdfs_to_graph(nodes_gdf, edges_gdf)
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



# centroids for dual analysis


def dual_gdf(nodes_gdf, edges_gdf, crs):
    """
    It creates two dataframes that are supposed to generate the dual graph of a street network. The nodes_dual gdf contains edges 
    centroids, the edges_dual gdf contains instead links between the street segment centroids. Those dual edges link real street segment 
    that share a junction.
    The centroids are stored with the original edge streetID, while the dual edges are associated with several attributes computed on the 
    original street segments (distance between centroids, deflection angle).
    
    Parameters
    ----------
    nodes_gdf: geopandas dataframe
    edges_gdf: geopandas dataframe  
    crs: dictionary
    """
    
    centroids_gdf = edges_gdf.copy()
    centroids_gdf['centroid'] = centroids_gdf['geometry'].centroid
    centroids_gdf['intersecting'] = 'NA'
    
    index_u = centroids_gdf.columns.get_loc("u")+1
    index_v = centroids_gdf.columns.get_loc("v")+1
    index_streetID = centroids_gdf.columns.get_loc("streetID")+1
         
    #find_intersecting segments
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
                intersections.append(p[index_streetID])  #appending streetID
                processed.append((p[0],c[0]))
    
        centroids_gdf.set_value(c[0],'intersecting', intersections)
        
    #creating vertexes representing street segments (centroids)
        
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
    index_geometry = nodes_dual.columns.get_loc("geometry")+1

    for row in nodes_dual.itertuples():
        
        streetID = row[index_streetID_nd] #streetID of the relative segment
        length = row[index_length]

        for i in list(row[index_intersecting]): #intersecting segments

            # i is the streetID
            length_i =  nodes_dual['length'][nodes_dual.streetID == i][i]
            distance = (length+length_i)/2
        
            # adding a row with u-v, key fixed as 0, Linestring geometry 
            # from the first centroid to the centroid intersecting segment 
            ls = LineString([row[index_geometry], nodes_dual.loc[i]['geometry']])
            
            edges_dual.loc[-1] = [streetID, i, 0, ls, distance] 
            edges_dual.index = edges_dual.index + 1
            
    edges_dual = edges_dual.sort_index(axis=0)
    geometry = edges_dual['geometry']
    edges_dual = gpd.GeoDataFrame(edges_dual[['u','v', 'key', 'length']], crs=crs, geometry=geometry)

    #computing deflection angle
    
    index_line = edges_dual.columns.get_loc("u")+1
    index_lineC = edges_dual.columns.get_loc("v")+1
    
    for row in edges_dual.itertuples():

        #retrieveing original lines from/to
        geo = edges_gdf[edges_gdf.index == row[index_line]].geometry.iloc[0]
        geoC = edges_gdf[edges_gdf.index == row[index_lineC]].geometry.iloc[0]
        
        deflection = ang(lineA, lineB, degree = True)
        deflection_rad = ang(lineA, lineB)
           
        edges_dual.set_value(row[0],'deg', deflection)
        edges_dual.set_value(row[0],'rad', deflection_rad)
        
    return(nodes_dual, edges_dual)

def get_dual_graph(nodes_dual, edges_dual, edge_costs):
    """
    The function generates a NetworkX graph from nodes and edges dual geopandas dataframes.
    
        
    Parameters
    ----------
    nodes_dual: geopandas dataframe
    edges_dual: geopandas dataframe  

    """
   
    nodes_dual.gdf_name = 'Dual_list'
    Gr = ox.gdfs_to_graph(nodes_dual, edges_dual)
    
    n = Gr.nodes()
    pos = {}

    for l, item in enumerate(n): pos[l] = (n[l]['x'], n[l]['y'], n[l]['streetID'])
        
    DG = nx.Graph() #Empty graph
    DG = DG.to_undirected()
    DG.add_nodes_from(pos.keys()) #Add nodes preserving coordinates
    
    for i, item in enumerate(DG.nodes()):
        DG.node[i]['x']=pos[i][0]
        DG.node[i]['y']=pos[i][1]
        DG.node[i]['streetID']=pos[i][2]
        
    for i, item in enumerate(Gr.edges()):
        DG.add_edge(item[0], item[1])
        for cost in costs: DG[item[0]][item[1]][cost] = Gr[item[0]][item[1]][0][cost]
        
    return(DG)

	
# natural roads extraction: the function has to be called twice for each segment, to check in both directions

def natural_roads(streetID, natural_id, direction, nodes_gdf, edges_gdf,): 
    """
    This function takes a direction "to" or "from" and two geopandas dataframes, one for roads one for nodes (or junctions).
    The dataframes are supposed to be simplified and can be obtained via the functions "get_fromOSM(place)" or "get_fromSHP(directory, 
    epsg)". 
    
    Natural road is a concept presented here and here and regards the cognitive perception/representation of road entities, regardless 
    changes in names or regardless interruptions. Rather, different street segments are mentally merged according to continuity rules (here 
    based on the deflection angle and the egoistic choice process, see:...)
    
    It takes the ID of the segment processed, the ID of the natural road (has to be controlled in a for loop, see example below)
    
    
    Parameters
    ----------
    streetID: int
    natural_id: int
    direction: string
    edges_gdf: geopandas dataframe
    nodes_gdf: geopandas dataframe
    """
        
    angles = {}
    directions_dict = {}
      
    index_geometry = edges_gdf.columns.get_loc("geometry")+1
    index_u = edges_gdf.columns.get_loc("geometry")+1
    index_v = edges_gdf.columns.get_loc("geometry")+1
    
    to_node = edges_gdf.loc[streetID]['u']
    from_node = edges_gdf.loc[streetID]['v']
    geo = edges_gdf.loc[streetID]['v']
    
        
    #assuming nodeID = ID dataframe
    x_t = float("{0:.10f}".format(nodes_gdf.loc[to_node]['x']))
    y_t = float("{0:.10f}".format(nodes_gdf.loc[to_node]['y']))
    x_f = float("{0:.10f}".format(nodes_gdf.loc[from_node]['x']))
    y_f = float("{0:.10f}".format(nodes_gdf.loc[from_node]['y']))
      
    #continue from the to_node     
    if (direction == "to"): intersecting = edges_gdf[(roads_gdf['u'] == to_node) | (edges_gdf['v'] == to_node)]
    
    #continue from the from_node
    else: intersecting = edges_gdf[(edges_gdf['u'] == from_node) | (edges_gdf['v'] == from_node)]

    if (len(intersecting) == 0): return
    
    for row_F in intersecting.itertuples():
        if ((streetID == row_F[0]) | (row_F[-1] != 'NA')): continue
        
        to_node_F = row_F[edges_gdf.columns.get_loc("u")]
        from_node_F = row_F[edges_gdf.columns.get_loc("v")]
        geo_F = row_F[edges_gdf.columns.get_loc("v")]


        if (to_node == to_node_F): towards = "fr"
        elif (to_node == from_node_F): towards = "to"
        elif (from_node == from_node_F): towards = "to"
        else: towards = "fr"

        deflection = ang(geo, geo_F, degree = True)

        if (deflection >= 45): continue
        else:
            angles[row_F[0]] = deflection
            directions_dict[row_F[0]] = towards
    
    #no natural continuations
    if (len(angles) == 0): 
        edges_gdf.set_value(streetID, 'natural_id', natural_id)
        return
    else:
        angles_sorted = sorted(angles, key = angles.get)
        matchID = angles_sorted[0]
        edges_gdf.set_value(id_road, 'natural_id', natural_id)
        natural_roads(matchID, natural_id, directions_dict[matchID], nodes_gdf, edges_gdf)

    
# centrality functions

def straightness_centrality(G, weight, normalized = True):
    """
        
    Parameters
    ----------
    """
    path_length = functools.partial(nx.single_source_dijkstra_path_length, weight=weight)

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
               
            if normalized:
                straightness_centrality[n] = straightness * (1.0/(len(G)-1.0) )

        else:
            straightness_centrality[n]=0.0

    return straightness_centrality

def reach_centrality(G, weight, radius):
    """
        
    Parameters
    ----------
    """
    path_length = functools.partial(nx.single_source_dijkstra_path_length, weight=weight)

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
                    weight_target = G.node[target]['weight']
                    reach = reach + weight_target
                        

            reach_centrality[n] = reach

        else:               
            reach_centrality[n]=0.0

    return reach_centrality

def local_betweenness_centrality(G, weight, radius):
    """
        
    Parameters
    ----------
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
	


    
