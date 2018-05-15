import osmnx as ox, networkx as nx, matplotlib.cm as cm, pandas as pd, numpy as np
import geopandas as gpd
import functools
import community
import math

from scipy import sparse
from scipy.sparse import linalg

from shapely.geometry import Point, LineString, Polygon, MultiPolygon, mapping, MultiLineString
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

def ang(lineA, lineB):
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
        
def ang_rad(lineA, lineB):
    """
    to get the angle in radian
    
    """
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

    return angle

def euclidean_distance(xs, ys, xt, yt):
    """ xs stands for x source and xt for x target """
    return sqrt((xs - xt)**2 + (ys - yt)**2)
	
# preparation functions
	
def get_fromOSM(place): 
    """
    
    The function downloads and creates a simplified OSMNx graph for a selected area.
    Afterwards, geopandas dataframes for nodes and edges are created, assignind new nodeID and streeID identifiers.
    osmid are indeed heavy and confusing.
        
    Parameters
    place: string, name of cities or areas in OSM
    ----------
    """
    G = ox.graph_from_place(place, network_type='all', simplify=True)
    G = ox.project_graph(G)
    
    for i, item in enumerate(G.edges()):
        if isinstance(G[item[0]][item[1]][0]['osmid'], (list,)):
            G[item[0]][item[1]][0]['osmid'] = G[item[0]][item[1]][0]['osmid'][0]
            
            
    nodes = ox.graph_to_gdfs(G, nodes=True, edges=False, node_geometry=True, fill_edge_geometry=False)
    nodes = nodes.drop(nodes[['highway', 'ref', 'lat', 'lon']], axis=1)
    edges = ox.graph_to_gdfs(G, nodes=False, edges=True, node_geometry=False, fill_edge_geometry=True)
    edges = edges[['geometry', 'length', 'osmid', 'u','v', 'highway','key' 'oneway', 'maxspeed','name']]
    
    edges = edges.rename(columns = {'u':'old_u'})
    edges = edges.rename(columns = {'v':'old_v'})
    
    nodes = nodes.reset_index(drop=True)
    nodes['old_nodeID'] = nodes.osmid.astype('int64')
    nodes['nodeID'] = nodes.index.values.astype(int)

    edges = pd.merge(edges, nodes[['old_nodeID', 'nodeID']], how='left', left_on="old_u", right_on="old_nodeID")
    edges = edges.rename(columns = {'nodeID':'u'})
    edges = pd.merge(edges, nodes[['old_nodeID', 'nodeID']], how='left', left_on="old_v", right_on="old_nodeID")
    edges = edges.rename(columns = {'nodeID':'v'})
    
    edges = edges.reset_index(drop=True)
    edges['streetID'] = edges.index.values.astype(int)
    
    nodes = nodes[['nodeID','x','y','geometry']]
    nodes.gdf_name = 'Nodes_gdf' #for OSMNx
    edges = edges[['streetID','u','v','key','geometry', 'length', 'highway','oneway', 'maxspeed','name']]
    
    return(nodes, edges)

def get_fromSHP(directory, epsg, crs, simplify = False, area = None,
                roadType_field = None, direction_field = None, speed_field = None, name_field = None):
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

    streets_gdf = gpd.read_file(directory)
    streets_gdf = streets_gdf.to_crs(epsg=epsg)
    
    if (area != None):
        cn = streets_gdf.geometry.unary_union.centroid
        buffer = cn.buffer(area) 
        streets_gdf = streets_gdf[streets_gdf.geometry.within(buffer)]
        
    columns = [roadType_field, direction_field, speed_field, name_field]
    new_columns = ['highway','oneway', 'maxspeed','name']
    streets_gdf['from'] = "NaN"
    streets_gdf['to'] = "NaN"
    
    for n, i in enumerate(columns):
        if (i is not None): streets_gdf[new_columns[n]] = streets_gdf[i]
     
    standard_columns = ['geometry', 'from', 'to']
    streets_gdf = streets_gdf[standard_columns + 
                              [new_columns[n] for n, i in enumerate(columns) if i is not None]]
    
    index_geometry = streets_gdf.columns.get_loc("geometry")+1
    
    for row in streets_gdf.itertuples():
        line = []
        
        # to remove Z coordinates (assuming a planar graph)
        coord = list(row[index_geometry].coords)
        from_node = coord[0][0:2]
        to_node = coord[-1][0:2]
                
        for i in range(0, len(coord)):
            point = coord[i][0:2]
            line.append(point)

        t = LineString([coor for coor in line])
        
        streets_gdf.set_value(row[0],'geometry', t)
        streets_gdf.set_value(row[0], 'from', from_node)
        streets_gdf.set_value(row[0], 'to', to_node)
        
    streets_gdf = streets_gdf.loc[streets_gdf['from'] != streets_gdf['to']]
    unique_nodes_tmp = list(streets_gdf['to'].unique()) + list(streets_gdf['from'].unique())
    unique_nodes = list(set(unique_nodes_tmp))
    
    if (simplify == True):
        for i in unique_nodes:
            tmp = streets_gdf[(streets_gdf['from'] == i) | (streets_gdf['to'] == i)]
            
            if(len(tmp) != 2): continue
            
            else: 
        
                index_first = tmp.iloc[0].name
                index_second = tmp.iloc[1].name

                if (tmp.iloc[0]['from'] == tmp.iloc[1]['from']):
                    streets_gdf.loc[index_first]['from'] = tmp.iloc[1]['to']
                elif (tmp.iloc[0]['from'] == tmp.iloc[1]['to']):
                    streets_gdf.loc[index_first]['from'] = tmp.iloc[1]['from']
                elif (tmp.iloc[0]['to'] == tmp.iloc[1]['from']):
                    streets_gdf.loc[index_first]['to'] = tmp.iloc[1]['to'] 
                else: #(tmp.iloc[0]['to'] == tmp.iloc[1]['to'])
                    streets_gdf.loc[index_first]['to'] = tmp.iloc[1]['from']

                multi_line = MultiLineString([tmp.iloc[0]['geometry'], tmp.iloc[1]['geometry']])
                merged_line = linemerge(multi_line)
                streets_gdf.loc[index_first]['geometry'] = merged_line
                streets_gdf.drop([index_second], axis=0)
       
    # getting unique nodes again
    streets_gdf = streets_gdf.loc[streets_gdf['from'] != streets_gdf['to']]
    unique_nodes_tmp = list(streets_gdf['to'].unique()) + list(streets_gdf['from'].unique())
    unique_nodes = list(set(unique_nodes_tmp))
    
    # assigning indexes
    streets_gdf.reset_index(inplace=True, drop=True)
    streets_gdf['streetID'] = streets_gdf.index.values.astype(int) 
    
    #preparing nodes geodataframe
    nodes_data = pd.DataFrame.from_records(unique_nodes, columns=['x', 'y']).astype('float')
    geometry = [Point(xy) for xy in zip(nodes_data.x, nodes_data.y)]
    nodes = gpd.GeoDataFrame(nodes_data, crs=crs, geometry=geometry)
    nodes.reset_index(drop=True, inplace = True)
    nodes['nodeID'] = nodes.index.values.astype(int)
    nodes['coordinates'] = list(zip(nodes.x, nodes.y))

    edges_tmp = pd.merge(streets_gdf, nodes[['nodeID','coordinates']], how='left', left_on="from", right_on="coordinates")
    edges_tmp.drop(edges_tmp[['coordinates']], axis = 1, inplace = True)
    edges_tmp.rename(columns = {'nodeID':'u'}, inplace = True)
    
    edges = pd.merge(edges_tmp, nodes[['nodeID','coordinates']], how='left', left_on="to", right_on="coordinates")
    edges = edges.drop(edges[['coordinates', 'from', 'to']], axis = 1)
    edges = edges.rename(columns = {'nodeID':'v'})
    edges['key']=0 #for OSMNx
    edges['length'] = gpd.GeoSeries(edges['geometry'].length)
    nodes.drop(['coordinates'], axis = 1, inplace = True)
    nodes.gdf_name = 'Nodes_gdf' #for OSMNx
        
    return(nodes, edges)
	
def graph_fromGDF(nodes, edges):
    """
    It creates from due geopandas dataframes (street junctions and street segments) a NetworkX graph, passing by a OSMnx function.
    The length of the segments and their ID is stored in the edges attributes, as well as the nodeID.
    
    Parameters
    ----------
    nodes: geopandas dataframe
    edges: geopandas dataframe    
    """
    G = ox.gdfs_to_graph(nodes, edges)
    t = G.nodes()
    pos = {}

    for l, item in enumerate(t): pos[item] = (t[l]['x'],t[l]['y'], t[item]['nodeID'])

    Ng = nx.Graph() #Empty graph
    Ng = Ng.to_undirected()
    Ng.add_nodes_from(pos.keys()) #Add nodes preserving coordinates

    for i, item in enumerate(Ng.nodes()):
        Ng.node[item]['x']=pos[item][0]
        Ng.node[item]['y']=pos[item][1]
        Ng.node[item]['nodeID']=pos[item][2]

    for i, item in enumerate(G.edges()):
        Ng.add_edge(item[0], item[1])
        Ng[item[0]][item[1]]['length']=G[item[0]][item[1]][0]['length']
        Ng[item[0]][item[1]]['streetID']=G[item[0]][item[1]][0]['streetID']
        
    return(Ng)

# centroids for dual analysis

def dual_gdf(nodes, edges, crs):
    """
    It creates two dataframes that are supposed to generate the dual graph of a street network. The nodes_dual gdf contains edges 
    centroids, the edges_dual gdf contains instead links between the street segment centroids. Those dual edges link real street segment 
    that share a junction.
    The centroids are stored with the original edge streetID, while the dual edges are associated with several attributes computed on the 
    original street segments (distance between centroids, deflection angle).
    
    Parameters
    ----------
    nodes: geopandas dataframe
    edges: geopandas dataframe  
    crs: dictionary
    """
    
    centroids_gdf = edges.copy()
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
    
   
    for row in nodes_dual.itertuples():
        
        streetID = row[1] #streetID of the relative segment
        length = row[3]

        for i in list(row[2]): #intersecting segments

            # i is the streetID
            length_i =  centroids_gdf['length'][centroids_gdf.streetID == i][i]
            distance = (length+length_i)/2
        
            # adding a row with u-v, key fixed as 0, Linestring geometry 
            # from the first centroid to the centroid intersecting segment 
            ls = LineString([row[4], nodes_dual.loc[streetID]['geometry']])
            
            edges_dual.loc[-1] = [streetID, i, 0, ls, distance] 
            edges_dual.index = edges_dual.index + 1
            
    edges_dual = edges_dual.sort_index(axis=0)
    geometry = edges_dual['geometry']
    edges_dual = gpd.GeoDataFrame(edges_dual[['u','v', 'key', 'length']], crs=crs, geometry=geometry)

    #computing deflection angle
    
    for row in edges_dual.itertuples():

        #retrieveing original lines from/to
        from_node = nodes.loc[edges['u'].loc[edges.streetID == row[1]]].index.tolist()[0]
        to_node = nodes.loc[edges['v'].loc[edges.streetID == row[1]]].index.tolist()[0]
                                            
        from_node2 =  nodes.loc[edges['u'].loc[edges.streetID == row[2]]].index.tolist()[0]           
        to_node2 =  nodes.loc[edges['v'].loc[edges.streetID == row[2]]].index.tolist()[0]
    
        if ((from_node == from_node2) & (to_node == to_node2) | (from_node == to_node2) & (to_node == from_node2)):
            deflection = 0
            deflection_rad=0
    
        else:
         
            try:  
        
                x_f = float("{0:.10f}".format(nodes.loc[from_node]['x']))
                y_f = float("{0:.10f}".format(nodes.loc[from_node]['y']))
                x_t = float("{0:.10f}".format(nodes.loc[to_node]['x']))
                y_t = float("{0:.10f}".format(nodes.loc[to_node]['y']))
            
                x_f2 = float("{0:.10f}".format(nodes.loc[from_node2]['x']))
                y_f2 = float("{0:.10f}".format(nodes.loc[from_node2]['y']))    
                x_t2 = float("{0:.10f}".format(nodes.loc[to_node2]['x']))
                y_t2 = float("{0:.10f}".format(nodes.loc[to_node2]['y']))          
                             
                if (to_node == to_node2):
                    lineA = ((x_f, y_f),(x_t,y_t))
                    lineB = ((x_t2, y_t2),(x_f2, y_f2))
    
                elif (to_node == from_node2):
                    lineA = ((x_f, y_f),(x_t,y_t))
                    lineB = ((x_f2, y_f2),(x_t2, y_t2))

                elif (from_node == from_node2):
                    lineA = ((x_t, y_t),(x_f,y_f))
                    lineB = ((x_f2, y_f2),(x_t2, y_t2))

                else: #(from_node == to_node2)
                    lineA = ((x_t, y_t),(x_f,y_f))
                    lineB = ((x_t2, y_t2),(x_f2, y_f2))
        
                deflection = ang(lineA, lineB)
                deflection_rad = ang_rad(lineA, lineB)
            
            except:
                deflection = 0
                deflection_rad = 0
    
        edges_dual.set_value(row[0],'deg', deflection)
        edges_dual.set_value(row[0],'rad', deflection_rad)
        
    return(nodes_dual, edges_dual)

def get_dual_graph(nodes_dual, edges_dual):
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
        DG[item[0]][item[1]]['length'] = Gr[item[0]][item[1]][0]['length']
        DG[item[0]][item[1]]['deg'] = Gr[item[0]][item[1]][0]['deg']
        DG[item[0]][item[1]]['rad'] = Gr[item[0]][item[1]][0]['rad']
        
    return(DG)

	
# natural roads extraction: the function has to be called twice for each segment, to check in both directions

def natural_roads(streetID, natural_id, direction, roads_gdf, nodes_gdf): 
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
    roads_gdf: geopandas dataframe
    nodes_gdf: geopandas dataframe
    """
        
    angles = {}
    directions_dict = {}
        
    to_node = df_roads.loc[streetID]['u']
    from_node = df_roads.loc[streetID]['v']
        
    #assuming nodeID = ID dataframe
    x_t = float("{0:.10f}".format(nodes_gdf.loc[to_node]['x']))
    y_t = float("{0:.10f}".format(nodes_gdf.loc[to_node]['y']))
    x_f = float("{0:.10f}".format(nodes_gdf.loc[from_node]['x']))
    y_f = float("{0:.10f}".format(nodes_gdf.loc[from_node]['y']))
      
    #continue from the to_node     
    if (direction == "to"): intersecting = roads_gdf[(roads_gdf['u']== to_node) | (roads_gdf['v']== to_node)]
    
    #continue from the from_node
    else: intersecting = roads_gdf[(roads_gdf['u'] == from_node) | (roads_gdf['v'] == from_node)]

    if (len(intersecting) == 0): return
    
    for row_F in intersecting.itertuples():
        if ((streetID == row_F[0]) | (row_F[-1] != 'NA')): continue
        
        to_node_F = row_F[intersecting.columns.get_loc("u")]
        from_node_F = row_F[intersecting.columns.get_loc("v")]

        x_tf = float("{0:.10f}".format(nodes_gdf.loc[to_node_F]['x']))
        y_tf = float("{0:.10f}".format(nodes_gdf.loc[to_node_F]['y']))
        x_ff = float("{0:.10f}".format(nodes_gdf.loc[from_node_F]['x']))
        y_ff = float("{0:.10f}".format(nodes_gdf.loc[from_node_F]['y']))

        if (to_node == to_node_F):
            lineA = ((x_f, y_f),(x_t,y_t))
            lineB = ((x_tf, y_tf),(x_ff, y_ff))
            towards = "fr"
        elif (to_node == from_node_F):
            lineA = ((x_f, y_f),(x_t,y_t))
            lineB = ((x_ff, y_ff),(x_tf, y_tf))
            towards = "to"
        elif (from_node == from_node_F):
            lineA = ((x_t, y_t),(x_f,y_f))
            lineB = ((x_ff, y_ff),(x_tf, y_tf))
            towards = "to"
        else: #(from_node == to_node_F)
            lineA = ((x_t, y_t),(x_f,y_f))
            lineB = ((x_tf, y_tf),(x_ff, y_ff))
            towards = "fr"

        deflection = ang(lineA, lineB)

        if (deflection >= 45): continue
        else:
            angles[row_F[0]] = deflection
            directions_dict[row_F[0]] = towards

    if (len(angles) == 0): #no natural continuations
        roads_gdf.set_value(streetID, 'natural_id', natural_id)
        return
    else:
        angles_sorted = sorted(angles, key = angles.get)
        matchID = angles_sorted[0]
        roads_gdf.set_value(id_road, 'natural_id', natural_id)
        natural_roads(matchID, natural_id, directions_dict[matchID], roads_gdf, nodes_gdf)
    
    
    
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
	

# landmarks

def select_buildings(city_buildings, area_to_clip, base = None):
    buildings = city_buildings[city_buildings.geometry.within(area_to_clip.geometry.loc[0])]
    obstructions = city_buildings[city_buildings.geometry.within(area_to_clip.geometry.loc[0].buffer(200))]
    
    buildings["area"] = buildings['geometry'].area
    if (base == None): buildings["base"] = 0
    else: buildings["base"] = buildings[base]
        
    buildings = buildings[buildings['area']>199]
    buildings = buildings[['height', 'base','geometry', 'area']]
    buildings['buildingID'] = buildings.index.values.astype(int)
    
    return(buildings, obstructions)


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
        #precise matches not necessary here
    
        if(len(possible_matches_index) < 1): continue 
        polygon = MultiPolygon([pol.buffer(11) for pol in possible_matches['geometry']])
        area_max = (b200.difference(cascaded_union(polygon))).area
        buildings_gdf.set_value(row[0], 'prom', area_max) #prominence
        
        possible_neigh_index = list(sindex.intersection(b150.bounds))
        possible_neigh = obstructions.iloc[possible_neigh_index]
        precise_neigh = possible_neigh[possible_neigh.intersects(b150)]
        buildings_gdf.set_value(row[0], 'neigh', len(precise_neigh)) #neighbours
        
        dist = g.distance(street_network)
        buildings_gdf.set_value(row[0], 'road', dist) #distance road

        
    buildings_gdf['ext'] = buildings_gdf.area*(buildings_gdf.length/buildings_gdf.width) #extension
    buildings_gdf['fac'] = buildings_gdf['height']*(buildings_gdf.width) #facade area
    buildings_gdf.drop(['width','length'], axis=1, inplace = True)
    
    return buildings_gdf


def visibility(buildings_gdf, sight_lines):

    avg = (sight_lines[['buildingID', 'Shape_Leng']].groupby(['buildingID'], 
                                            as_index=False)['Shape_Leng'].mean())
    
    avg.rename(columns={'Shape_Leng': 'mean_length'}, inplace=True)
    
    count =(sight_lines[['buildingID', 'Shape_Leng']].groupby(['buildingID'], 
                                            as_index=False)['Shape_Leng'].count())
    count.rename(columns={'Shape_Leng': 'n_lines'}, inplace=True)

    tmp = sight_lines.set_index('buildingID')
    distant = tmp.groupby('buildingID').agg(lambda copy: copy.values[copy['Shape_Leng'].values.argmax()])
    distant = distant.reset_index()

    visibility = pd.merge(distant, avg, left_on='buildingID', right_on='buildingID')
    visibility.drop(['DIST_ALONG','Visibility', 'geometry'], axis=1, inplace=True)
    visibility.rename(columns = {'Shape_Leng':'distance','mean_length':'mean_distance'}, inplace=True) 
    
    tmp = pd.merge(buildings_gdf, visibility[['buildingID','distance', 'mean_distance']], on='buildingID', how='left')  
    tmp['distance'].fillna((tmp['distance'].min()), inplace=True)
    tmp['mean_distance'].fillna((tmp['mean_distance'].min()), inplace=True)
    tmp.rename(columns={'distance': 'vis', 'mean_distance': 'mean_vis'}, inplace=True)
    
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
        precise_matches = buildings_gdf[buildings_gdf.intersects(b)]

        neigh = precise_matches.groupby(['land_use'], as_index=True)['nr'].sum()
        Nj = neigh.loc[use]
        #Pj = Nj/N

        Pj = 1-(Nj/precise_matches['nr'].sum())
        buildings_gdf.set_value(row[0], 'prag', Pj) #pragmatic meaning
        
    return buildings_gdf
        
        
def compute_scores(buildings_gdf):

    col = ['area', 'ext', 'fac', 'height', 'prag', 'prom','cult', 'vis']
    for i in col: scaling_columnDF(buildings_gdf, i)

    col = ['neigh', 'road']
    for i in col: scaling_columnDF(buildings_gdf, i, inverse = True)     

    buildings_gdf['vScore']= (buildings_gdf['fac_sc']*30 + buildings_gdf['height_sc']*20 + buildings_gdf['vis_sc']*50)/100
    buildings_gdf['sScore']= (buildings_gdf['ext_sc']*30 + buildings_gdf['neigh_sc']*20 + buildings_gdf['prom_sc']*30 
                              + buildings_gdf['road_sc']*20)/100

    col = ['vScore', 'sScore']
    for i in col: scaling_columnDF(buildings_gdf, i)

    buildings_gdf['gScore']=(buildings_gdf['vScore_sc']*50 + buildings_gdf['sScore_sc']*30 + buildings_gdf['cult_sc']*10 
                       + buildings_gdf['prag_sc']*10)/100

    scaling_columnDF(buildings_gdf, 'gScore')
    
    return buildings_gdf


def local_scores(landmarks_gdf, point, buffer_extension):
    
    LL = landmarks_gdf[landmarks_gdf.geometry.within(point.buffer(buffer_extension))]
    
    col = ['area', 'ext', 'fac', 'height', 'prag', 'prom','cult', 'vis']
    for i in col: scaling_columnDF(LL, i)

    col = ['neigh', 'road']
    for i in col: scaling_columnDF(LL, i, inverse = True)
        
    LL['vScore']= (LL['fac_sc']*30 + LL['height_sc']*20 + LL['vis_sc']*50)/100
    LL['sScore']= (LL['ext_sc']*30 + LL['neigh_sc']*20 + LL['prom_sc']*30 + LL['road_sc']*20)/100

    col = ['vScore', 'sScore']
    for i in col: scaling_columnDF(LL, i)

    LL['lScore']=(LL['vScore_sc']*20 + LL['sScore_sc']*30 + LL['cult_sc']*10 + LL['prag_sc']*40)/100
    scaling_columnDF(LL, 'lScore')
    
    return LL


def decision_score(nodes_gdf, landmarks_gdf):
    
        
    # assigning a decision point score per each node
    # each node get the score of the best local landmark

    spatial_index = landmarks_gdf.sindex
    nodes_gdf['DE'] = 0.0
    index_geometry = nodes_gdf.columns.get_loc("geometry")+1

    for row in nodes_gdf.itertuples():
        g = row[index_geometry] #geometry
        b = g.buffer(50)    
        local_landmarks = local_scores(landmarks_gdf, g, 800)
        
        spatial_index = local_landmarks.sindex
        possible_matches_index = list(spatial_index.intersection(b.bounds))
        possible_matches = local_landmarks.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(b)]
        
        if (len(precise_matches)==0): continue
    
        precise_matches = precise_matches.sort_values(by='lScore', ascending=False).reset_index()
        score = precise_matches['lScore'].loc[0]
        nodes_gdf.set_value(row[0], 'DE', score)
        
    return nodes_gdf


def distant_score(nodes_gdf, landmarks_gdf, sight_lines):
    
    # keeping relevant landmarks and the sight lines that point at them

    relevant = landmarks_gdf[landmarks_gdf.gScore_sc > 0.5]
    index_relevant = landmarks_gdf['buildingID'].values.astype(int) 
    sightLine_to_relevant = sight_lines[sight_lines['buildingID'].isin(index_relevant)]

    # per each node, the sight lines to the relevant landmarks are extracted. 
    # The visibility score of each node is the sum of the score of the visible landmarks visible from it, regardless the direction

    nodes_gdf['DD'] = 0.0
    nodes_gdf['visible_landmarks'] = 'NA' 

    index_nodeID = nodes_gdf.columns.get_loc("nodeID")+1
    
    for row in nodes_gdf.itertuples():
        sight_node = sightLine_to_relevant[sightLine_to_relevant['nodeID'] == row[index_nodeID]] 
        index_relevant_fromNode = list(sight_node['buildingID'].values.astype(int))
        relevant_fromNode = relevant[relevant['buildingID'].isin(index_relevant_fromNode)] 
        
        score = relevant_fromNode['vis_sc'].sum()
        nodes_gdf.set_value(row[0], 'DD', score)   
        nodes_gdf.set_value(row[0], 'visible_landmarks', list(relevant_fromNode['buildingID'].values.astype(int)))
        
    return nodes_gdf 
    
