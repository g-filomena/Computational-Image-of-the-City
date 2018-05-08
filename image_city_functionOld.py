import osmnx as ox, networkx as nx, matplotlib.cm as cm, pandas as pd, numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import functools
import community
import math

from scipy import sparse
from scipy.sparse import linalg

from shapely.geometry import Point, LineString, Polygon, MultiPolygon, mapping
from math import sqrt
import pandas as pd
from shapely.ops import cascaded_union
pd.set_option('precision', 10)


# other functions
	
def scaling_columnDF(df, i):
    df[i+'_sc'] = (df[i]-df[i].min())/(df[i].max()-df[i].min())
	
def nodes_dict(G):
    nodes_list = G.nodes()
    nodes_dict = {}

    for i, item in enumerate(nodes_list):
        cod = item
        x = nodes_list[item]['x']
        y = nodes_list[item]['y']
        nodes_dict[cod] = (x,y)
    
    return nodes_dict

def dict_to_df(list_dict, list_col):
    
    df = pd.DataFrame(list_dict).T
    df.columns = ['d{}'.format(i) for i, col in enumerate(df, 1)]
    df.columns = list_col
    
    return(df)

# math functions

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
    G = ox.graph_from_place(place, network_type='all', simplify=True)
    G = ox.project_graph(G)
    
    for i, item in enumerate(G.edges()):
        if isinstance(G[item[0]][item[1]][0]['osmid'], (list,)):
            G[item[0]][item[1]][0]['osmid'] = G[item[0]][item[1]][0]['osmid'][0]
            
            
    nodes = ox.graph_to_gdfs(G, nodes=True, edges=False, node_geometry=True, fill_edge_geometry=False)
    nodes = nodes.drop(nodes[['highway', 'ref', 'lat', 'lon']], axis=1)
    edges = ox.graph_to_gdfs(G, nodes=False, edges=True, node_geometry=False, fill_edge_geometry=True)
    edges = edges[['geometry', 'length', 'osmid', 'u','v']]
    
    edges = edges.rename(columns = {'u':'old_u'})
    edges = edges.rename(columns = {'v':'old_v'})
    
    nodes = nodes.reset_index(drop=True)
    print(nodes.head())
    nodes['old_nodeID'] = nodes.osmid.astype('int64')
    print(nodes.head())
    nodes['nodeID'] = nodes.index.values.astype(int)

    edges = pd.merge(edges, nodes[['old_nodeID', 'nodeID']], how='left', left_on="old_u", right_on="old_nodeID")
    edges = edges.rename(columns = {'nodeID':'u'})
    edges = pd.merge(edges, nodes[['old_nodeID', 'nodeID']], how='left', left_on="old_v", right_on="old_nodeID")
    edges = edges.rename(columns = {'nodeID':'v'})
    
    edges = edges.reset_index(drop=True)
    edges['streetID'] = edges.index.values.astype(int)
    
    nodes = nodes[['nodeID','x','y','geometry']]
    nodes.gdf_name = 'Nodes_gdf' #for OSMNx
    edges = edges[['streetID','u','v','key','geometry']]
    
    return(nodes, edges)

def get_fromSHP(directory, epsg):
    
    #try reading street network

    streets_gdf = gpd.read_file(directory)
    streets_gdf = streets_gdf.to_crs(epsg=epsg)
    
    streets_gdf['from'] = "NaN"
    streets_gdf['to'] = "NaN"
    streets_gdf = streets_gdf[['geometry', 'from', 'to']]
    
    for index, row in streets_gdf.iterrows():
        line = []
        
        # to remove Z coordinates (assuming a planar graph)
        coord = list(row['geometry'].coords)
        from_node = coord[0][0:2]
        to_node = coord[-1][0:2]
    
        for i in range(0, len(coord)):
            point = coord[i][0:2]
            line.append(point)

        t = LineString([coor for coor in line])
        
        streets_gdf.set_value(index,'geometry', t)
        streets_gdf.set_value(index, 'from', from_node)
        streets_gdf.set_value(index, 'to', to_node)
        
    #removing pseudo-lines and assigning Index
    streets_gdf = streets_gdf.loc[streets_gdf['from'] != streets_gdf['to']]
    streets_gdf.reset_index(inplace=True, drop=True)
    streets_gdf['streetID'] = streets_gdf.index.values.astype(int) 
    

    #getting unique nodes
    unique_nodes_tmp = list(streets_gdf['to'].unique()) + list(streets_gdf['from'].unique())
    unique_nodes = list(set(unique_nodes_tmp))
    
    #preparing nodes geodataframe
    nodes_data = pd.DataFrame.from_records(unique_nodes, columns=['x', 'y']).astype('float')
    geometry = [Point(xy) for xy in zip(nodes_data.x, nodes_data.y)]
    nodes = gpd.GeoDataFrame(nodes_data, crs=crs, geometry=geometry)
    nodes = nodes.reset_index(drop=True)
    nodes['nodeID'] = nodes.index.values.astype(int)
    nodes['coordinates'] = list(zip(nodes.x, nodes.y))
    nodes.gdf_name = 'Nodes_gdf' #for OSMNx
    
    edges_tmp = pd.merge(streets_gdf, nodes[['nodeID','coordinates']], how='left', left_on="from", right_on="coordinates")
    edges_tmp = edges_tmp.drop(edges_tmp[['coordinates']], axis = 1)
    edges_tmp = edges_tmp.rename(columns = {'nodeID':'u'})
    
    edges = pd.merge(edges_tmp, nodes[['nodeID','coordinates']], how='left', left_on="to", right_on="coordinates")
    edges = edges.drop(edges[['coordinates', 'from', 'to']], axis = 1)
    edges = edges.rename(columns = {'nodeID':'v'})
    edges['key']=0 #for OSMNx
    edges['length'] = gpd.GeoSeries(edges['geometry'].length)
        
    return(nodes, edges)
	
def graph_fromGDF(nodes, edges):
    G = ox.gdfs_to_graph(nodes, edges)
    t = G.nodes()
    pos = {}
    idx = {}

    for l, item in enumerate(t): 
        pos[item] = (t[l]['x'],t[l]['y'])
        idx[item] = t[item]['nodeID']

    Ng = nx.Graph() #Empty graph
    Ng = Ng.to_undirected()
    Ng.add_nodes_from(pos.keys()) #Add nodes preserving coordinates

    for i, item in enumerate(Ng.nodes()):
        Ng.node[item]['x']=pos[item][0]
        Ng.node[item]['y']=pos[item][1]
        Ng.node[item]['nodeID']=idx[item]

    for i, item in enumerate(G.edges()):
        Ng.add_edge(item[0], item[1])
        Ng[item[0]][item[1]]['length']=G[item[0]][item[1]][0]['length']
        Ng[item[0]][item[1]]['streetID']=G[item[0]][item[1]][0]['streetID']
        
    return(Ng)

# centroids for dual analysis

def dual_gdf(edges, nodes, crs):
    
    centroids_gdf = edges.copy()
    centroids_gdf['centroid'] = gpd.GeoSeries(centroids_gdf['geometry'].centroid)
    centroids_gdf['intersecting'] = centroids_gdf['intersecting'].astype(object)
    
    index_u = centroids_gdf.columns.get_loc("u")
    index_v = centroids_gdf.columns.get_loc("v")
    index_streetID = centroids_gdf.columns.get_loc("streetID")
         
    #find_intersecting segments
    processed = []    
    for c in centroids_gdf.itertuples():
        
        intersections = []
        from_node = e[index_u]
        to_node = e[index_v]
    
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
            length_i =  centroids_gdf['length'][centroids_gdf.streetID == i][streetID]
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


	
# centrality functions

def straightness_centrality(G, weight, normalized = True):

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

            straightness_centrality[n]= straightness
               
            if normalized:
                straightness_centrality[n] = straightness * (1.0/(len(G)-1.0) )

        else:
            straightness_centrality[n]=0.0

    return straightness_centrality

def reach_centrality(G, weight, radius, normalized=True):
  
    path_length = functools.partial(nx.single_source_dijkstra_path_length, weight=weight)

    nodes = G.nodes()
    reach_centrality = {}
    coord_nodes = nodes_dict(G)

    for n in nodes:
        reach = 0
        sp = path_length(G,n)
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

def local_betweenness_centrality(G, w, radius, distance=True):

    if distance is True:
        weight=w
        
    path_length = functools.partial(nx.single_source_dijkstra_path_length, weight=weight)

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
	

