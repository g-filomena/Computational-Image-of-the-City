# def is_t_junciont(nodes_gdf, edges_gdf, nodeID):
#     tmp = edges_gdf[edges_gdf.u == nodeID | edges_gdf.v == nodeID].copy()
#     if len(tmp) != 3: return(False)
    
#     if ((isparallel(tmp.iloc[0].geometry, tmp.iloc[1].geometry) == True) 
#         & (isperpendicular(tmp.iloc[0].geometry, tmp.iloc[2].geometry) == True))
         
#     elif ((isparallel(tmp.iloc[0].geometry, tmp.iloc[2].geometry) == True) & 
#          (isperpendicular(tmp.iloc[0].geometry, tmp.iloc[1].geometry) == True)):
#     elif ((isparallel(tmp.iloc[1].geometry, tmp.iloc[2].geometry) == True) &
#          (isperpendicular(tmp.iloc[0].geometry, tmp.iloc[1].geometry) == True)):

# math functions for angle computations
# from Abhinav Ramakrishnan answer in https://stackoverflow.com/a/28261304/7375309


def euclidean_distance(xs, ys, xt, yt):
    """ xs stands for x source and xt for x target """
    return sqrt((xs - xt)**2 + (ys - yt)**2)


def isparallel(geolineA, geolineB, hard = False):
    
    angle = uf.ang_geoline(geolineA, geolineB, degree = True)
    if ((angle <= 20) | (angle >= 160)): return True
        
    line_coordsA = list(geolineA.line_coords)
    line_coordsB = list(geolineB.line_coords)
    
    if ((len(line_coordsA) == 2) | (len(line_coordsB) == 2)): return False
    
    spA = Point(line_coordsA[0])
    spB = Point(line_coordsB[0])
    epB = Point(line_coordsB[-1])
    if spA.distance(spB) > spA.distance(epB): line_coordsB.reverse()
    
    if hard == False:
        # remove 1st coordinates (A,B)
        geolineA = LineString([coor for coor in line_coordsA[1:]])
        geolineB = LineString([coor for coor in line_coordsB[1:]])
        angle = uf.ang_geoline(geolineA, geolineB, degree = True)
        if ((angle <= 20) | (angle >= 160)): return True
        
        # remove 1st (A) and last (B)
        geolineB = LineString([coor for coor in line_coordsB[:-1]])
        angle = uf.ang_geoline(geolineA, geolineB, degree = True)
        if ((angle <= 20) | (angle >= 160)): return True
        
        # remove last (A) and 1st (B)
        geolineA = LineString([coor for coor in line_coordsA[:-1]])
        geolineB = LineString([coor for coor in line_coordsB[1:]])
        angle = uf.ang_geoline(geolineA, geolineB, degree = True)
        if ((angle <= 20) | (angle >= 160)): return True
        
        # remove last coordinates (A, B)
        geolineA = LineString([coor for coor in line_coordsA[:-1]])
        geolineB = LineString([coor for coor in line_coordsB[:-1]])
        angle = uf.ang_geoline(geolineA, geolineB, degree = True)
        if ((angle <= 20) | (angle >= 160)): return True
        
        if ((len(line_coordsA) == 3) | (len(line_coordsB) == 3)): return False
        geolineA = LineString([coor for coor in line_coordsA[1:-1]])
        geolineB = LineString([coor for coor in line_coordsB[1:-1]])

        angle = uf.ang_geoline(geolineA, geolineB, degree = True)
        if ((angle <= 20) | (angle >= 160)): return True
    return False






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
        coming_from = v
        possible_matches = edges_gdf[(edges_gdf.u == v)].copy()
    else: 
        u = edges_gdf.loc[index_line]['u']
        coming_from = u
        possible_matches = edges_gdf[(edges_gdf.v == u)].copy()
    
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
    
    found = False
    while (found is False):
        if len(possible_matches) == 0: 
            found = True
            return(None, None, None, None, None)
        
        for p in possible_matches.itertuples():
            if isparallel(line, p[index_geometry], hard = True) is False:
                possible_matches.drop(p[0], axis = 0, inplace = True)
                continue
            else:
                uCP, vCP = p[index_u], p[index_v]
                
                if search_direction == 'v': 
                    cluster = nodes_gdf.loc[vCP].cluster
                    coming_from = vCP
                    if vCP in nodes_encountered: 
                        possible_matches = possible_matches[0:0]
                        break
                else:  
                    cluster = nodes_gdf.loc[uCP].cluster
                    if uCP in nodes_encountered: 
                        possible_matches = possible_matches[0:0]
                        break
                                    
                if cluster == 'NA':
                    lines_traversed.append(p[0])
                    
                    if search_direction == 'v': 
                        possible_matches = edges_gdf[(edges_gdf.u == vCP)].copy()
                        nodes_encountered.append(uCP) 
                        line_coords = line_coords + list(p[index_geometry].coords)
                    else:
                        possible_matches = edges_gdf[(edges_gdf.v == uCP)].copy()
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
        geo_lineA = LineString([coor for coor in new_line_A])
        geo_lineB = LineString([coor for coor in new_line_B])

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
        geo_lineA = LineString([coor for coor in new_line_A])
        geo_lineB = LineString([coor for coor in new_line_B])
    
    return((lineA, lineB), np)

def interpolate(u, center_line, last_node, nodes_to_interpolate, nodes_gdf, edges_gdf, index_line):
    print('interpolating1', index_line)
    geo_line = center_line    
    new_index = index_line
    
    for counter, node in enumerate(nodes_to_interpolate):

        result, np = split_line_interpolation(node, geo_line, nodes_gdf, edges_gdf)
              
        #first part of the segment, adjusting node coordinates
        nodes_gdf.set_value(node, 'x', np.coords[0][0])
        nodes_gdf.set_value(node, 'y', np.coords[0][1])
        nodes_gdf.set_value(node, 'geometry', np)
        
        if counter == 0: edges_gdf.set_value(new_index, 'u', u)
            
        edges_gdf.set_value(new_index, 'geometry', result[0])
        edges_gdf.set_value(new_index, 'v', node)
                
        # second part of the segment
        new_index = max(edges_gdf.index)+1
        print(new_index, 'new_index') 
        edges_gdf.loc[new_index] = edges_gdf.loc[index_line]
        edges_gdf.set_value(new_index, 'geometry', result[1])
        edges_gdf.set_value(new_index, 'u', node)
        edges_gdf.set_value(new_index, 'v', last_node)
        edges_gdf.set_value(new_index, 'streetID', new_index) 
        geo_line = result[1]
            
    return(nodes_gdf, edges_gdf)

def interpolate_multi(u, center_line, last_node, list_nodes, nodes_gdf, edges_gdf, index_line):
    print('interpolating', index_line)
    geo_line = center_line  
    new_index = index_line
    distances = {}
    
    for node in list_nodes:
        distance = nodes_gdf.loc[node]['geometry'].distance(Point(center_line.coords[0]))
        distances[node] = distance

    distance_sorted = sorted(distances.items(), key=lambda kv: kv[1])               
    for counter, node in enumerate(distance_sorted):
        
        node = distance_sorted[counter][0]
        result, np = split_line_interpolation(node, geo_line, nodes_gdf, edges_gdf)
        
        #first part of the segment, adjusting node coordinates
        nodes_gdf.set_value(node, 'x', np.coords[0][0])
        nodes_gdf.set_value(node, 'y', np.coords[0][1])
        nodes_gdf.set_value(node, 'geometry', np)

        if counter == 0: edges_gdf.set_value(new_index, 'u', u)
        print(new_index, 'new_index')   
        edges_gdf.set_value(new_index, 'geometry', result[0])
        edges_gdf.set_value(new_index, 'v', node)
         
        # second part of the segment
        new_index = max(edges_gdf.index)+1
        
        edges_gdf.loc[new_index] = edges_gdf.loc[index_line]
        edges_gdf.set_value(new_index, 'geometry', result[1])
        edges_gdf.set_value(new_index, 'u', node)
        edges_gdf.set_value(new_index, 'v', last_node) 
        edges_gdf.set_value(new_index, 'streetID', new_index) 
        geo_line = result[1]
            
    return(nodes_gdf, edges_gdf)

def simplify_dual_lines(nodes_gdf, edges_gdf, junctions_gdf):
    
    edges_gdf = edges_gdf.copy()
    nodes_gdf = nodes_gdf.copy()
    list_cluster = junctions_gdf.index.values.tolist()  
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['cluster', 'nodeID']], how = 'left', left_on= "u", right_on = "nodeID")
    edges_gdf = edges_gdf.rename(columns = {'cluster':'cluster_u'})
    edges_gdf = pd.merge(edges_gdf, nodes_gdf[['cluster', 'nodeID']], how = 'left', left_on= "v", right_on = "nodeID")
    edges_gdf = edges_gdf.rename(columns = {'cluster':'cluster_v'})  
    edges_gdf.set_index('streetID', drop = False, append = False, inplace = True)
    del edges_gdf.index.name
    
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
                        line_coords = list(g[index_geometry].coords)
                        line_coords.reverse() 
                        new_geo_line = LineString([coor for coor in line_coords])
                        old_u = g[index_u]
                        old_cluster_u = g[index_cluster_u]
                        old_cluster_uR = g[index_cluster_uR]

                        group.set_value(g[0],'geometry', new_geo_line)
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
                    geo_line, geo_lineC = group.iloc[0]['geometry'], group.iloc[1]['geometry']
                    dr, drC = group.iloc[0]['direction'], group.iloc[1]['direction']
                    index_line, index_lineC  = group.iloc[0].name, group.iloc[1].name

                    if (c_v == c_vC) & (c_v != 'NA'):
                        destination = c_v
                        cl = center_line(u, v, uC, vC, geo_line, geo_lineC)
                        edges_gdf.drop(index_lineC, axis = 0, inplace = True)
                        if dr == 'u':
                            line_coords = list(cl.coords)
                            line_coords.reverse() 
                            cl = LineString([coor for coor in line_coords])
                        
                        edges_gdf.set_value(index_line, 'geometry', cl)
                        processed = processed + [index_line, index_lineC]
                        junctions_gdf.set_value(destination, 'keep', True)
                        break # next group

                    ######################################################## 
                    ## SUB-OPTION 2: only one reaches another cluster:

                    elif (c_v != 'NA') | (c_vC != 'NA'):

                        if c_v != 'NA': 
                            destination = c_v
                            found, geo_lineC, lines_t, nodes_en, vC = find_next_cluster(nodes_gdf,edges_gdf, index_lineC, drC)
                            last_node = vC
                        else: 
                            destination = c_vC
                            found, geo_line, lines_t, nodes_en, v = find_next_cluster(nodes_gdf,edges_gdf, index_line, dr)
                            last_node = v
                        
                        cl =  center_line_cluster(geo_line, geo_lineC, nodes_gdf, junctions_gdf, u, destination, one_cluster = True)
                        nodes_gdf, edges_gdf = interpolate(u, cl, last_node, nodes_en, nodes_gdf, edges_gdf, index_line)                                       
                        processed = processed + [index_line, index_lineC] + lines_t
                        lines_t.append(index_lineC)
                        edges_gdf.drop(lines_t, axis = 0, inplace = True, errors = 'ignore')
                        junctions_gdf.set_value(destination, 'keep', True)
                        break # next group                          

                    ####################################################### 
                    # SUB-OPTION 3: none reaches a cluster directly; comparing the first reached cluster
                    else:                    
                        dest, geo_line, lines_t, nodes_en, v = find_next_cluster(nodes_gdf, edges_gdf, index_line, dr)
                        destC, geo_lineC, lines_tC, nodes_enC, vC = find_next_cluster(nodes_gdf, edges_gdf, index_lineC, drC)    

                        # the center line is built in relation to the variable cluster as 'u', or from_node --> to_node
                        cl =  center_line_ocluster(geo_line, geo_lineC, nodes_gdf, junctions_gdf, u, dest, one_cluster = True)

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
                elif ((isparallel(row[index_geometry], rowC[index_geometry], hard = False) is True) |
                      (row[index_name] == rowC[index_name])): continue
                else: group.drop(rowC[0], axis = 0, inplace = True)

            # does the line considered in the loop reach a cluster? if not straight away, at some point?
            
            group['direction'] = 'v'
            # orientate everything from "u" to "v"
            for rowC in group.itertuples():
                if rowC[index_cluster_v] == cluster:
                    line_coords = list(rowC[index_geometry].coords)
                    line_coords.reverse() 
                    new_geo_line = LineString([coor for coor in line_coords])
                    old_u = rowC[index_u]
                    old_cluster_u = rowC[index_cluster_u]
                    old_cluster_uR = rowC[index_cluster_uR]

                    group.set_value(rowC[0],'geometry', new_geo_line)
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
                geo_line, geo_lineC = group.iloc[0]['geometry'], group.iloc[1]['geometry']
                index_line, index_lineC  = group.iloc[0].name, group.iloc[1].name
                
                ######################################################## 
                ## SUB-OPTION 1: they all reach another cluster:
                    
                if (c_v == c_vC) & (c_v != 'NA'):
                    destination = c_v
                    cl = center_line_cluster(geo_line, geo_lineC, nodes_gdf, junctions_gdf, cluster, destination)
                    edges_gdf.drop(index_lineC, axis = 0, inplace = True)
                    if dr == 'u':
                        line_coords = list(cl.coords)
                        line_coords.reverse() 
                        cl = LineString([coor for coor in line_coords])
                    
                    edges_gdf.set_value(index_line, 'geometry', cl)
                    processed = processed + [index_line, index_lineC]
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break # next group

                ######################################################## 
                ## SUB-OPTION 2: only one reaches another cluster:
                    
                elif (c_v != 'NA') | (c_vC != 'NA'):
                    
                    if c_v != 'NA': 
                        destination = c_v
                        found, geo_lineC, lines_t, nodes_en, vC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)
                        last_node = vC
                    else: 
                        destination = c_vC
                        found, geo_line, lines_t, nodes_en, v = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                        last_node = v
                    
                    cl = center_line_cluster(geo_line, geo_lineC, nodes_gdf, junctions_gdf, cluster, destination)
                    nodes_gdf, edges_gdf = interpolate(u, cl, last_node, nodes_en, nodes_gdf, edges_gdf, index_line)                                       
                    processed = processed + [index_line, index_lineC] + lines_t
                    lines_t.append(index_lineC)
                    edges_gdf.drop(lines_t, axis = 0, inplace = True, errors = 'ignore')
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break # next group                          
                    
                ####################################################### 
                # SUB-OPTION 3: none reaches a cluster directly; comparing the first reached cluster
                else: 
                    dest, geo_line, lines_t, nodes_en, v = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                    destC, geo_lineC, lines_tC, nodes_enC, vC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)    
                        
                    # the center line is built in relation to the variable cluster as 'u', or from_node --> to_node
                    cl = center_line_cluster(geo_line, geo_line, nodes_gdf, junctions_gdf, cluster, dest)
                        
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
                geo_line, geo_line, geo_lineCC = group.iloc[0]['geometry'], group.iloc[1]['geometry'], group.iloc[2]['geometry']
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
                    geo_line, geo_lineC = group.iloc[0]['geometry'], group.iloc[1]['geometry']
                    index_line, index_lineC  = group.iloc[0].name, group.iloc[1].name
                    
                    destination, line, lines_t, nodes_en, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)
                    cl = geo_line                                    
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
                        destC, geo_lineC, lines_tC, nodes_enC, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)
                        destCC, geo_lineCC, lines_tCC, nodes_enCC, last_nodeCC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineCC, drCC)
                        lines_t, nodes_en = [], []
                    elif (c_vC != 'NA'):
                        dest, geo_line, lines_t, nodes_en, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                        destCC, geo_lineCC, lines_tCC, nodes_enCC, last_nodeCC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineCC, drCC)
                        lines_tC, nodes_enC = [], []
                    else:
                        dest, geo_line, lines_t, nodes_en, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                        destC, geo_lineC, lines_tC, nodes_enC, last_nodeC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)                   
                        lines_tCC, nodes_enCC = [], []
                        
                        # exclude the 2 lines furter away                          
                    
                    max_dist = 0
                    dict_lines = {index_line: geo_line, index_lineC: geo_lineC, index_lineCC: geo_lineCC}
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

                    if central == index_line: cl = geo_line
                    elif central == index_lineC: cl = geo_line
                    else: cl = geo_line
                    
                    to_drop = secondary_lines + lines_t + lines_tC + lines_tCC
                    nodes_gdf, edges_gdf = interpolate_multi(u, cl, last_node, list_nodes, nodes_gdf, edges_gdf, index_line)     
                    edges_gdf.drop(to_drop, axis = 0, inplace = True, errors = 'ignore')
                    processed = processed + to_drop + [central]
                    junctions_gdf.set_value(cluster, 'keep', True)
                    break
                 
                ########################################################  
                ## SUB-OPTION 4: none reaches a cluster:
                else: 
                    dest, geo_line, lines_t, nodes_en, last_node = find_next_cluster(nodes_gdf, old_edges_gdf, index_line, dr)
                    destC, geo_lineC, lines_tC, nodes_enC, last_nodeC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineC, drC)
                    destCC, geo_lineCC, lines_tCC, nodes_enCC, last_nodeCC = find_next_cluster(nodes_gdf, old_edges_gdf, index_lineCC, drCC)  
                    # exclude the 2 lines furter away                          
                    
                    max_dist = 0
                    dict_lines = {index_line: geo_line, index_lineC: geo_lineC, index_lineCC: geo_lineCC}
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

                    if central == index_line: cl = geo_line
                    elif central == index_lineC: cl = geo_lineC
                    else: cl = geo_lineCC
                    
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

def reassign_edges(nodes_gdf, edges_gdf, junctions_gdf):
    
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
#     edges_gdf = edges_gdf[(edges_gdf.cluster_u != edges_gdf.cluster_v) | ((edges_gdf.cluster_u == 'NA') & (edges_gdf.cluster_v == 'NA'))]
    
    for row in edges_gdf.itertuples():

        line_coords = list(row[index_geometry].coords)

        u = nodes_gdf.loc[row[index_old_u]]["cluster"]
        v = nodes_gdf.loc[row[index_old_v]]["cluster"]
        old_u = row[index_old_u]
        old_v = row[index_old_v]
               
        if ((u != 'NA') & (v !='NA')):  # change starting and ending node in the list of coordinates for the line
                if (junctions_gdf.loc[u].keep == False) & (junctions_gdf.loc[v].keep == False): 
                    u = old_u
                    v = old_v
                    line_coords[0] = (nodes_gdf.loc[old_u]['x'], nodes_gdf.loc[old_u]['y'])
                    line_coords[-1] = (nodes_gdf.loc[old_v]['x'], nodes_gdf.loc[old_v]['y'])
                elif junctions_gdf.loc[v].keep == False:
                    v = old_v
                    line_coords[0] = (junctions_gdf.loc[u]['x'], junctions_gdf.loc[u]['y'])
                    line_coords[-1] = (nodes_gdf.loc[old_v]['x'], nodes_gdf.loc[old_v]['y'])
                elif junctions_gdf.loc[u].keep == False:  
                    u = old_u    
                    line_coords[0] = (nodes_gdf.loc[old_u]['x'], nodes_gdf.loc[old_u]['y'])
                    line_coords[-1] = (junctions_gdf.loc[v]['x'], junctions_gdf.loc[v]['y'])
                else:
                    line_coords[0] = (junctions_gdf.loc[u]['x'], junctions_gdf.loc[u]['y'])
                    line_coords[-1] = (junctions_gdf.loc[v]['x'], junctions_gdf.loc[v]['y'])

        elif ((u == 'NA') & (v =='NA')):  # maintain old_u and old_v
                u = old_u
                v = old_v
                line_coords[0] = (nodes_gdf.loc[old_u]['x'], nodes_gdf.loc[old_u]['y'])
                line_coords[-1] = (nodes_gdf.loc[old_v]['x'], nodes_gdf.loc[old_v]['y'])

        elif ((u == 'NA') & (v != 'NA')) : # maintain old_u
                u = old_u
                line_coords[0] = (nodes_gdf.loc[old_u]['x'], nodes_gdf.loc[old_u]['y'])
                
                if junctions_gdf.loc[v].keep == False:
                    v = old_v
                    line_coords[-1] = (nodes_gdf.loc[v]['x'], nodes_gdf.loc[v]['y'])
                else:
                    line_coords[-1] = (junctions_gdf.loc[v]['x'], junctions_gdf.loc[v]['y'])

        else: #(( u =! 'NA') & (v == 'NA') !: # maintain old_v
                v = old_v
                line_coords[-1] = (nodes_gdf.loc[old_v]['x'], nodes_gdf.loc[old_v]['y'])
                if junctions_gdf.loc[u].keep == False:
                    u = old_u
                    line_coords[0] = (nodes_gdf.loc[u]['x'], nodes_gdf.loc[u]['y'])
                else:
                    line_coords[0] = (junctions_gdf.loc[u]['x'], junctions_gdf.loc[u]['y'])

        geo_line = (LineString([coor for coor in line_coords]))
        
        if u == v: 
            edges_gdf.drop(row[0], axis = 0, inplace = True)
            continue
            
        edges_gdf.set_value(row[0],"u", u)
        edges_gdf.set_value(row[0],"v", v)
        edges_gdf.set_value(row[0],"geometry", geo_line)

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