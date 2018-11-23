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

"""
Series of utilies for plotting LineString, Points or Polygons geodataframes, and other operations.

"""

## Plot
	
    
def plot_points(gdf, column, classes = 7, ms = 0.9, ms_col = None, title = 'Plot', scheme = 'fisher_jenks', cmap = 'Greys_r', legend = False, bb = True, f = 10):
    
    gdf.sort_values(by = column, ascending = True, inplace = True)    
    f, ax = plt.subplots(1, figsize=(f, f))
    
    # background black or white
    if bb == True: tcolor = 'white'
    else: tcolor = 'black'
    rect = f.patch    
    if bb == True: rect.set_facecolor('black')
    else: rect.set_facecolor('white')
    
    f.suptitle(title, color = tcolor) 
    plt.axis('equal')
    ax.set_axis_off()
    if (ms_col != None): ms = gdf[ms_col]

    if scheme != None:
        gdf.plot(ax = ax, column = column, k = classes, cmap = cmap, s = ms, scheme = scheme, legend = legend)
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(vmin = gdf[column].min(), vmax = gdf[column].max()))
        if legend == True:
            leg = ax.get_legend()  
            leg.get_frame().set_linewidth(0.0)
                
    else: gdf.plot(ax = ax, linewidth = lw, color = 'black')

    plt.show()
    
    
    
    
def plot_lines(gdf, classes = 7, lw = 0.9, column = None, title = 'Plot', scheme = None, cmap = 'Greys_r', legend = False, bb = True, f=10):
    
    if column != None: gdf.sort_values(by = column, ascending = True, inplace = True)    
    f, ax = plt.subplots(1, figsize=(f, f))
    if bb == True: tcolor = 'white'
    else: tcolor = 'black'
    rect = f.patch    
    if bb == True: rect.set_facecolor('black')
    else: rect.set_facecolor('white')

    f.suptitle(title, color = tcolor) 
    plt.axis('equal')
    ax.set_axis_off()
    
    if (column != None) & (scheme == None):
        gdf.plot(ax = ax, column = column, cmap = cmap, linewidth = lw, legend = legend) 
    
    elif scheme == "LynchBreaks":
        bins = [0.12, 0.25, 0.50, 0.75, 1.00]
        cl = ps.User_Defined(gdf[column], bins)
        gdf.assign(cl = cl.yb).plot(ax = ax, column= 'cl', categorical = True, k = 5, cmap = cmap, linewidth = lw, legend=True)
    
        leg = ax.get_legend()
        leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
        leg.get_texts()[0].set_text('0.00 - 0.12')
        leg.get_texts()[1].set_text('0.12 - 0.25')
        leg.get_texts()[2].set_text('0.25 - 0.50')
        leg.get_texts()[3].set_text('0.50 - 0.75')
        leg.get_texts()[4].set_text('0.75 - 1.00')
    
    elif scheme != None:
        gdf.plot(ax = ax, column = column, k = classes, cmap = cmap, linewidth = lw, scheme = scheme, legend = legend)
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(vmin = gdf[column].min(), vmax = gdf[column].max()))
        if legend == True:
            leg = ax.get_legend()  
            leg.get_frame().set_linewidth(0.0)
        
#         sm._A = []
#         cb=plt.colorbar()
#         f.colorbar(sm)
        
    else: gdf.plot(ax = ax, linewidth = lw, color = 'black')
    
     
    plt.show()

def plot_lines_aside(gdf, classes = 7, lw = 0.9, column = None, column_a = None, title = 'Plot', scheme = 'fisher_jenks', cmap = 'Greys_r', fcolor = 'white', legend = False, black_back = True):
        
    if column != None: gdf.sort_values(by = column, ascending = True, inplace = True)    
    if column_a != None: gdf.sort_values(by = column, ascending = True, inplace = True)
    f, (ax1, ax2) = plt.subplots(ncols = 2, figsize=(15, 8), facecolor = fcolor)
    if black_back == True: tcolor = 'white'
    else: tcolor = 'black'
    rect = f.patch    
    if black_back == True: rect.set_facecolor('black')
    else: rect.set_facecolor('white') 
        
    f.suptitle(title, color = tcolor) 
    plt.axis('equal')
    ax2.set_axis_off()
    ax1.set_axis_off()
        
    col = [column, column_a]
    for n, i in enumerate([ax1, ax2]):
        i.set_aspect('equal')
        if (col[n] != None) & (scheme == None):
            gdf.plot(ax = i, column = col[n], linewidth = lw, legend = legend) 

        elif scheme == "LynchBreaks":
            bins = [0.12, 0.25, 0.50, 0.75, 1.00]
            cl = ps.User_Defined(gdf[col[n]], bins)
            gdf.assign(cl = cl.yb).plot(ax = i, column= 'cl', categorical = True, k = 5, cmap = cmap, linewidth = lw, legend=True)

            leg = i.get_legend()
            leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
            leg.get_texts()[0].set_text('0.00 - 0.12')
            leg.get_texts()[1].set_text('0.12 - 0.25')
            leg.get_texts()[2].set_text('0.25 - 0.50')
            leg.get_texts()[3].set_text('0.50 - 0.75')
            leg.get_texts()[4].set_text('0.75 - 1.00')
            i.legend(loc=1)

        elif scheme != None:
            gdf.plot(ax = i, column = col[n], k = classes, cmap = cmap, linewidth = lw, scheme = scheme, legend = legend)
            sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(vmin = gdf[col[n]].min(), vmax = gdf[col[n]].max()))
            if legend == True:
                leg = i.get_legend()  
                leg.get_frame().set_linewidth(0.0)
                leg.set_bbox_to_anchor((0.0, 0.0, 0.2, 0.2))

    #         sm._A = []
    #         cb=plt.colorbar()
    #         f.colorbar(sm)

        else:
            gdf.plot(ax = i, linewidth = lw, color = 'black')    
    
#     ax1.legend(loc="upper right")
#     ax2.legend(loc="upper right")
    plt.axis('equal')
    plt.show()
    
    
    
        
def plot_polygons(gdf, classes = 7, column = None, title = 'Plot', scheme = None, cmap = 'Greens_r', black_back = True, legend = False):
    
    f, ax = plt.subplots(1, figsize=(10,10))
    ax.set_axis_off()
    plt.axis('equal')
    
    if black_back == True: tcolor = 'white'
    else: tcolor = 'black'
    rect = f.patch    
    if black_back == True: rect.set_facecolor('black')
    else: rect.set_facecolor('white')
    
    if (column != None) & (scheme == None):
        gdf.plot(ax = ax, column = column)    
    
    elif scheme == "LynchBreaks":
        bins = [0.12, 0.25, 0.50, 0.75, 1.00]
        cl = ps.User_Defined(gdf[column], bins)
        gdf.assign(cl = cl.yb).plot(ax = ax, column= 'cl', categorical = True, k = 5, cmap = cmap, legend=True)
    
        leg = ax.get_legend()
        leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
        leg.get_texts()[0].set_text('0.00 - 0.12')
        leg.get_texts()[1].set_text('0.12 - 0.25')
        leg.get_texts()[2].set_text('0.25 - 0.50')
        leg.get_texts()[3].set_text('0.50 - 0.75')
        leg.get_texts()[4].set_text('0.75 - 1.00')    
    
    elif scheme != None:
        gdf.plot(ax = ax, column = column, k = classes, cmap = cmap,  scheme = scheme, legend = legend)
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(vmin = gdf[column].min(), vmax = gdf[column].max()))
        if legend == True:
            leg = ax.get_legend()  
            leg.get_frame().set_linewidth(0.0)    
    else: gdf.plot(ax = ax, color = 'orange')
    f.suptitle(title, color = tcolor)  
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

# math functions for angle computations
# from Abhinav Ramakrishnan answer in https://stackoverflow.com/a/28261304/7375309

def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]

def ang(lineA, lineB, degree = True):
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
    angle_rad = math.acos(dot_prod/magB/magA)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle_rad)%360
     #if ang_deg-180<=0:
     #    return 180 - ang_deg
    # else:
     #    return ang_deg

    if degree == True: return angle_deg
    else: return angle_rad

def ang_geoline(geolineA, geolineB, degree = False):

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
        

def euclidean_distance(xs, ys, xt, yt):
    """ xs stands for x source and xt for x target """
    return sqrt((xs - xt)**2 + (ys - yt)**2)
	
# preparation functions
	
# landmarks


def get_coord_angle(origin, distance, angle):

    # calculate offsets with light trig
    (disp_x, disp_y) = (distance * math.sin(math.radians(angle)),
                        distance * math.cos(math.radians(angle)))

    return (origin[0] + disp_x, origin[1] + disp_y)



    
