import matplotlib as mp, pandas as pd, numpy as np, geopandas as gpd
import functools
import math
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pysal as ps
import random
import pylab

from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Point, LineString, MultiLineString

from scipy import sparse
from scipy.sparse import linalg
pd.set_option('precision', 10)

"""
Series of utilies for plotting LineString, Points or Polygons geodataframes, and other operations.

"""

## Plotting
	
    
def plot_points(gdf, column, classes = 7, ms = 0.9, ms_col = None, title = 'Plot', scheme = 'fisher_jenks', cmap = 'Greys_r', legend = False, bb = True, f = 10):
        
    """
    It creates a plot from a points GDF. 
    When column and scheme are not 'None' it plots the distribution over value and geographical space of variable 'column using scheme
    'scheme'. If only 'column' is provided, a categorical map is depicted.
    
    Parameters
    ----------
    gdf: GeoDataFrame
    column: string, column on which the plot is based
    classes: = 7, classes for visualising when scheme is not None and different from 'LynchBreaks'
    ms: float, markersize value (fixed)
    ms_col: String, column name in the GDF'S column where markersize values are stored
    title: string, title of the graph
    scheme: {'fisher_jenks', 'quantiles', 'equal_interval', 'LynchBreaks}
    cmap: String, see matplotlib colormaps for a list of possible values
    legend: boolean
    bb: boolean, black background or white
    f: float, size figure extent
    
    """
    
    # background black or white - basic settings
    if column != None: gdf.sort_values(by = column, ascending = True, inplace = True)    
    fig, ax = plt.subplots(1, figsize=(f, f))
    if scheme == "LynchBreaks": legend = True
    if bb == True: tcolor = 'white'
    else: tcolor = 'black'
    rect = fig.patch    
    if bb == True: rect.set_facecolor('black')
    else: rect.set_facecolor('white')
    fs = f+5 # font-size   
    fig.suptitle(title, color = tcolor, fontsize=fs)
    plt.axis('equal')
    ax.set_axis_off()
    
    # markers size from column
    if (ms_col != None): ms = gdf[ms_col]
        
    if scheme == "LynchBreaks":
        bins = [0.12, 0.25, 0.50, 0.75, 1.00]
        cl = ps.User_Defined(gdf[column], bins)
        gdf.assign(cl = cl.yb).plot(ax = ax, column= 'cl', categorical = True, k = 5, cmap = cmap, linewidth = lw, legend=True)
        
        # Lynch legend
        leg = ax.get_legend()
        leg.get_texts()[0].set_text('0.00 - 0.12')
        leg.get_texts()[1].set_text('0.12 - 0.25')
        leg.get_texts()[2].set_text('0.25 - 0.50')
        leg.get_texts()[3].set_text('0.50 - 0.75')
        leg.get_texts()[4].set_text('0.75 - 1.00')

    elif scheme != None: gdf.plot(ax = ax, column = column, k = classes, cmap = cmap, s = ms, scheme = scheme, legend = legend)
    else: gdf.plot(ax = ax, linewidth = lw, color = 'black')
        
    if legend == True:
        leg = ax.get_legend()  
        leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
        leg.get_frame().set_linewidth(0.0) # remove legend border
        leg.set_zorder(102)
        if bb == True:
            for text in leg.get_texts(): text.set_color("white")   

    plt.show()
    
    
def plot_lines(gdf, column = None, classes = 7, lw = 1.1, title = 'Plot', scheme = None, cmap = 'Greys_r', legend = False, cb = False,     bb = True, f=10):
    """
    It creates a plot from a lineStrings GDF.
    When column and scheme are not 'None' it plots the distribution over value and geographical space of variable 'column using scheme
    'scheme'. If only 'column' is provided, a categorical map is depicted.
    
    Parameters
    ----------
    gdf: GeoDataFrame
    column: string, column on which the plot is based
    classes: = 7, classes for visualising when scheme is not None and different from 'LynchBreaks'
    lw: float, line width (fixed)
    ms_col: String, column name in the GDF'S column where markersize values are stored
    title: string, title of the graph
    scheme: {'fisher_jenks', 'quantiles', 'equal_interval', 'LynchBreaks}
    cmap: String, see matplotlib colormaps for a list of possible values
    legend: boolean
    cb: boolean, show color bar
    bb: boolean, black background or white
    f: float, size figure extent
    
    """
    
    # background black or white - basic settings    
    if column != None: gdf.sort_values(by = column, ascending = True, inplace = True)    
    fig, ax = plt.subplots(1, figsize=(f, f))
    if scheme == "LynchBreaks": legend = True
    if bb == True: tcolor = 'white'
    else: tcolor = 'black'
    rect = fig.patch  
    rect.set_zorder(0)
    if bb == True: rect.set_facecolor('black')
    else: rect.set_facecolor('white')
    fs = f+5 # font-size    
    fig.suptitle(title, color = tcolor, fontsize=fs)
    plt.axis('equal')
    ax.set_axis_off()
    
    # categorigal variable
    if (column != None) & (scheme == None):
        if cmap == None: # boolean map
            colors = ['white', 'black']
            gdf.plot(ax = ax, categorical = True, column = column, color = colors, linewidth = lw, legend = legend) 
        else: gdf.plot(ax = ax, column = column, cmap = cmap, linewidth = lw, legend = legend) # cateogircal map
    
    # LynchBreaks scheme
    elif scheme == "LynchBreaks":
        bins = [0.12, 0.25, 0.50, 0.75, 1.00]
        cl = ps.User_Defined(gdf[column], bins)
        gdf.assign(cl = cl.yb).plot(ax = ax, column= 'cl', categorical = True, k = 5, cmap = cmap, linewidth = lw, legend=True)
        
        # Lynch legend
        leg = ax.get_legend()
        leg.get_texts()[0].set_text('0.00 - 0.12')
        leg.get_texts()[1].set_text('0.12 - 0.25')
        leg.get_texts()[2].set_text('0.25 - 0.50')
        leg.get_texts()[3].set_text('0.50 - 0.75')
        leg.get_texts()[4].set_text('0.75 - 1.00')
    
    elif scheme != None:
        gdf.plot(ax = ax, column = column, k = classes, cmap = cmap, linewidth = lw, scheme = scheme, legend = legend)
    else: gdf.plot(ax = ax, linewidth = lw, color = 'black') # plain map
                      
    if legend == True:
        leg = ax.get_legend()  
        leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
        leg.get_frame().set_linewidth(0.0) # remove legend border
        leg.set_zorder(102)
        if bb == True:
            for text in leg.get_texts(): text.set_color("white")
                
    if cb == True:
        # create color bar
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(vmin = gdf[column].min(), vmax = gdf[column].max()))
        fig = ax.get_figure()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        sm._A = []
        cbar = fig.colorbar(sm, cax=cax, ticks=[gdf[column].min(), gdf[column].max()]) 
        cbar.ax.set_yticklabels(['min', 'max'])
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color=tcolor, fontsize=(fs-5))   
        
    plt.show()

    
        
def plot_polygons(gdf, column = None, classes = 7, title = 'Plot', scheme = None, cmap = 'Greens_r', legend = False, cb = False, bb = True,  f = 10):
    """
    It creates a plot from a polygons GDF.
    When column and scheme are not 'None' it plots the distribution over value and geographical space of variable 'column using scheme
    'scheme'. If only 'column' is provided, a categorical map is depicted.
    
    
    Parameters
    ----------
    gdf: GeoDataFrame
    column: string, column on which the plot is based
    classes: = 7, classes for visualising when scheme is not None and different from 'LynchBreaks'
    title: string, title of the graph
    scheme: {'fisher_jenks', 'quantiles', 'equal_interval', 'LynchBreaks}
    cmap: String, see matplotlib colormaps for a list of possible values
    legend: boolean
    cb: boolean, show color bar
    bb: boolean, black background or white
    f: float, size figure extent
    
    """
    
    # background black or white - basic settings   
    fig, ax = plt.subplots(1, figsize=(f, f))
    if scheme == "LynchBreaks": legend = True
    if bb == True: tcolor = 'white'
    else: tcolor = 'black'
    rect = fig.patch  
    rect.set_zorder(0)
    if bb == True: rect.set_facecolor('black')
    else: rect.set_facecolor('white')
    fs = f+5 # font-size    
    fig.suptitle(title, color = tcolor, fontsize=fs)
    plt.axis('equal')
    ax.set_axis_off()
    
    # categorigal variable
    if (column != None) & (scheme == None): gdf.plot(ax = ax, column = column, cmap = cmap, categorical = True, legend = legend)       
    
    # LynchBreaks scheme
    elif scheme == "LynchBreaks":
        bins = [0.12, 0.25, 0.50, 0.75, 1.00]
        cl = ps.User_Defined(gdf[column], bins)
        gdf.assign(cl = cl.yb).plot(ax = ax, column= 'cl', categorical = True, k = 5, cmap = cmap, legend=True)
        
        # Lynch Legend
        leg = ax.get_legend()
        leg.get_texts()[0].set_text('0.00 - 0.12')
        leg.get_texts()[1].set_text('0.12 - 0.25')
        leg.get_texts()[2].set_text('0.25 - 0.50')
        leg.get_texts()[3].set_text('0.50 - 0.75')
        leg.get_texts()[4].set_text('0.75 - 1.00')    
    
    # Other schemes
    elif scheme != None:
        if cb == True: legend = False # if colorbar is activated, don't show legend
        gdf.plot(ax = ax, column = column, k = classes, cmap = cmap,  scheme = scheme, legend = legend)
        
    else: gdf.plot(ax = ax, color = 'orange')  # plain map
    
    if legend == True:
        leg = ax.get_legend()  
        leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
        leg.get_frame().set_linewidth(0.0) # remove legend border
        leg.set_zorder(102)
        if bb == True:
            for text in leg.get_texts(): text.set_color("white")
                
    if cb == True:
        # create color bar
        sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(vmin = gdf[column].min(), vmax = gdf[column].max()))
        fig = ax.get_figure()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        sm._A = []
        cbar = fig.colorbar(sm, cax=cax, ticks=[gdf[column].min(), gdf[column].max()]) 
        cbar.ax.set_yticklabels(['min', 'max'])
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color=tcolor, fontsize=(fs-5))   
       
    plt.show()    
    
def plot_lines_aside(gdf, classes = 7, lw = 0.9, column = None, column_a = None, title = 'Plot', scheme = 'fisher_jenks', cmap = 'Greys_r', fcolor = 'white', legend = False, bb = True):
    
    # background black or white - basic settings   
    if column != None: gdf.sort_values(by = column, ascending = True, inplace = True)    
    if column_a != None: gdf.sort_values(by = column, ascending = True, inplace = True)
    fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize=(15, 8), facecolor = fcolor)
    if bb == True: tcolor = 'white'
    else: tcolor = 'black'
    rect = fig.patch    
    if bb == True: rect.set_facecolor('black')
    else: rect.set_facecolor('white')   
    fig.suptitle(title, color = tcolor)
    plt.axis('equal')
    ax2.set_axis_off()
    ax1.set_axis_off()
        
    col = [column, column_a]
    
    for n, i in enumerate([ax1, ax2]):
        i.set_aspect('equal')
        if (col[n] != None) & (scheme == None): gdf.plot(ax = i, column = col[n], linewidth = lw, legend = legend) # categorical map
        
        # LynchBreaks scheme
        elif scheme == "LynchBreaks":
            bins = [0.12, 0.25, 0.50, 0.75, 1.00]
            cl = ps.User_Defined(gdf[col[n]], bins)
            gdf.assign(cl = cl.yb).plot(ax = i, column= 'cl', categorical = True, k = 5, cmap = cmap, linewidth = lw, legend=True)

            leg = i.get_legend()
            leg.get_texts()[0].set_text('0.00 - 0.12')
            leg.get_texts()[1].set_text('0.12 - 0.25')
            leg.get_texts()[2].set_text('0.25 - 0.50')
            leg.get_texts()[3].set_text('0.50 - 0.75')
            leg.get_texts()[4].set_text('0.75 - 1.00')
        
        # other schemes
        elif scheme != None: gdf.plot(ax = i, column = col[n], k = classes, cmap = cmap, linewidth = lw, scheme = scheme, legend = legend)
        else: gdf.plot(ax = i, linewidth = lw, color = 'black')  # plain map     
    
        if legend == True:
#             i.legend(loc=1)
            leg = i.get_legend()
            leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
            leg.set_zorder(102)
#             leg.get_frame().set_linewidth(0.0) # remove legend border
            if bb == True:
                for text in leg.get_texts(): text.set_color("white")

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
    It takes a list of dictionary and merge them in a df, as columns.
    
    ----------
    Parameters
    list_dict: list of dictionaries that will become df's columns
    list_col: list of the names that will be used as colums heading
    
    Returns:
    ----------
    DataFrame
    """
    
    df = pd.DataFrame(list_dict).T
    df.columns = ['d{}'.format(i) for i, col in enumerate(df, 1)]
    df.columns = list_col
    
    return(df)

	
# Angle/math functions
	
def euclidean_distance(xs, ys, xt, yt):
    """ xs stands for x source and xt for x target """
    return sqrt((xs - xt)**2 + (ys - yt)**2)

"""
math functions for angle computations
readapted for LineStrings from Abhinav Ramakrishnan answer in https://stackoverflow.com/a/28261304/7375309
"""

def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]


def get_coord_angle(origin, distance, angle):
    """
    Given two LineStrings it computes the deflection angle between them. Returns value in degrees or radians.
    
    ----------
    Parameters
    geolineA, geolineB: LineString geometries
    degree: Boolean
    
    Returns:
    ----------
    float
    """

    (disp_x, disp_y) = (distance * math.sin(math.radians(angle)), distance * math.cos(math.radians(angle)))

    return (origin[0] + disp_x, origin[1] + disp_y)

def ang_geoline(geolineA, geolineB, degree = False):
    
    """
    Given two LineStrings it computes the deflection angle between them. Returns value in degrees or radians.
    
    ----------
    Parameters
    geolineA, geolineB: LineString geometries
    degree: Boolean
    
    Returns:
    ----------
    float
    """
    
    # extracting coordinates and creates lines
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
    
def ang_geoline_n(geolineA, geolineB, degree = False):
    
    """
    Given two LineStrings it computes the deflection angle between them. Returns value in degrees or radians.
    
    ----------
    Parameters
    geolineA, geolineB: LineString geometries
    degree: Boolean
    
    Returns:
    ----------
    float
    """
    
    # extracting coordinates and creates lines
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
        lineA = ((x_tA,y_tA), (x_fA, y_fA))
        lineB = ((x_tB, y_tB),(x_fB, y_fB))

    elif ((x_tA, y_tA) == (x_fB, y_fB)):
        lineA = ((x_tA,y_tA), (x_fA, y_fA))
        lineB = ((x_tB, y_tB),(x_fB, y_fB))

    elif ((x_fA, y_fA) == (x_fB, y_fB)):
        lineA = ((x_fA,y_fA),(x_tA, y_tA))
        lineB = ((x_fB, y_fB),(x_tB, y_tB))

    else: #(from_node == to_node2)
        lineA = ((x_fA,y_fA),(x_tA, y_tA))
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
#         if angle_deg-180 >= 0: angle_deg = 360 - angle_deg
    except:
        angle_deg = 0
        angle_rad = 0
    

        
    if degree == True: return angle_deg
    else: return angle_rad



    
