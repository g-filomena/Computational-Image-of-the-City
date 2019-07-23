import matplotlib as mp, pandas as pd, numpy as np, geopandas as gpd
import functools
import math
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pysal as ps
import random
import pylab
import matplotlib.colors as cols
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Point, LineString, MultiLineString
from numpy.random import randn
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.patches as mpatches
import sys
from time import sleep
pd.set_option('precision', 10)

"""
Series of utilies for plotting LineString, Points or Polygons geodataframes, and other operations.

"""

## Plotting
	
    
def plot_points(gdf, column, classes = 7, ms = 0.9, ms_col = None, title = 'Plot', scheme = 'fisher_jenks', cmap = 'Greys_r', 
                legend = False, bb = True, line_back = pd.DataFrame({'a' : []}), f = 10):
        
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
    
    if line_back.empty == False:
        line_back.plot(ax = ax, color = 'white', linewidth = 1.1, alpha = 0.3)
    
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
        
    elif scheme != None: gdf.plot(ax = ax, column = column, k = classes, cmap = cmap, s = ms, scheme = scheme, legend = legend, alpha = 1)
    else: gdf.plot(ax = ax, linewidth = lw, color = 'black')
        
    if legend == True:
        leg = ax.get_legend()  
        leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
        leg.get_frame().set_linewidth(0.0) # remove legend border
        leg.set_zorder(102)
        if bb == True:
            for text in leg.get_texts(): text.set_color("white")   

    plt.show() 
    
    
def plot_lines(gdf, column = None, classes = 7, lw = 1.1, title = 'Plot', scheme = None, cmap = 'Greys_r', 
               legend = False, cb = False, bb = True, f=10):
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
            colors = ['white', 'red']
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
        
    elif scheme == "pedestrians":
        bins = [10, 30, 50, 70, 90, 110, 130, 150, 170]
        cl = ps.User_Defined(gdf[column], bins)
        gdf.assign(cl = cl.yb).plot(ax = ax, column= 'cl', categorical = True, k = 9, cmap = cmap, linewidth = lw, legend=True)
        
        # Lynch legend
        leg = ax.get_legend()
        leg.get_texts()[0].set_text('0 - 10')
        leg.get_texts()[1].set_text('10 - 30')
        leg.get_texts()[2].set_text('30 - 50')
        leg.get_texts()[3].set_text('50 - 70')
        leg.get_texts()[4].set_text('70 - 90')
        leg.get_texts()[5].set_text('90 - 110')
        leg.get_texts()[6].set_text('110 - 130')
        leg.get_texts()[7].set_text('130 - 150')
        leg.get_texts()[8].set_text('150 - 170')
    
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

    
        
def plot_polygons(gdf, column = None, classes = 7, title = 'Plot', scheme = None, cmap = 'Greens_r', legend = False, 
                  cb = False, bb = True,  f = 10):
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
    
def plot_lines_aside(gdf, gdf_c = None, classes = 7, lw = 0.9, columnA = None, columnB = None, title = 'Plot', 
                     scheme = 'fisher_jenks', bins = None, cmap = 'Greys_r', legend = False, bb = True, f = 15):
    
    # background black or white - basic settings 
    
    if columnA != None: gdf.sort_values(by = columnA, ascending = True, inplace = True)    

    fcolor = 'white'
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
    col = [columnA, columnB]
    
    for n, i in enumerate([ax1, ax2]):
        if (n == 1) & (columnB != None): gdf.sort_values(by = columnB, ascending = True, inplace = True)  
        i.set_aspect('equal')
        
        if (col[n] != None) & (scheme == None): gdf.plot(ax = i, linewidth = lw, color = 'black') # categorical map
              
        # LynchBreaks scheme
        elif scheme == "LynchBreaks":
            bins = [0.12, 0.25, 0.50, 0.75, 1.00]
            cl = ps.User_Defined(gdf[col[n]], bins)
            print(cl.yb)
            gdf.assign(cl = cl.yb).plot(ax = i, column = 'cl', categorical = True, k = 5, cmap = cmap, linewidth = lw, legend=True)

            leg = i.get_legend()
            leg.get_texts()[0].set_text('0.00 - 0.12')
            leg.get_texts()[1].set_text('0.12 - 0.25')
            leg.get_texts()[2].set_text('0.25 - 0.50')
            leg.get_texts()[3].set_text('0.50 - 0.75')
            leg.get_texts()[4].set_text('0.75 - 1.00')
        
        elif scheme == "pedestrians":
                       
            gdf['cat'] = 0
            text = []
            
            for en, b in enumerate(bins):
                if en == 0: 
                    text.append("0")
                    continue
                gdf['cat'][(gdf[col[n]] > bins[en-1]) & (gdf[col[n]] < b)] = en
                text.append(str(bins[en-1])+"-"+str(b))
                if en == len(bins)-1:
                    gdf['cat'][(gdf[col[n]] >= b)] = en+1
                    text.append(">"+str(b))
            
            color0 = '#000000'
            color1 = '#800026'
            color2 = '#e31a1c'
            color3 = '#fc4e2a'
            color4 = '#fd8d3c'
            color5 = '#feb24c'
            color6 = '#fed976'
            color7 = '#ffffcc'

            categories = [0,1,2,3,4,5,6,7]
            colors = [color0, color1, color2, color3, color4, color5, color6, color7]
            colordict = dict(zip(categories, colors))
            gdf["Color"] = gdf['cat'].apply(lambda x: colordict[x])
         
#             cl = ps.User_Defined(gdf[col[n]], bins)
            gdf.plot(ax = i, categorical = True, color=gdf.Color,  linewidth = lw, legend=True)
           
            legend_dict = dict(zip(text, colors))
            patchList = []
            for key in legend_dict:
                data_key = mpatches.Patch(color=legend_dict[key], label=key)
                patchList.append(data_key)
            leg = plt.legend(handles=patchList, loc = 3)
            for text in leg.get_texts():
                plt.setp(text, color = 'w')
   
 
        # other schemes
        elif scheme != None: 

            gdf.plot(ax = i, column = col[n], k = classes, cmap = cmap, linewidth = lw, scheme = scheme, legend = legend)
            leg = i.get_legend()
            leg.set_bbox_to_anchor((0., 0., 0.2, 0.2))
        
        else: gdf.plot(ax = i, linewidth = lw, color = 'black')  # plain map     
    
#         if legend == True:
#             i.legend(loc=1)
#             leg = i.get_legend()
#             
#             leg.set_zorder(102)
#             leg.get_frame().set_linewidth(0.0) # remove legend border
#             if bb == True:
#                 for text in leg.get_texts(): text.set_color("white")

    plt.show()
    
def overlap_lines(gdf, gdf_c, lw = 0.9, title = 'Plot', f = 15):
    
    # background black or white - basic settings 
#     fig, ax = plt.subplots(ncols = 1, figsize=(f, f), facecolor = fcolor)
#     fig.suptitle(title, color = tcolor)
        # background black or white - basic settings 
    fig, ax = plt.subplots(ncols = 1, figsize=(f, f),)
    fig.suptitle(title)
    plt.axis('equal')
    ax.set_axis_off()
    gdf.plot(ax = ax, linewidth = lw, color = 'black') # categorical map
    gdf_c.plot(ax = ax, linewidth = lw, color = 'red')
    
#     if bb == True: tcolor = 'white'
#     else: tcolor = 'black'
#     rect = fig.patch    
#     if bb == True: rect.set_facecolor('black')
#     else: rect.set_facecolor('white')   

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

def ang_geoline(geolineA, geolineB, degree = False, deflection = False):
    
    """
    Given two LineStrings it computes the deflection angle between them. Returns value in degrees or radians.
    
    ----------
    Parameters
    geolineA, geolineB: LineString geometries
    degree: Boolean
    deflection: Boolean, if true computes angle of incidence, otherwise angle formed by the two vectors as diverging from the common node.
    
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
    
    if deflection == True:
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
    else:
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
        
    except:
        angle_deg = 0
        angle_rad = 0
        
    if degree == True: return angle_deg
    else: return angle_rad

def ang(lineA, lineB, degree = True):
    # Get nicer vector form
    
    vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
    print(vA, vB)
    # Get dot prod
    dot_prod = dot(vA, vB)

    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    print(dot_prod)
    # Get cosine value
    cos_ = dot_prod/magA/magB
    # Get angle in radians and then convert to degrees
    angle_rad = math.acos(dot_prod/magB/magA)
    print(angle_rad)
    # Basically doing angle <- angle mod 360
    angle_deg = math.degrees(angle_rad)%360
    if degree == True: return angle_deg
    else: return angle_rad
    
def center_line(line_geo, line_geoC): 
    
    line_coordsA = list(line_geo.coords)
    line_coordsB = list(line_geoC.coords)
        
    if ((line_coordsA[0] == line_coordsB[-1]) | (line_coordsA[-1] == line_coordsB[0])): line_coordsB.reverse()  
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


def dist_to_gdf(point, gpd):
    gpd = gpd.copy()
    gpd['Dist'] = gpd.apply(lambda row:  point.distance(row.geometry),axis=1)
    geoseries = gpd.loc[gpd['Dist'].argmin()]
    distance  = geoseries.Dist
    index = geoseries.name
    return distance, index


def merge_lines(geolines):
    
    first = list(geolines[0].coords)
    second = list(geolines[1].coords)
    coords = []

    reverse = False
    if first[0] == second[0]: 
        reverse = True
        first.reverse()
    if first[-1] == second[-1]: second.reverse()
    if first[0] == second[-1]:
        first.reverse()
        second.reverse()
        reverse = True
    
    coords = first + second
    last = second
    for n,i in enumerate(geolines):
        if n < 2: continue
        next_coords = list(i.coords)
        if (next_coords[-1] == last[-1]) :
            next_coords.reverse()
            last = next_coords
            
    if reverse == True: coords.reverse()
    geoline = LineString([coor for coor in coords])
    return(geoline)
            
def print_row(index_column):
    sys.stdout.write('\r')
    sys.stdout.write("at row: "+ str(index_column))
    sys.stdout.flush()
    sleep(0.05)
    
    
            
            
            
            
            
            


    
