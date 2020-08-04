#!/usr/bin/env python

# Author: Dogacan S. Ozturk

# Function to calculate integrals.
def calculate_point_integral(numberOfPoints, delta, val):
    intVal = 0.0
    for i in range(numberOfPoints):
        intVal = intVal+(i*delta*val)/(numberOfPoints-1)
    return intVal

# Function to calculate distance between two points on Earth.
def calc_distance(lat1,lat2,lon1,lon2):
    '''
    Function returns the distance between two points on Earth based on their
    geographic coordinates.
    
    Parameters:
    ===========
    lat1: Float latitude of the first point in degrees.
    lat2: Float latitude of the second point in degrees.
    lon1: Float longitude of the first point in degrees.
    lon2: Float longitude of the second point in degrees.
    
    Returns:
    ========
    dist: Float distance between the two points in km.
    
    Example:
    ========
    >> from hime_helper_functions import calc_distance
    >> distance = calc_distance(lat1, lat2, lon1, lon2)

    '''
    
    import numpy as np
    
    a = np.sin(np.deg2rad(lat1-lat2)/2.)**2+np.cos(np.deg2rad(lat1))*np.cos(np.deg2rad(lat2))*np.sin(np.deg2rad(lon1-lon2)/2)**2
    c = 2*np.arcsin(a**0.5)
    dist = c*6371
    return dist

# Function to calculate averages.
def smooth(y,box_pts):
    import numpy as np
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
