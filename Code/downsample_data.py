#!/usr/bin/env python

# Author: Dogacan S. Ozturk

def downsample_pfisr_data(filename, gridRes):

    '''
    Given a filename and the grid resolution, this function downsamples
    the electric field estimates and returns the new electric fields as well as
    electric field potentials on the new uniform grid.
    
    The function increases the grid size and fills the boundary with zeros to
    prevent loss of data.
    
    The function uses a Forward Euler method to calculate the potentials and
    a central differencing algorithm to calculate electric fields. The function
    returns these new electric fields as well as the downsampled electric fields
    to provide the user with the opportunity to inspect the output.
    
    Parameters:
    ===========
    filename: Full path of the PFISR 2D VEF file.
    gridRes:  Desired gridResolution, should be a float and larger than the
              resolution of the PFISR 2D VEF grids.
    
    Returns:
    ========
    PhiX:            Numpy array of Potential change in the X direction on the
                     new grid.
    PhiY:            Numpy array of Potential change in the Y direction on the
                     new grid.
    Ex_downsampled:  Numpy array of Ex values downsampled on the new grid.
    Ey_downsampled:  Numpy array of Ey values downsampled on the new grid.
    Ex_calculated:   Numpy array of Ex values calculated on the downsampled grid.
    Ey_calculated:   Numpy array of Ey values calculated on the downsampled grid.
    xxNew:           Numpy array of longitude values for the new grid.
    yyNew:           Numpy array of latitude values for the new grid.
    experimentTimes: Datetime array corresponding to the PFISR sampling times.
    
    Example:
    ========
    >> import downsample_pfisr_data
    >> PhiX, PhiY, Ex_downsampled, Ey_downsampled, Ex_calculated, Ey_calculated,
       xxNew, yyNew, experimentTimes = downsample_pfisr_data('./Examples/pfisr', 0.3)
    '''
    
    # Import Python libraries.
    import datetime as dt
    import tables
    import numpy as np
    import sys
    
    # Import custom libraries.
    from hime_helper_functions import calc_distance
    
    startOfTimes = dt.datetime(1970,1,1) # UnixTime start.
    
    h5file = tables.open_file(filename) # Read in the file.
    pfisrData = h5file.get_node('/')

    # Read the electric field data.
    Ex_geo = pfisrData['Fit2D']['Ex_geo'].read()
    Exg = np.nan_to_num(Ex_geo)

    Ey_geo = pfisrData['Fit2D']['Ey_geo'].read()
    Eyg = np.nan_to_num(Ey_geo)

    Ex_mag = pfisrData['Fit2D']['Ex'].read()
    Exm = np.nan_to_num(Ex_mag)

    Ey_mag = pfisrData['Fit2D']['Ey'].read()
    Eym = np.nan_to_num(Ey_mag)

    # Read the grid data.
    xgeo = np.nan_to_num(pfisrData['Grid']['X_geo'].read())
    ygeo = np.nan_to_num(pfisrData['Grid']['Y_geo'].read())

    xmag = np.nan_to_num(pfisrData['Grid']['X'].read())
    ymag = np.nan_to_num(pfisrData['Grid']['Y'].read())

    # Get the dimensions.
    nTimes = np.shape(Exm)[0]
    nX = np.shape(Exm)[1]
    nY = np.shape(Exm)[2]
    
    #                        START DOWNSAMPLING                                #

    # Define adapted grids.
    naX = int((int(xmag.max())-int(xmag.min()))/gridRes+1)+2
    naY = int((int(ymag.max())-int(ymag.min()))/gridRes+1)+2

    XadaptGrids = np.linspace(int(xmag.min()-gridRes),int(xmag.max()+gridRes),naX)
    YadaptGrids = np.linspace(int(ymag.min()-gridRes),int(ymag.max()+gridRes),naY)

    # Allocate arrays to store potential and electric field values.
    PhiX = np.zeros((nTimes, naX, naY))
    PhiY = np.zeros((nTimes, naX, naY))

    Exm_adapted = np.zeros((nTimes, naX, naY))
    Eym_adapted = np.zeros((nTimes, naX, naY))
    
    # Downsample and calculate the new potentials.
    for t in range(nTimes):
        for i in range(1,naX-1):
            for j in range(1,naY-1):
                dx = calc_distance(YadaptGrids[j],YadaptGrids[j], XadaptGrids[i+1],XadaptGrids[i-1])
                dy = calc_distance(YadaptGrids[j+1],YadaptGrids[j-1], XadaptGrids[i],XadaptGrids[i])
                indY = np.where(np.abs(xmag[0,:]-XadaptGrids[i]) == np.abs(xmag[0,:]-XadaptGrids[i]).min())[0][0]
                indX = np.where(np.abs(ymag[:,0]-YadaptGrids[j]) == np.abs(ymag[:,0]-YadaptGrids[j]).min())[0][0]
                
                # Forward Euler to calculate potentials.
                PhiX[t,i+1,j] = PhiX[t,i-1,j]-(Exm[t,indX,indY])*dx
                PhiY[t,i,j+1] = PhiY[t,i,j-1]-(Eym[t,indX,indY])*dy
                
                # Store electric field values on the new grid.
                Exm_adapted[t,i,j] = Exm[t,indX,indY]
                Eym_adapted[t,i,j] = Eym[t,indX,indY]


        # Boundary conditions
        PhiX[t,0,:] = PhiX[t,1,:]
        PhiX[t,-1,:] = PhiX[t,-2,:]
        PhiY[t,:,0] = PhiY[t,:,1]
        PhiY[t,:,-1] = PhiY[t,:,-2]
        
    # Calculate the electric fields back from the potentials on the downsampled
    # grid.
    Ex_calculated = np.zeros((nTimes, naX, naY))
    Ey_calculated = np.zeros((nTimes, naX, naY))

    for t in range(nTimes):
        for i in range(1,naX-1):
            for j in range(1,naY-1):
                dx = calc_distance(YadaptGrids[j],YadaptGrids[j], XadaptGrids[i+1],XadaptGrids[i-1])
                dy = calc_distance(YadaptGrids[j+1],YadaptGrids[j-1], XadaptGrids[i],XadaptGrids[i])
                
                # Central differencing to calculate the new Electric Fields.
                Ex_calculated[t,i,j] = -(PhiX[t,i+1,j] - PhiX[t,i-1,j])/dx
                Ey_calculated[t,i,j] = -(PhiY[t,i,j+1] - PhiY[t,i,j-1])/dy

    # Store the values by averaging the experiment times.
    pfisrTimes = np.zeros(nTimes,dtype=object)
    for k in range(nTimes):
        pfisrTimes[k] = startOfTimes+dt.timedelta(seconds=int(pfisrData['Time']['UnixTime'][k][0]))
        
    return(PhiX, PhiY, Exm_adapted, Eym_adapted, Ex_calculated, Ey_calculated, XadaptGrids, YadaptGrids, pfisrTimes)
