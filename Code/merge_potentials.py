#!/usr/bin/env python

# Author: Dogacan S. Ozturk

def merge_pfisr_with_gitm_potentials(PhiX, PhiY, XadaptGrids, YadaptGrids, pfisrTimes, weimerSimulationList, gridRes, mergeParameter, plotPotentials, savePotentials):
    
    '''
    This function merges the potentials calculated from PFISR estimates with the
    global potential patterns obtained from Weimer driven GITM simulations.
    
    To merge potentials, the user is required to provide a global potential pattern.
    These potentials can be obtained by directly using the Global Ionosphere
    Thermosphere Model (GITM) output driven with Weimer Potentials. This code can work
    with 3DALL, 3DION, 3DUSR, or 3DHME output from GITM, as long as the output
    files have the potential and grid values. If the output doesn't have the
    PotentialY, we suggest creating a dictionary item called 'PotentialY'
    which is a copy of the 'Potential' values, such as:
    gitmSimulationResults['PotentialY'] = gitmSimulationResults['Potential']
    
    User passes a mergeParameter that is the scalar value of the standard
    deviation for Gaussian kernel, between 0.1 to 1.0, 0.1 indicating maximum
    and 1.0 indicating minimum smoothing during merging of the local and global
    parameters.
    
    Parameters:
    ===========
    PhiX:                   A numpy array of calculated potential differences in
                            longitudinal direction, provided in Volts.
    PhiY:                   A numpy array of calculated potential differences in
                            latitudinal direction, provided in Volts.
    XadaptGrids:            A numpy array containing the longitude coordinates
                            of the downsampled grid points in degrees.
    YadaptGrids:            A numpy array containing the latitude coordinates of
                            the downsampled grid points in degrees.
    PfisrTimes:             An array of datetime objects storing the values of
                            averaged PFISR experiments.
    weimerSimulationsList:  String for full path of the simulation results. The
                            files can be 3DALL, 3DUSR, 3DION, or 3DHME as long
                            as the results containt 'Potential' and 'PotentialY'
    gridRes:                A float value of the new uniform grid spacing.
    mergeParameter:         A float value between 0.1 to 1.0 that in summary,
                            defines the degree of smoothing 0.1 being maximum,
                            while 1.0 is the minimum smoothing applied to
                            estimated potentials.
    plotPotentials:         Logical, if True will save plots of the potentials.
    savePotentials:         Logical, if True will save the output files in a
                            format similar to AMIE procedure.
    
    Returns:
    ========
    phiXhime:  Numpy array of merged potentials in longitudinal direction. [V]
    phiYhime:  Numpy array of merged potentials in latitudinal direction. [V]
    himeEx:    Numpy array of electric field values in longitudinal direction
               calculated from the new HIME potentials. [mV/m]
    himeEy:    Numpy array of electric field values in latitudinal direction
               calculated from the new HIME potentials. [mV/m]
    xHimeMesh: Numpy array of longitudinal values of the new global grid. [deg]
    yHimeMesh: Numpy array of latitudinal values of the new global grid. [deg]
    himeTimes: Datetime array of merged potentials.
    
    Example:
    ========
    >> from merge_potentials import merge_pfisr_with_gitm_potentials
    >> weimerSimulationList = glob('../Examples/Files/Simulations/3D*.bin')
    >> mergeParameter = 0.6
    >> plotPotentials = True
    >> savePotentials = True
    >> phiXhime, phiYhime, himeEx, himeEy, xHimeMesh, yHimeMesh, himeTimes =
       merge_pfisr_with_gitm_potentials(PhiX, PhiY, XnewGrids, YnewGrids,
       experimentTimes, weimerSimulationList, gridRes, mergeParameter,
       plotPotentials, savePotentials)
        
    '''
    
    # Import default python libraries.
    import numpy as np
    import datetime as dt
    from apexpy import Apex
    from scipy.ndimage.filters import gaussian_filter
    from scipy.interpolate import griddata
    
    # Import custom python libraries.
    from spacepy.pybats import gitm
    from hime_helper_functions import calc_distance
    from hime_processing_functions import plot_hime_output
    from hime_processing_functions import write_hime_output
    
    # Sort the simulation files.
    weimerSimulationList.sort()
    nFiles = len(weimerSimulationList)
    
    # Define adapted grids.
    naX = np.shape(XadaptGrids)[0]
    naY = np.shape(YadaptGrids)[0]
    
    # Allocate arrays for the new values of electric fields.
    ExMerged = np.zeros((nFiles, naX, naY))
    EyMerged = np.zeros((nFiles, naX, naY))
    
    # Create a big mesh and a median mesh in magnetic local time.
    nres2deg = gridRes*24./360.
    xBigMeshPts = np.arange(0.0, 24.0 + nres2deg, nres2deg)
    xMedMeshPts = np.arange(0.0, 30.0 + nres2deg, nres2deg)
    yBigMeshPts = np.arange(0.0, 90., gridRes)
    
    nxBigMesh = len(xBigMeshPts)
    nyBigMesh = len(yBigMeshPts)
    
    nxMedMesh = len(xMedMeshPts)
    nyMedMesh = len(yBigMeshPts)
    
    yBigMesh, xBigMesh = np.meshgrid(yBigMeshPts, xBigMeshPts)
    yMedMesh, xMedMesh = np.meshgrid(yBigMeshPts, xMedMeshPts)
    
    # Allocate arrays for the HIME potentials and times.
    himePotX = np.zeros((nFiles, nxBigMesh, nyBigMesh))
    himePotY = np.zeros((nFiles, nxBigMesh, nyBigMesh))
    himeTimes = np.zeros(nFiles, dtype=object)
    
    for file in range(nFiles):
        gitmOutput = gitm.GitmBin(weimerSimulationList[file]) # Read the HME output.
        gitmTime = gitmOutput['time'] # Set time.
        himeTimes[file] = gitmTime
        
        # Find matching PFISR experiment time.
        pfisrTimeInd = np.where(np.abs(gitmTime-pfisrTimes)==np.abs(gitmTime-pfisrTimes).min())[0][0]
        
        # Read in simulation mesh.
        gitmLats = np.rad2deg(gitmOutput['Latitude'][:,:,-1])
        gitmLons = np.rad2deg(gitmOutput['Longitude'][:,:,-1])
        
        # Get the dimensions of the simulation mesh.
        nGitmLons = np.shape(gitmLats)[0]
        nGitmLats = np.shape(gitmLats)[1]
        
        # Calculate number of points in simulation grids.
        nPoints = nGitmLats*nGitmLons
        
        # Obtain the magnetic local time coordinates of the grid points.
        dmlon = np.reshape(gitmOutput['Magnetic Longitude'][:,:,-1], [nPoints,1])
        mlt_time = Apex(date=gitmTime)
        dmlon = mlt_time.mlon2mlt(dmlon,gitmTime)
        dmlon[np.where(dmlon<=6.)] = dmlon[np.where(dmlon<=6.)]+24.
        dmlat = np.reshape(gitmOutput['Magnetic Latitude'][:,:,-1], [nPoints,1])
        
        mltPoints = np.zeros([nPoints,2])
        for n in range(nPoints):
            mltPoints[n,0] = dmlon[n]
            mltPoints[n,1] = dmlat[n]
            
        # Obtain values in geographic coordinates.
        geopotx = np.reshape(gitmOutput['Potential'][:,:,-1], [nPoints, 1])
        geopoty = np.reshape(gitmOutput['PotentialY'][:,:,-1], [nPoints, 1])
            
        # Interpolate the new data on the median mesh.
        gitmXonMedMesh = griddata(mltPoints,geopotx,(xMedMesh,yMedMesh), method='linear')[:,:,0]
        gitmYonMedMesh = griddata(mltPoints,geopoty,(xMedMesh,yMedMesh), method='linear')[:,:,0]
        
        # Create a new mesh to store values of the interpolated potentials.
        gitmXonBigMesh = np.zeros([nxBigMesh,nyBigMesh])
        gitmYonBigMesh = np.zeros([nxBigMesh,nyBigMesh])
        
        for lon in range(nxMedMesh):
            for lat in range(nyMedMesh):
                if(xMedMesh[lon,lat]>=24.):
                    xActual = xMedMesh[lon,lat]-24.0
                else:
                    xActual = xMedMesh[lon,lat]
                    
                indX = np.where(np.abs(xActual - xBigMeshPts) == np.abs(xActual - xBigMeshPts).min())
                gitmXonBigMesh[indX,lat] = gitmXonMedMesh[lon,lat]
                gitmYonBigMesh[indX,lat] = gitmYonMedMesh[lon,lat]
        
        # Place the potentials estimated from PFISR measurements on the
        # global simulation mesh.
        
        pfisrXonBigMesh = np.zeros([nxBigMesh, nyBigMesh])
        pfisrYonBigMesh = np.zeros([nxBigMesh, nyBigMesh])
        
        for lon in range(naX):
            for lat in range(naY):
                pointLat = YadaptGrids[lat]
                pointLon = mlt_time.mlon2mlt(XadaptGrids[lon], gitmTime)
                
                indLat = np.where(np.abs(pointLat - yBigMeshPts) == np.abs(pointLat - yBigMeshPts).min())[0]
                
                indLon = np.where(np.abs(pointLon - xBigMeshPts) == np.abs(pointLon - xBigMeshPts).min())[0]
                
                pfisrXonBigMesh[indLon, indLat] = PhiX[pfisrTimeInd, lon, lat]
                pfisrYonBigMesh[indLon, indLat] = PhiY[pfisrTimeInd, lon, lat]
        
        # Merge the calculated potentials with simulated potentials, both
        # on the simulation mesh.
        mergedPotX = gaussian_filter(gitmXonBigMesh+pfisrXonBigMesh,mergeParameter)
        mergedPotY = gaussian_filter(gitmYonBigMesh+pfisrYonBigMesh,mergeParameter)
        
        # Clean the nan values that arise from interpolation.
        mergedPotX[np.isnan(mergedPotX)] = 0.0
        mergedPotY[np.isnan(mergedPotY)] = 0.0
        
        # Plot potentials if user set plotPotentials to true.
        if(plotPotentials):
            plot_hime_output(xBigMesh, yBigMesh, mergedPotX, mergedPotY, gitmTime)
        
        # Save potentials if user set savePotentials to true.
        if(savePotentials):
            write_hime_output(xBigMesh, yBigMesh, mergedPotX, mergedPotY, gitmTime)
        
        # Store the values of the merged potentials.
        himePotX[file] = mergedPotX
        himePotY[file] = mergedPotY
        
        # Calculate the electric fields using central differencing at this step
        # again for inspecting the results.
        for i in range(1,naX-1):
            for j in range(1,naY-1):
                pointLat = YadaptGrids[j]
                pointLon = mlt_time.mlon2mlt(XadaptGrids[i], gitmTime)

                indLat = np.where(np.abs(pointLat - yBigMeshPts) == np.abs(pointLat - yBigMeshPts).min())[0]
            
                indLon = np.where(np.abs(pointLon - xBigMeshPts) == np.abs(pointLon - xBigMeshPts).min())[0]
                
                dx = calc_distance(YadaptGrids[j],YadaptGrids[j],XadaptGrids[i+1],XadaptGrids[i-1])
                dy = calc_distance(YadaptGrids[j+1],YadaptGrids[j-1],XadaptGrids[i],XadaptGrids[i])
            
                if (indLon >= (nxBigMesh-1)):
                    ExMerged[file, i, j] = -(mergedPotX[indLon, indLat] - mergedPotX[indLon-1, indLat])/(dx/2.)
                elif (indLat >= (nyBigMesh-1)):
                    EyMerged[pfisrTimeInd, i, j] = -(mergedPotY[indLon, indLat+1] - mergedPotY[indLon, indLat-1])/(dy/2.)
                else:
                    ExMerged[file, i, j] = -(mergedPotX[indLon+1, indLat] - mergedPotX[indLon-1, indLat])/dx
                    EyMerged[file, i, j] = -(mergedPotY[indLon, indLat+1] - mergedPotY[indLon, indLat-1])/dy
                    
    return(himePotX, himePotY, ExMerged, EyMerged, xBigMesh, yBigMesh, himeTimes)
    
