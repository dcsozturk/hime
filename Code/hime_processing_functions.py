#!/usr/bin/env python

# Author: Dogacan S. Ozturk

def write_hime_output(xmesh,ymesh,potx,poty,timestep):

    '''
    Writes the grid and potential information at each time step in a format
    similar to the AMIE procedure.
    
    After these files are generated, users would be required to use a
    pre-processing script to convert these files to input files for each
    respective GCM. For GITM swmf_to_amie.pro converts these output files to
    GITM compatible input.
    
    Parameters:
    ===========
    xmesh:      Numpy array of longitude values in degrees.
    ymesh:      Numpy array of latitude values in degrees.
    potx:       Numpy array of potential differences in X-dir. in Volts.
    poty:       Numpy array of potential differences in Y-dir. in Volts.
    timestep:   Datetime object for the timestamp of the values.
    
    Returns:
    =========
    The function doesn't return anything.
    
    Example:
    ========
    >> from hime_processing_functions import write_hime_output
    >> write_hime_output(xx, yy, potx, poty, timestep)
    
    '''
    # Import Python libraries.
    import numpy as np
    import os
    
    # Create folders to save the output files.
    savedir = './potentials'
    if not os.path.exists(savedir):
        os.mkdir(savedir)
    
    outdir = savedir + '/hime_input'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    # Get the dimensions of the mesh.
    nYmesh = np.shape(ymesh)[1]
    nXmesh = np.shape(xmesh)[0]
    
    # Format for writing the values
    fmtVals = '{0:13.5f}{1:13.5f}{2:13.5f}{3:13.5f}{4:13.5f}{5:13.5f}\n'
    
    # Name the file for output.
    savefile = 'it{0:%Y%m%d}'.format(timestep)+'_{0:%H%M%S}'.format(timestep)+'_000.idl'
    
    # Start writing the file.
    outputFile = open(outdir+'/'+savefile, 'w')
    outputFile.write('TITLE\n')
    outputFile.write('\t HIME output for PFISR + Weimer potentials at {0:%Y%m%d_%H:%M:%S}'.format(timestep)+'\n')
    outputFile.write('\n')
    outputFile.write('NUMERICAL VALUES\n')
    outputFile.write('\t'+str(6)+' nvars\n')
    outputFile.write('\t'+str(nYmesh)+' nTheta\n')
    outputFile.write('\t'+str(nXmesh)+' nPsi\n')
    outputFile.write('\n')
    outputFile.write('VARIABLE LIST\n')
    outputFile.write('\t1'+' Theta [deg]\n')
    outputFile.write('\t2'+' Psi [deg]\n')
    outputFile.write('\t3'+' PhiX [kV]\n')
    outputFile.write('\t4'+' PhiY [kV]\n')
    outputFile.write('\t5'+' E-Flux [W/m2]\n')
    outputFile.write('\t6'+' Ave-E [eV]\n')
    outputFile.write('\n')
    outputFile.write('TIME\n')
    outputFile.write('\t{0:%Y}'.format(timestep)+' Year\n')
    outputFile.write('\t{0:%m}'.format(timestep)+' Month\n')
    outputFile.write('\t{0:%d}'.format(timestep)+' Day\n')
    outputFile.write('\t{0:%H}'.format(timestep)+' Hour\n')
    outputFile.write('\t{0:%M}'.format(timestep)+' Minute\n')
    outputFile.write('\t{0:%S}'.format(timestep)+' Second\n')
    outputFile.write('\t00 '+'Milisecond\n')
    outputFile.write('\n')
    outputFile.write('BEGIN NORTHERN HEMISPHERE\n')
    for j in range(nXmesh):
        for k in range(nYmesh):
            outputFile.write(fmtVals.format(90.0-ymesh[j][k],
                                            xmesh[j][k],
                                            potx[j][k]/1000.,
                                            poty[j][k]/1000.,
                                            0.0,
                                            0.0))
    outputFile.close()


def plot_hime_output(xmesh,ymesh,potx,poty,timestep):
    
    '''
    Plots the HIME X and Y potential differentces at the Northern hemisphere as
    polar plots in magnetic local time coordinates.
    
    The plots are saved in the hime_plots directory.
    
    Parameters:
    ===========
    xmesh:      Numpy array of longitude values in degrees.
    ymesh:      Numpy array of latitude values in degrees.
    potx:       Numpy array of potential differences in X-dir. in Volts.
    poty:       Numpy array of potential differences in Y-dir. in Volts.
    timestep:   Datetime object for the timestamp of the values.
    
    Returns:
    =========
    The function doesn't return anything.
    
    Example:
    ========
    >> from hime_processing_functions import plot_hime_output
    >> plot_hime_output(xx, yy, potx, poty, timestep)
    
    '''
    
    # Import Python libraries.
    import os
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    
    # Create folders to save the output files.
    savedir = './potentials'
    if not os.path.exists(savedir):
        os.mkdir(savedir)
    
    plotdir = savedir + '/hime_plots'
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)
    
    # Plot specific parameters.
    plimit = 60.
    nlevels= 25
    plevels = np.linspace(-plimit, plimit,nlevels)
    cMapPotential2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#760797","#C2837D","#F2DCC8","#6EB0A6","#2D61A2"])
    
    # Start plotting.
    fig = plt.figure(figsize=(14,7))
    a1 = plt.subplot(121,polar=True)
    potXcolour = a1.contourf(np.deg2rad(xmesh*15.0), 90.0-ymesh, potx/1000., levels=plevels, extend='both', cmap=cMapPotential2)
    potXline = a1.contour(np.deg2rad(xmesh*15.0), 90.0-ymesh, potx/1000., levels=plevels[::1], colors='k')
    plt.clabel(potXline, inline=1,fontsize=10)
    cbar = fig.colorbar(potXcolour, ax=a1, shrink=0.7)
    cbar.set_label('$\Phi_{X}$ [kV]', fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    a1.set_theta_zero_location('S')
    a1.set_xticks(np.deg2rad([0,90,180,270]))
    a1.set_xticklabels(['00','06','12','18'], fontsize=12)
    a1.set_yticks([0,10,20,30,40])
    a1.set_yticklabels(['o', '$80^{\circ}$', '$70^{\circ}$', '$60^{\circ}$', '$50^{\circ}$'], fontsize=12)
    a1.set_ylim([0,40])
    a1.grid(color='gray', linestyle='-.', linewidth=1)
    
    a2 = plt.subplot(122, polar=True)
    potYcolour = a2.contourf(np.deg2rad(xmesh*15.0), 90.0-ymesh, poty/1000., levels=plevels, extend='both', cmap=cMapPotential2)
    potYline = a2.contour(np.deg2rad(xmesh*15.0), 90.0-ymesh, poty/1000., levels=plevels[::2], extend='both', colors='k')
    plt.clabel(potYline, inline=1,fontsize=10)
    cbar = fig.colorbar(potYcolour, ax=a2, shrink=0.7)
    cbar.set_label('$\Phi_{Y}$ [kV]', fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    a2.set_theta_zero_location('S')
    a2.set_xticks(np.deg2rad([0,90,180,270]))
    a2.set_xticklabels(['00','06','12','18'], fontsize=12)
    a2.set_yticks([0,10,20,30,40])
    a2.set_yticklabels(['o', '$80^{\circ}$', '$70^{\circ}$', '$60^{\circ}$', '$50^{\circ}$'], fontsize=12)
    a2.set_ylim([0,40])
    a2.grid(color='gray', linestyle='-.', linewidth=1)
    
    plt.suptitle('Input Potentials at {0:%H:%M}'.format(timestep), size=14)
    plt.tight_layout()

    plt.savefig(plotdir + '/hime_patterns_{0:%Y%m%d%H%M}.png'.format(timestep), bbox_inches='tight',dpi=200)
    plt.close()
