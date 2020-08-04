<img src="./Examples/Images/hime_logo.png" width=250 height=250 />

# HIME: High-latitude Input for Meso-scale Electrodynamics 





## Overview

Understanding space weather effects and forecasting the dynamics of the Earth’s upper atmosphere (the ionosphere-thermosphere system) requires advanced modelling efforts and better specifications of inputs and boundary conditions for existing numerical models. The High-latitude Input for Meso-scale electrodynamics (HIME) framework is developed in order to determine the upper boundary conditions (drivers) for global ionosphere-thermosphere models. It is based on two dimensional estimates of electric fields derived directly from regional observations. HIME currently builds upon the existing global empirical model of ionospheric electric field potentials (i.e.: Weimer05). The local estimates of electric fields from observations are converted to potentials in HIME and can be merged with any global electric potential information provided in a gridded structure.

HIME uses a novel technique developed by Ozturk et al. (2020) that allows incorporating measurements by incoherent scatter radars and local electric field estimates into global modelling of the ionosphere-thermosphere system. The HIME framework allows for driving numerical models with realistic drivers over selected regions, rather than statistical averages that are traditionally used for modelling. This approach enables modelling of ionospheric electrodynamics with better spatial and temporal resolutions based on better resolved drivers. Our first results demonstrated improved estimates of local heating that is very important for understanding space weather.

## License

Copyright (c) 2020, California Institute of Technology ("Caltech") and University of Michigan. U.S. Government sponsorship acknowledged. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list      of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of Caltech nor its operating division, the Jet Propulsion Laboratory, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Developers

Author: Dogacan S. Ozturk <br>
Co-authors: Xing Meng, Olga P. Verkhoglyadova, Aaron J. Ridley

Contact: dcsoztrk@umich.edu <br>
Institutional contacts: Olga Verkhoglyadova (olga.verkhoglyadova@jpl.nasa.gov), Xing Meng (xing.meng@jpl.nasa.gov)

## Dependencies
The requirements for running the HIME Framework are listed in Requirements.txt. <br>

#### For MAC users:
The dependencies are provided in python_install_dependencies_for_mac.sh. Running: 

sudo -H ./python_install_dependencies_for_mac.sh 

in a Terminal window will install the dependencies.

## Usage

The detailed information can be found in the HIME Tutorial provided inside the Examples folder.

## Acknowledgements
Authors gratefully acknowledge various resources used for the HIME Framework. 

### AMISR
This release of the HIME Framework uses data from the Advanced Modular Incoherent Scatter Radar operated by SRI with a cooperative agreement with NSF (AGS-1840962). The data is provided by Roger Varney and Ashton Reimer. For further information on the data use and access, we refer users to the AMISR website.

References: <br>
[1] Nicolls, M. J., Cosgrove, R., & Bahcivan, H. (2014). Estimating the vector electric field using monostatic, multibeam incoherent scatter radar measurements. Radio Science, 49, 1124-1139. doi: 10.1002/2014RS005519

URL for the AMISR Webpage: https://amisr.com/amisr/

### GITM
This release of the HIME Framework predominantly uses the University of Michigan Global Ionosphere Thermosphere Model (GITM) output. For further information on GITM and its rules of the road, we refer users to the references below.

References: <br>
[1] Ridley, A. J., Deng, Y., & Toth, G. (2006). The global ionosphere-thermosphere model. Journal of Atmospheric and Solar-Terrestrial Physics. doi: 10.1016/j.jastp.2006.01.008 <br>

[2] Zhu, J., Ridley, A. J., & Deng, Y. (2016). Simulating electron and ion temperature in a global ionosphere thermosphere model: Validation and modeling an idealized sub- storm. Journal of Atmospheric and Solar-Terrestrial Physics, 138-139, 243-260. doi: 10.1016/j.jastp.2016.01.005 <br>

URL for the GITM Github repository: https://github.com/aaronjridley/GITM

### Weimer Model
This release of the HIME Framework uses the Fortran version of the Weimer Model [2005] in GITM for obtaining the electric potentials. For further information on the Weimer Model and its rules of the road, we refer users to the references below.

References: <br>
[1] Weimer, D. R., Improved ionospheric electrodynamic models and application to calculating Joule heating rates,  J. Geophys. Res., 110, A05306, doi:10.1029/2004JA010884, 2005. <br>

[2] Weimer, D. R., Predicting Surface Geomagnetic Variations Using Ionospheric Electrodynamic Models,  J. Geophys. Res., 110, A12307, doi:10.1029/2005JA011270, 2005. <br>

URL for the Weimer model: https://zenodo.org/record/2530324

### Spacepy
This release of the HIME Framework uses the Spacepy package to read the binary GITM output. For further information on the Spacepy package and its rules of the road, we refer users to the references below.

References: <br>
[1] Morley, S. K., Koller, J., Welling, D. T., Larsen, B. A., Henderson, M. G., and Niehof, J. T. (2011), Spacepy - A Python-based library of tools for the space sciences, Proceedings of the 9th Python in science conference (SciPy 2010), Austin, TX <br>

URL for the Spacepy Github repository: https://github.com/spacepy/spacepy

### Apexpy
This release of the HIME Framework uses the Apexpy package for coordinate transformation. For further information on the Apexpy package and its rules of the road, we refer users to the references below.

References: <br>
[1] Emmert, J. T., A. D. Richmond, and D. P. Drob (2010), A computationally compact representation of Magnetic-Apex and Quasi-Dipole coordinates with smooth base vectors, J. Geophys. Res., 115(A8), A08322, doi:10.1029/2010JA015326. <br>

[2] Richmond, A. D. (1995), Ionospheric Electrodynamics Using Magnetic Apex Coordinates, Journal of geomagnetism and geoelectricity, 47(2), 191–212, doi:10.5636/jgg.47.191. <br>

URL for the Apexpy Github repository: https://github.com/aburrell/apexpy

### NASA High End Computing Capability
Resources supporting this work were provided by the NASA High-End Computing (HEC) Program through the NASA Advanced Supercomputing (NAS) Division at the Ames Research Center.

URL for the NASA HECC Resources: https://www.nas.nasa.gov/hecc/#url

### Python Packages
This release of the HIME Framework uses various Python packages for numerical calculations and visualization. For further information on these packages and their rules of the road, we refer users to the references below.

References: <br>
[1] K. Jarrod Millman and Michael Aivazis. Python for Scientists and Engineers, Computing in Science & Engineering, 13, 9-12 (2011), doi:10.1109/MCSE.2011.36 <br>

[2] Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011), doi:10.1109/MCSE.2011.37 <br>

[3] John D. Hunter. Matplotlib: A 2D Graphics Environment, Computing in Science & Engineering, 9, 90-95 (2007), doi:10.1109/MCSE.2007.55 <br>


URL for the Scipy package: https://www.scipy.org/ <br>

URL for the Numpy package: https://numpy.org/ <br>

URL for the Matplotlib package: https://matplotlib.org/ <br>

## Attribution

For publishing any work using the HIME Framework, please provide appropriate credit to the developers via citation or acknowledgement.

Ozturk, D. S., Verkhoglyadova, O., Meng, X., Semeter, J., Varney, R., Reimer, A. (2020), A new framework to incorporate high-latitude input for meso-scale electrodynamics: HIME, Journal of Geophysical Research: Space Physics, doi.org/10.1029/2019JA027562

## Issues and Bug Reporting


For issues, custom requests, and bug reporting contact dogacan.s.ozturk@jpl.nasa.gov. Please provide as much as information as possible. Bug reports can include screenshots or copies of error messages.
