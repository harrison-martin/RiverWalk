# RiverWalk-Strat - River Avulsion Simulation (Now with Stratigraphy!)
RiverWalk-Strat v1.0.0 - River Avulsion Simulation (Now with Stratigraphy!)

Caitlin Sifuentes, Harrison Martin, and Doug Edmonds

_Correspondence to/Uploaded by:_

Harrison Martin
Postdoctoral Scholar, Division of Geological and Planetary Sciences
California Institute of Technology 
hkm@caltech.edu 

_or_

Douglas Edmonds
Professor, Department of Earth and Atmospheric Sciences
Indiana University
edmondsd@iu.edu

This document was last updated: August 2024

# Introduction 
This MATLAB code is intented to accompany, and reproduce the results of, a manuscript submitted to ESurf by Caitlin Sifuentes[1,2], Harrison K. Martin[1,3], Kyle M. Straub[4], Elizabeth A. Hajek[5], and Douglas A. Edmonds[1].

[1]Department of Earth and Atmospheric Sciences, Indiana University, Bloomington, IN, USA

[2]Geosyntec Consultants, Long Beach, CA, USA

[3]Division of Geological and Planetary Sciences, California Institute of Technology, Pasadena, CA, USA.

[4]Department of Earth and Environmental Sciences, Tulane University, New Orleans, LA, USA

[5]Department of Geosciences, Penn State University, University Park, PA, USA


This code is titled RiverWalk-Strat and it is built on a modified form of RiverWalk, a model written by myself [Harrison Martin] and Doug Edmonds as part of my PhD, which can be found at: https://github.com/harrison-martin/RiverWalk. That model was described in two publications: https://doi.org/10.5194/esurf-10-555-2022 and 
https://doi.org/10.1130/G51138.1. This new version, RiverWalk-Strat is the creation of Caitlin Sifuentes, myself, and Doug Edmonds. Beginning with an earlier version of RiverWalk, it consists of two scripts that run in MATLAB. The first, RiverWalk, generates landscapes as a result of an avulsing river over stratigraphic timescales. This version has been modified to save and export those landscapes at each timestep in a format that can be read by the second script, StratCode. StratCode takes the outputs from RiverWalk and constructs synthetic stratigraphy along specified strike- or dip-oriented cross-sections. It also contains tools for statistically analysing along-strike river position and stratigraphic compensation. 

Please feel free to reach out if you'd like help understanding or using the code, or for any other questions. We are also always happy to chat if you have an interesting research problem that you think the model may be able to help solve. You're also free to modify the code as you wish, in line with the GNU GPL v3.0 license (link below).           

Thanks for reading, and we hope you enjoy the code!
 - Harrison  

# License statement
GNU GPL (General Public License) v3.0.                                  
https://choosealicense.com/licenses/gpl-3.0/#                           
https://csdms.colorado.edu/wiki/License

RiverWalk-Strat v1.0.0 - River Avulsion Simulation (Now with Stratigraphy!) (MATLAB)
Copyright (C) 2021 Caitlin Sifuentes, Harrison Martin, and Douglas Edmonds
Correspondence can be addressed to hkm@caltech.edu or edmondsd@iu.edu.
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
