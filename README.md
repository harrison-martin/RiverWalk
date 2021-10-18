# RiverWalk - River Avulsion Simulation
RiverWalk v1.0.0 - River Avulsion Simulation

Harrison Martin

PhD Candidate, Department of Earth & Atmospheric Sciences 

Indiana University 

hkmartin@iu.edu 

This document was last updated: October 2021 

# Introduction 
This MATLAB code is intented to accompany, and reproduce the results of, a manuscript submitted to ESurf by myself and my advisor Dr. Douglas A. Edmonds (edmondsd@iu.edu), who also contributed to the code.          

This release is intended to be as understandable as possible, but as with any large(r) model, there are some components that might seem opaque or idiosyncratic (because sometimes they are). I have tried to comment everything as clearly as possible, and usually used intuitive variable names, at least for stuff a user would usually interact with. I also tried to ensure that everything in the code works as advertised and is stable, but can't warrant perfect performance in all cases. As such, I removed some old, unusued components or setting options that are unnecessary or irrelevant to reproducing the manuscript's results, because I wasn't able to test old components sufficiently to have a stable release.   

All said, I'm more than happy to provide whatever help I can if you want to use or understand this code! Feel free to reach out to me if you have any questions. I'm also always happy to chat if you have an interesting research problem that you think the model may be able to help solve. You're also free to modify the code as you wish, in line with the GNU GPL v3.0 license (link below). If you're interested in adding or modifying components, I'd love to chat! Collaboration is great. Also, there's a fair chance that I may have tried to implement some version of it at a previous point and have results or code that I would be more than happy to share.            

Thanks for reading, and I hope you enjoy the code!
 - Harrison  

# Pre-requisites

Luckily, the code is pretty plug-and-play. Aside from MATLAB, no other software or files need to be downloaded. You can optionally download a custom colour bar that I created that is the same as parula, but with the first few values set to white instead of deep blue. This just provides more visual contrast betwen low values and zero values. It's included in this directory under the name "whitetip4". If you don't want to use it, just modify the code such that UseCustomColormap = 0 (under codified constants).

The only other step you need to do ahead of time is to specify the directory in which you will work. Modify the exportDirectory = '' line such that it points to an actual directory on your computer. Then, in that directory, create a new folder for each of your entries in the "TitlingForRuns" variable at the top under experimental variables. Then you should be able to hit Run and start generating results!

# do a paragraph describing the code, summary, etc. no dependency. matlab code.

# License statement
GNU GPL (General Public License) v3.0.                                  
https://choosealicense.com/licenses/gpl-3.0/#                           
https://csdms.colorado.edu/wiki/License

RiverWalk v1.0.0 - River Avulsion Simulation (MATLAB)
Copyright (C) 2021 Harrison Martin
Developer can be contacted at hkmartin@iu.edu.
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
