# ElechemR
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
This repository contains the source code of the R package ElechemR, available for download at https://cran.r-project.org/package=EleChemr.
Digital simulation of electrochemical processes. Each function allows for implicit and explicit solution of the differential equation using methods like Euler, Backwards implicit, Runge Kutta 4, Crank Nicholson and Backward differentiation formula as well as different number of points for derivative approximation. Several electrochemical processes can be simulated such as: Chronoamperometry, Potential Step, Linear Sweep, Cyclic Voltammetry, Cyclic Voltammetry with electrochemical reaction followed by chemical reaction (EC mechanism) and CV with two following electrochemical reaction (EE mechanism). In update 1.1.0 has been added a general purpose CV function that allow to simulate up to 4 EE mechanism combined with chemical reaction for each species.Update 1.2.0 improved the accuracy of the measurements and allow personalized data resolution for simulation. 

Bibliography regarding this methods can be found in the following texts. 
- Dieter Britz, Jorg Strutwolf (2016) <ISBN:978-3-319-30292-8>. 
- Allen J. Bard, Larry R. Faulkner (2000) <ISBN:978-0-471-04372-0>.
