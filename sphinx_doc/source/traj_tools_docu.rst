
.. _traj_tools_docu:

========================================================
Log for traj_tools module
========================================================

This module is designed to do the following tasks:

* Read Model Output (Trajectories), in bin format, save it as Pickle file (this step will be unnecessary once 
new output will be in NetCDF format) --> loadbin
* Filter trajectories by various criteria, e.g. ascent in certain time, or location --> filters
* Compute statistics (not yet implemented or needed)
* Create plots: XY plots of trajectories, histograms of ascent times, and others later --> plots

* common module contains helper functions to load case specifics, etc. 



Version 1.0 - 6.5.2013
-------------------------
Starting Version numbering. 
















To Do
-----------------

**Major Problems**

* Module Structure, Filtering trajectories for use in plots is not efficient or convenient, eg sorting by start indices
* Eliminate global variables, think of good way of using case varibales!


**Smaller Issues**

* Use COSMO Variable names for selecting trajectory tracers
* Save tracers in seperate files to speed up reading process
* Memory leak could be multiple loading of global variables in LoadCaseSpec()
