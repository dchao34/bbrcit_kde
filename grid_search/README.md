Grid search for 2 dimensional data
===

Overview
---

Input: Text file containing of rows of points, each of which is three space separated coluns. 
  + Column 1: `x` coordinate.
  + Column 2: `y` coordinate.
  + Column 3: Point weight `w`; e.g. `1.0` if all contribute equally. 

  See `example.dat` for an example of such a file. 

Output: Depends on configuration file settings. See `template.cfg` for some explanations. 

Procedure
---

1. First create the executable by running `make`. 
2. Set configuration parameters in a config file. See `template.cfg` for an example. 
3. `./grid_search template.cfg`
4. Visualize results. See `example.ipynb` for an example. 
