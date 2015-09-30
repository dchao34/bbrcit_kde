Example
-------

+ `make`. This will create all the unit tests. 

+ Generate some data. Currently a bimodal gaussian (see source for details):

  ```
  ./generate_bimodal_gauss 10000
  ```

  this will produce `bimodal_gaussian.csv`. 

+ Bandwidth cross validation:

  ```
  ./prodkde2d_fftcv bimodal_gaussian.csv 1
  ./prodkde2d_fftcv bimodal_gaussian.csv 2
  ```

  I got 0.074 for the first dimension and 0.022 for the second. 

+ Scan the product kde's:

  ```
  python call_kde_scan.py bimodal_gaussian.csv 0.074 0.022
  python call_akde_scan.py bimodal_gaussian.csv 0.074 0.022
  ```

  this will produce `prodkde2d_scan2d.dat`, `prodkde2d_scan1d.dat`, `prodadakde2d_scan2d.dat`, and `prodadakde2d_scan2d.dat`.

+ `jupyter notebook`. Then run all of `kde_visual.ipynb`. 
