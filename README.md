Library for kernel density estimation. 
===

Although a CUDA capable GPU is not required to use this library, we highly recommend it for large problem sizes. 

Supported features and algorithms
---

+ Direct kernel density computation. 
+ Single and dual-tree evaluation with guaranteed relative and absolute tolerance. The algorithm is based on A. Gray and A. Moore's [paper](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.215.5303), and the implementation is based on a Kdtree. 
+ GPU accelerated for sections that require N-body computations. The algorithm used here is similar to those described in [GPU Gems 3](http://http.developer.nvidia.com/GPUGems3/gpugems3_ch31.html). 
+ Adaptive kernel density. The algorithm used here is the same as those described in section 5.3 of [Silverman's book](http://http.developer.nvidia.com/GPUGems3/gpugems3_ch31.html). 
+ Likelihood cross validation. 
+ Least squares cross validation. Two flavors:
  1. Convolution kernel. 
  2. Numerical integration. Currently only for 2d data. 
+ Grid search tool for least squares cross validation. 
+ Arbitrary dimensions kernel densities. For GPU evaluations, one must provide the appropriate Point`N`d specializations if it is not already provided. At present, 1D and 2D estimators work out of the box.  
+ Arbitrary kernels, but any new additions must conform to the specified API. At present, we provide the following out of the box:
  + EpanechnikovKernel.h: Single bandwidth Epanechnikov kernel in arbitrary dimensions. 
  + GaussianKernel.h: Single bandwidth Gaussian kernel in arbitrary dimensions. 
  + EpanechnikovProductKernel2d.h: Epanechnikov product kernel in 2 dimensions. 
  + GaussianProductKernel2d.h: Gaussian product kernel in 2 dimensions. 
+ Weighted contributions of reference points. 
+ Simulating points from the kernel density. 
