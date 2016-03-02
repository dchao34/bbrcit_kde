Library for kernel density estimation. 
===

Notable features
---

+ Single and dual-tree evaluation with guaranteed relative and absolute tolerance. The algorithm is based on A. Gray and A. Moore's [paper](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.215.5303), and the implementation is based on a Kdtree. 
+ Gauranteed relative and absolute tolerance. 
+ Arbitrary dimensions. 
+ Arbitrary single bandwidth kernels, but any new additions must conform to the specified API. 
+ Product kernels. In this case, one can choose different bandwidths for the different dimensions. 
+ Weighted reference points. 
