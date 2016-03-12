Unit Tests
==========

The following is a listing of all the unit tests available in this directory. 

Each line below is associated with a single unit test. The binary executable that performs the test appears to the left of the `:`, and a description of the test appears to the right.

For the CPU tests, the binary executables are all produced using the single source file whose name is the same except for the extension `.cc`. For example, `test_kde1` is produced with `test_kde1.cc`. 

For the GPU tests, the binary execuables all have the suffix `_cuda`. They are produced using the source file whose name is the same as the string preceding `_cuda`, but with the extension `.cc`. For example, `test_point2d_cuda` is produced using `test_point2d.cc`. 

A single source file can produce both a CPU as well as a GPU test. In this case, the cpu binary has the same name as the prefix of the source file, while the gpu binary has the same name as the prefix, except with a `cuda` extension. For example, `test_point2d.cc` produces the cpu test `test_point2d` as well as a gpu test `test_point2d_cuda`. 

CPU tests
---------
+ `test_kde10`: Basic use of KernelDensity<>. 
+ `test_kde11`: Runtime benchmark of KernelDensity<> for all points evaluation on the CPU. 
+ `test_kde12`: Same as `test_kde10`, but with product kernels. 
+ `test_kde13`: Basic use of adaptive kernels.
+ `test_kde14`: Same as `test_kde10`, but for 1D kernels. 
+ `test_kde16`: Full example 2d kde with cross validation and adaptive densities. 
+ `test_kde17`: Full example 2d kde with cross validation and adaptive densities, but with weighted points. 
+ `test_point2d`:
+ `test_kernels`:
+ `test_kde_cppthread`:
+ `test_kde7`:
+ `test_kde6`:
+ `test_kde5`:
+ `test_kde1`:
+ `test_kde2`:
+ `test_kde3`:
+ `test_kde4`:
+ `test_kdtree3`:
+ `test_kdtree2`:
+ `test_kdtree1`:
+ `test_rectangle`:
+ `test_kde_attributes`:
+ `test_interval`:
+ `test_float_utils`:
+ `test_point_weight_attributes`:
+ `test_decorated_point`:
+ `test_point `:
+ `test_fft`:
+ `generate_bimodal_gauss`:
+ `prodkde2d_scan`:
+ `prodkde2d_fftcv`:
+ `prodadakde2d_scan`:

GPU tests
---------
+ `test_kde10_cuda`: Same test as `test_kde10`, but with the GPU. 
+ `test_kde11_cuda`: Same test as `test_kde11`, but with the GPU. 
+ `test_kde12_cuda`: Same test as `test_kde12`, but with the GPU. 
+ `test_kde13_cuda`: Same test as `test_kde13`, but with the GPU. 
+ `test_kde14_cuda`: Same test as `test_kde14`, but with the GPU. 
+ `test_kde16_cuda`: Same test as `test_kde16`, but with the GPU. 
+ `test_kde17_cuda`: Same test as `test_kde17`, but with the GPU. 
+ `test_point2d_cuda`:
+ `test_cukde0_cuda`: Basic usage of CudaDirectKde.
+ `test_cukde1_cuda`: Tests non kde evaluating methods (e.g. constructors, copy-control, etc.). Largely incomplete!
+ `test_cukde2_cuda`: Tests kde evaluating methods. 
