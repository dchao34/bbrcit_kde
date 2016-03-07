#ifndef BBRCITKDE_CUDADIRECTKDE_H__
#define BBRCITKDE_CUDADIRECTKDE_H__

#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <vector>
#include <KdeTraits.h>
#include <Kernels/EpanechnikovKernel.h>

namespace bbrcit {

template<int D, typename FloatT, typename KernelT> class CudaDirectKde;

template<int D, typename FloatT, typename KernelT> 
void swap(CudaDirectKde<D,FloatT,KernelT>&, 
          CudaDirectKde<D,FloatT,KernelT>&);

// CudaDirectKde<> computes the kernel density estimate 
// using the direct algorithm, but on a CUDA capable GPU. 
//
// At present, each CudaDirectKde<> object performs evaluations
// on a fixed set of reference and query points. 
//
// (consider adding feature to allow point configurations?)
template<int D, 
         typename FloatT=double,
         typename KernelT=EpanechnikovKernel<D,FloatT>>
class CudaDirectKde {

  public:
    using KernelType = KernelT;
    using FloatType = FloatT;
    using DevicePointType = typename DevicePointTraits<D,FloatT>::Type;

    friend void swap<>(CudaDirectKde<D,FloatT,KernelT>&,
                       CudaDirectKde<D,FloatT,KernelT>&);

    // default constructor gives an empty estimator.
    CudaDirectKde();

    // copy-control
    CudaDirectKde(const CudaDirectKde<D,FloatT,KernelT>&);
    CudaDirectKde(CudaDirectKde<D,FloatT,KernelT>&&);
    CudaDirectKde& operator=(CudaDirectKde<D,FloatT,KernelT>);
    ~CudaDirectKde();


    // construct a density estimator over points in `ref_pts` to evaluate 
    // at points in `query_pts`. 
    //
    // Note: evaluations are specified using index ranges into the 
    //       given points. the ordering of points are assumed to be the 
    //       same as those specified at the time of construction. 
    template<typename HostPointT>
      CudaDirectKde(const std::vector<HostPointT> &ref_pts, 
                    const std::vector<HostPointT> &query_pts);



    // returns the number of reference points
    size_t reference_size() const; 

    // returns the number of query points
    size_t query_size() const; 

    // return a reference to the kernel. allow users to configure it directly
    const KernelType& kernel() const;
    KernelType& kernel();



    // evaluate the density contributions due to reference points in the 
    // closed interval [r_i, r_j] for query points in the 
    // closed interval [q_i, q_j]. 
    //
    // the evaluation result for query point q_k in [q_i, q_j]
    // is stored in result[q_k-q_i]. 
    //
    // use `blocksize` to tune performance
    void eval(size_t r_i, size_t r_j, size_t q_i, size_t q_j, 
              std::vector<FloatT> &results, size_t block_size=128);



  private:

    // CudaDirectKde<>'s internal state: 
    //
    // + dev_ref_points_, n_ref_points_: pointer to device
    //       reference points and the number of points. 
    //
    // + dev_query_points_, n_query_points_: pointer to device
    //       query points and the number of points. 
    //
    // + dev_results_: pointer to staging area for computation 
    //       results. the number of elements is n_query_points_. 
    //
    // + kernel_: density kernel on the host. users configure this 
    //       directly; it is copied to the device everytime an evaluation 
    //       is requested.

    size_t n_ref_points_;
    size_t n_query_points_;

    DevicePointType *dev_ref_points_;
    DevicePointType *dev_query_points_;
    FloatType *dev_results_;

    KernelType kernel_;

    // constructor helper functions
    template<typename HostPointT>
      void alloc_copy_host_dev_points(
          DevicePointType**, const std::vector<HostPointT>&);

};


template<typename PointT, typename FloatT, typename KernelT>
__global__ void eval_kernel(
    size_t i_r, size_t j_r, size_t i_q, size_t j_q,
    PointT *ref_points, PointT *query_points, 
    FloatT *results, 
    KernelT *kernel) {

  extern __shared__ PointT ref_cache[];

  // cache the density kernel into a register
  KernelT kern = *kernel;

  // cache the query point assigned to this thread into a register
  size_t q_idx = i_q + blockIdx.x * blockDim.x + threadIdx.x;
  PointT query; 
  if (q_idx <= j_q) { query = query_points[q_idx]; }

  // main kde computation using the tiling algorithm. see 
  // chapter 31 of GPU Gems 3 for details
  FloatT sum = ConstantTraits<FloatT>::zero();
  for (int tile_start = i_r; tile_start <= j_r; tile_start+=blockDim.x) {

    // phase 1: cache reference points belonging to this tile into 
    //          shared memory. each thread in the block is responsible
    //          for fetching one point. 
    size_t r_idx = tile_start + threadIdx.x;
    if (r_idx <= j_r) { ref_cache[threadIdx.x] = ref_points[r_idx]; }
    __syncthreads();

    // phase 2: sum over contributions from the cached reference points
    if (q_idx <= j_q) {
      for (size_t i = 0; i < blockDim.x && i+tile_start <= j_r; ++i) {
        sum += ref_cache[i].w() * kern.unnormalized_eval(query, ref_cache[i]);
      }
    }
    __syncthreads();
  }

  // post process and write the result back to global memory. 
  if (q_idx <= j_q) {
    sum *= kern.normalization();
    results[q_idx - i_q] = sum;
  }

  return;

}


template<int D, typename FT, typename KT>
void CudaDirectKde<D,FT,KT>::eval(
    size_t i_r, size_t j_r, size_t i_q, size_t j_q, 
    std::vector<FloatType> &results, size_t block_size) {

  if (j_r >= n_ref_points_ || j_q >= n_query_points_) { 
    throw std::out_of_range(
        "CudaDirectKde<>: eval(): j_r and j_q should be at most "
        "the number of reference and query points. ");
  }

  if (i_r > j_r || i_q > j_q) { 
    throw std::invalid_argument(
        "CudaDirectKde<>: eval(): must have i_r <= j_r and i_q <= j_q. ");
  }

  // explicitly copy the density kernel to the device. this is because
  // users are allowed to configure kernels differently between 
  // separate calls
  KernelType *dev_kernel;
  cudaMalloc((void**)&dev_kernel, sizeof(KernelType));
  cudaMemcpy(dev_kernel, &kernel_, 
             sizeof(KernelType), cudaMemcpyHostToDevice);

  // configure and launch the CUDA kernel
  size_t n_queries = j_q-i_q+1;
  size_t grid_size = (n_queries+block_size-1)/block_size;
  size_t shared_memsize = block_size * sizeof(DevicePointType);

  eval_kernel<<<grid_size, block_size, shared_memsize>>>(
      i_r, j_r, i_q, j_q, 
      dev_ref_points_, dev_query_points_,
      dev_results_, dev_kernel);

  // copy the results to the host
  if (results.size() < n_queries) { results.resize(n_queries); }
  cudaMemcpy(&results[0], dev_results_, 
             sizeof(FloatType)*n_queries, cudaMemcpyDeviceToHost);

  cudaFree(dev_kernel);

  return; 
}


template<int D, typename FT, typename KT>
void swap(CudaDirectKde<D,FT,KT> &lhs, CudaDirectKde<D,FT,KT> &rhs) {
  using std::swap;

  swap(lhs.n_ref_points_, rhs.n_ref_points_);
  swap(lhs.dev_ref_points_, rhs.dev_ref_points_);

  swap(lhs.n_query_points_, rhs.n_query_points_);
  swap(lhs.dev_query_points_, rhs.dev_query_points_);

  swap(lhs.dev_results_, rhs.dev_results_);

  swap(lhs.kernel_, rhs.kernel_);
}


template<int D, typename FT, typename KT>
inline CudaDirectKde<D,FT,KT>& CudaDirectKde<D,FT,KT>::operator=(
    CudaDirectKde<D,FT,KT> rhs) {
  swap(*this, rhs); return *this;
}


template<int D, typename FT, typename KT>
CudaDirectKde<D,FT,KT>::CudaDirectKde(
  const CudaDirectKde<D,FT,KT> &rhs) : 
  n_ref_points_(rhs.n_ref_points_), 
  n_query_points_(rhs.n_query_points_),
  kernel_(rhs.kernel_) {

  cudaMalloc((void**)&dev_ref_points_, 
             sizeof(DevicePointType)*n_ref_points_);
  cudaMemcpy(dev_ref_points_, rhs.dev_ref_points_,
             sizeof(DevicePointType)*n_ref_points_, 
             cudaMemcpyDeviceToDevice);

  cudaMalloc((void**)&dev_query_points_, 
             sizeof(DevicePointType)*n_query_points_);
  cudaMemcpy(dev_query_points_, rhs.dev_query_points_,
             sizeof(DevicePointType)*n_query_points_, 
             cudaMemcpyDeviceToDevice);

}

template<int D, typename FT, typename KT>
CudaDirectKde<D,FT,KT>::CudaDirectKde(CudaDirectKde<D,FT,KT> &&rhs) 
  : n_ref_points_(std::move(rhs.n_ref_points_)), 
    n_query_points_(std::move(rhs.n_query_points_)),
    kernel_(std::move(rhs.kernel_)) {

  dev_ref_points_ = rhs.dev_ref_points_;
  rhs.dev_ref_points_ = nullptr;

  dev_query_points_ = rhs.dev_query_points_;
  rhs.dev_query_points_ = nullptr;
}

template<int D, typename FT, typename KT>
inline size_t CudaDirectKde<D,FT,KT>::reference_size() const {
  return n_ref_points_;
}

template<int D, typename FT, typename KT>
inline size_t CudaDirectKde<D,FT,KT>::query_size() const {
  return n_query_points_;
}

template<int D, typename FT, typename KT>
inline const typename CudaDirectKde<D,FT,KT>::KernelType&
CudaDirectKde<D,FT,KT>::kernel() const {
  return kernel_;
}

template<int D, typename FT, typename KT>
inline typename CudaDirectKde<D,FT,KT>::KernelType&
CudaDirectKde<D,FT,KT>::kernel() {
  return const_cast<KernelType&>(
           static_cast<const CudaDirectKde<D,FT,KT>&>(*this).kernel()
      );
}

template<int D, typename FT, typename KT>
CudaDirectKde<D,FT,KT>::CudaDirectKde() :
  n_ref_points_(0),
  n_query_points_(0),
  dev_ref_points_(nullptr),
  dev_query_points_(nullptr), 
  dev_results_(nullptr), 
  kernel_() {}


// Input: 
//   + host_pts: vector of HostPointT's on the host.  
//
// Output: 
//   + dev_p: pointer to an array of DevicePointType's where 
//            dev_p[i] is obtained by first copy-constructing host_pt[i]
//            and subsequently transferred to the device. 
template<int D, typename FT, typename KT>
  template<typename HostPointT>
void CudaDirectKde<D,FT,KT>::alloc_copy_host_dev_points(
    DevicePointType **dev_p, 
    const std::vector<HostPointT> &host_pts) {

  // allocate array to copy-construct DevicePointType's on the host
  size_t n_points = host_pts.size();
  DevicePointType *host_points_p = new DevicePointType[n_points];

  // copy construct
  for (size_t i = 0; i < n_points; ++i) { 
    host_points_p[i] = host_pts[i]; 
  }

  // allocate memory on the device and then transfer the copy-constructed
  // DevicePointType's on the host over
  cudaMalloc((void**)dev_p, sizeof(DevicePointType)*n_points);
  cudaMemcpy(*dev_p, host_points_p,
             sizeof(DevicePointType)*n_points, cudaMemcpyHostToDevice);

  delete[] host_points_p;

}

template<int D, typename FT, typename KT>
  template<typename HostPointT>
CudaDirectKde<D,FT,KT>::CudaDirectKde(
    const std::vector<HostPointT> &ref_pts,
    const std::vector<HostPointT> &query_pts) : kernel_() {

  // allocate device reference points
  n_ref_points_ = ref_pts.size();
  alloc_copy_host_dev_points(&dev_ref_points_, ref_pts);

  // allocate device query points
  n_query_points_ = query_pts.size();
  alloc_copy_host_dev_points(&dev_query_points_, query_pts);

  // allocate device memory to stage computed results
  cudaMalloc((void**)&dev_results_, sizeof(FloatType) * n_query_points_);

  cudaDeviceSynchronize();

}

template<int D, typename FT, typename KT>
CudaDirectKde<D,FT,KT>::~CudaDirectKde() {
  cudaFree(dev_ref_points_);
  cudaFree(dev_query_points_);
  cudaFree(dev_results_);
}

}

#endif
