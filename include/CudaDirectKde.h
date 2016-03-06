#ifndef BBRCITKDE_CUDADIRECTKDE_H__
#define BBRCITKDE_CUDADIRECTKDE_H__

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <KdeTraits.h>
#include <Kernels/EpanechnikovKernel.h>

namespace bbrcit {

template<int D, 
         typename FloatT=double,
         typename KernelT=EpanechnikovKernel<D,FloatT>>
class CudaDirectKde {

  public:
    using KernelType = KernelT;
    using FloatType = FloatT;
    using DevicePointType = typename DevicePointTraits<D,FloatT>::Type;

    CudaDirectKde();
    ~CudaDirectKde();

    template<typename HostPointT>
      CudaDirectKde(const std::vector<HostPointT> &ref_pts, 
                    const std::vector<HostPointT> &query_pts);

    const KernelType& kernel() const;
    KernelType& kernel();

    void eval(size_t r_i, size_t r_j, size_t q_i, size_t q_j, 
              std::vector<FloatT> &results, size_t block_size);

  private:
    size_t n_ref_points_;
    size_t n_query_points_;
    DevicePointType *dev_ref_points_;
    DevicePointType *dev_query_points_;
    FloatType *dev_results_;

    KernelType kernel_;

    template<typename HostPointT>
      void alloc_copy_host_dev_points(
          DevicePointType**, const std::vector<HostPointT>&);

};


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

template<typename PointT, typename FloatT, typename KernelT>
__global__ void eval_kernel(
    size_t i_r, size_t j_r, size_t i_q, size_t j_q,
    PointT *ref_points, 
    PointT *query_points, 
    FloatT *results, 
    KernelT *kernel) {

  extern __shared__ PointT ref_cache[];

  KernelT kern = *kernel;

  // load query point
  size_t q_idx = i_q + blockIdx.x * blockDim.x + threadIdx.x;
  PointT query; 
  if (q_idx <= j_q) { query = query_points[q_idx]; }

  FloatT sum = ConstantTraits<FloatT>::zero();
  for (int tile_start = i_r; tile_start <= j_r; tile_start+=blockDim.x) {

    // phase 1: load reference points
    size_t r_idx = tile_start + threadIdx.x;
    if (r_idx <= j_r) { ref_cache[threadIdx.x] = ref_points[r_idx]; }
    __syncthreads();

    // phase 2: advance sums
    if (q_idx <= j_q) {
      for (size_t i = 0; i < blockDim.x && i+tile_start <= j_r; ++i) {
        sum += ref_cache[i].w() * kern.unnormalized_eval(query, ref_cache[i]);
      }
    }
    __syncthreads();
  }

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

  KernelType *dev_kernel;
  cudaMalloc((void**)&dev_kernel, sizeof(KernelType));
  cudaMemcpy(dev_kernel, &kernel_, 
             sizeof(KernelType), cudaMemcpyHostToDevice);

  size_t n_queries = j_q-i_q+1;
  size_t grid_size = (n_queries+block_size-1)/block_size;
  size_t shared_memsize = block_size * sizeof(DevicePointType);

  eval_kernel<<<grid_size, block_size, shared_memsize>>>(
      i_r, j_r, i_q, j_q, 
      dev_ref_points_, dev_query_points_,
      dev_results_, dev_kernel);

  if (results.size() < n_queries) { results.resize(n_queries); }
  cudaMemcpy(&results[0], dev_results_, 
             sizeof(FloatType)*n_queries, cudaMemcpyDeviceToHost);

  cudaFree(dev_kernel);

  return; 
}

template<int D, typename FT, typename KT>
CudaDirectKde<D,FT,KT>::CudaDirectKde() :
  n_ref_points_(0),
  n_query_points_(0),
  dev_ref_points_(nullptr),
  dev_query_points_(nullptr), 
  dev_results_(nullptr), 
  kernel_() {}

template<int D, typename FT, typename KT>
  template<typename HostPointT>
void CudaDirectKde<D,FT,KT>::alloc_copy_host_dev_points(
    DevicePointType **dev_p, 
    const std::vector<HostPointT> &host_pts) {

  size_t n_points = host_pts.size();
  
  DevicePointType *host_points_p = new DevicePointType[n_points];

  // copy host point type into device point type on the host
  for (size_t i = 0; i < n_points; ++i) { 
    host_points_p[i] = host_pts[i]; 
  }

  // copy device point type from host to device
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

  // reference points
  n_ref_points_ = ref_pts.size();
  alloc_copy_host_dev_points(&dev_ref_points_, ref_pts);

  // query points
  n_query_points_ = query_pts.size();
  alloc_copy_host_dev_points(&dev_query_points_, query_pts);

  // memory to stage results on the device
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
