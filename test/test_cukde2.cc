#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

#include <DecoratedPoint.h>
#include <KdeTraits.h>
#include <Attributes/KdeAttributes.h>
#include <Kernels/EpanechnikovKernel.h>
#include <CudaDirectKde.h>
#include <FloatUtils.h>

namespace {
  using FloatType = float;
  using HostPointType = bbrcit::DecoratedPoint<2, bbrcit::KdeAttributes<FloatType>, FloatType>;
  using DevicePointType = bbrcit::DevicePointTraits<2,FloatType>::Type;
  using KernelType = bbrcit::EpanechnikovKernel<2, FloatType>;
}

void generate_data(
    size_t n_ref_pts, size_t n_query_pts, 
    std::vector<HostPointType> &ref_pts, 
    std::vector<HostPointType> &query_pts) {

  ref_pts.clear(); query_pts.clear();

  std::default_random_engine e;
  std::normal_distribution<FloatType> d;

  FloatType point_weight = 
    bbrcit::ConstantTraits<FloatType>::one() / n_ref_pts;

  for (size_t i = 0; i < n_ref_pts; ++i) {
    ref_pts.push_back({{d(e), d(e)}, {point_weight, 0.0f, 0.0f}});
  }

  for (size_t i = 0; i < n_query_pts; ++i) {
    query_pts.push_back({{d(e), d(e)}, {point_weight, 0.0f, 0.0f}});
  }
}

void cpu_kde(size_t i_r, size_t j_r, size_t i_q, size_t j_q,
             const std::vector<DevicePointType> &ref_pts, 
             const std::vector<DevicePointType> &query_pts, 
             std::vector<FloatType> &results, 
             const KernelType &kernel) {

  for (size_t i = i_q; i <= j_q; ++i) {
    FloatType sum = bbrcit::ConstantTraits<FloatType>::zero();
    for (size_t j = i_r; j <= j_r; ++j) {
      sum += ref_pts[j].w() * 
             kernel.unnormalized_eval(ref_pts[j], query_pts[i]);
    }
    sum *= kernel.normalization();
    results[i-i_q] = sum;
  }

}

bool is_cpu_gpu_consistent(
    size_t n_queries,
    const std::vector<FloatType> &cpu_results,
    const std::vector<FloatType> &gpu_results,
    FloatType rel_tol, FloatType abs_tol) {

  bool is_consistent = true;
  for (size_t i = 0; i < n_queries; ++i) {
    if (!bbrcit::approximately_equal(cpu_results[i], gpu_results[i],
                                     static_cast<FloatType>(rel_tol),
                                     static_cast<FloatType>(abs_tol))) {
      std::cout << "cpu vs gpu results disagree (i = " << i << "): ";
      std::cout << cpu_results[i] << " " << gpu_results[i] << std::endl;
      is_consistent = false;
    }
  }

  std::cout << "cpu results (i=0): ";
  std::cout << cpu_results[0] << std::endl;
  std::cout << "gpu results (i=0): ";
  std::cout << gpu_results[0] << std::endl;
  std::cout << std::endl;
  std::cout << "cpu results (i=" << n_queries-1 << "): ";
  std::cout << cpu_results[n_queries-1] << std::endl;
  std::cout << "gpu results (i=" << n_queries-1 << "): ";
  std::cout << gpu_results[n_queries-1] << std::endl;
  std::cout << std::endl;

  std::cout << "correctness status: ";
  if (is_consistent) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
  }


  return is_consistent;
}

void eval_test(
    size_t i_r, size_t j_r, size_t i_q, size_t j_q, 
    std::vector<DevicePointType> &cpu_refs,
    std::vector<DevicePointType> &cpu_queries,
    const KernelType &cpu_kernel,
    std::vector<FloatType> &cpu_results,
    bbrcit::CudaDirectKde<2,FloatType> &cuda_kde, 
    std::vector<FloatType> &gpu_results) {


  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double, std::milli> elapsed;

  // cpu
  start = std::chrono::high_resolution_clock::now();
  cpu_kde(i_r, j_r, i_q, j_q, 
          cpu_refs, cpu_queries, cpu_results, cpu_kernel);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  FloatType cpu_time = elapsed.count();
  std::cout << "cpu time: " << cpu_time << " ms" << std::endl;

  // gpu
  start = std::chrono::high_resolution_clock::now();
  cuda_kde.eval(i_r, j_r, i_q, j_q, gpu_results);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  FloatType gpu_time = elapsed.count();
  std::cout << "gpu time: " << gpu_time << " ms" << std::endl;

  std::cout << "speed up: " << cpu_time / gpu_time << "x " << std::endl;


}

int main() {

  size_t n_ref_pts = 1024; size_t n_query_pts = 1024;

  // generate data
  std::vector<HostPointType> ref_points, query_points;
  generate_data(n_ref_pts, n_query_pts, ref_points, query_points);

  // setup cpu kde
  std::vector<DevicePointType> cpu_refs(n_ref_pts), cpu_queries(n_query_pts);
  for (int i = 0; i < n_ref_pts; ++i) { cpu_refs[i] = ref_points[i]; }
  for (int i = 0; i < n_query_pts; ++i) { cpu_queries[i] = query_points[i]; }
  KernelType kernel; kernel.set_bandwidth(1.0);
  std::vector<FloatType> cpu_results(n_query_pts);

  // setup gpu kde
  bbrcit::CudaDirectKde<2,FloatType, KernelType> cuda_kde(
      ref_points, query_points);
  cuda_kde.kernel().set_bandwidth(1.0);
  std::vector<FloatType> gpu_results(n_query_pts);

  // tests
  // -----

  size_t i_r, j_r, i_q, j_q;
  size_t n_queries;
  
  std::cout << "+ all queries / all refs" << std::endl;
  std::cout << std::endl; 

  i_r = 0, j_r = n_ref_pts-1, i_q = 0, j_q = n_query_pts-1;
  n_queries = j_q-i_q+1;

  eval_test(i_r, j_r, i_q, j_q,
            cpu_refs, cpu_queries, kernel, 
            cpu_results, cuda_kde, gpu_results);
  
  std::cout << std::endl;

  is_cpu_gpu_consistent(n_queries, 
                        cpu_results, gpu_results, 
                        1e-6, 1e-8);

  std::cout << std::endl;
  std::cout << std::endl;
  
  std::cout << "+ first half queries / all refs" << std::endl;
  std::cout << std::endl; 

  i_r = 0, j_r = n_ref_pts-1, i_q = 0, j_q = n_query_pts/2-1;
  n_queries = j_q-i_q+1;

  eval_test(i_r, j_r, i_q, j_q,
            cpu_refs, cpu_queries, kernel, 
            cpu_results, cuda_kde, gpu_results);
  
  std::cout << std::endl;

  is_cpu_gpu_consistent(n_queries, 
                        cpu_results, gpu_results, 
                        1e-6, 1e-8);

  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "+ second half queries / all refs" << std::endl;
  std::cout << std::endl; 

  i_r = 0, j_r = n_ref_pts-1, i_q = n_query_pts/2, j_q = n_query_pts-1;
  n_queries = j_q-i_q+1;

  eval_test(i_r, j_r, i_q, j_q,
            cpu_refs, cpu_queries, kernel, 
            cpu_results, cuda_kde, gpu_results);
  
  std::cout << std::endl;

  is_cpu_gpu_consistent(n_queries, 
                        cpu_results, gpu_results, 
                        1e-6, 1e-8);

  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "+ all queries / first half refs" << std::endl;
  std::cout << std::endl; 

  i_r = 0, j_r = n_ref_pts/2-1, i_q = 0, j_q = n_query_pts-1;
  n_queries = j_q-i_q+1;

  eval_test(i_r, j_r, i_q, j_q,
            cpu_refs, cpu_queries, kernel, 
            cpu_results, cuda_kde, gpu_results);
  
  std::cout << std::endl;

  is_cpu_gpu_consistent(n_queries, 
                        cpu_results, gpu_results, 
                        1e-6, 1e-8);

  std::cout << std::endl;
  std::cout << std::endl;


  std::cout << "+ all queries / second half refs" << std::endl;
  std::cout << std::endl; 

  i_r = n_ref_pts/2, j_r = n_ref_pts-1, i_q = 0, j_q = n_query_pts-1;
  n_queries = j_q-i_q+1;

  eval_test(i_r, j_r, i_q, j_q,
            cpu_refs, cpu_queries, kernel, 
            cpu_results, cuda_kde, gpu_results);
  
  std::cout << std::endl;

  is_cpu_gpu_consistent(n_queries, 
                        cpu_results, gpu_results, 
                        1e-6, 1e-8);

  std::cout << std::endl;
  std::cout << std::endl;

  return 0;
}
