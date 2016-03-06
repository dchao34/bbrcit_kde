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

int main() {

  // generate data
  size_t n_ref_pts = 8192;
  size_t n_query_pts = 8192;

  std::vector<HostPointType> ref_points;
  std::vector<HostPointType> query_points;

  std::default_random_engine e;
  std::normal_distribution<FloatType> d;

  FloatType point_weight = 
    bbrcit::ConstantTraits<FloatType>::one() / n_ref_pts;

  for (size_t i = 0; i < n_ref_pts; ++i) {
    ref_points.push_back({{d(e), d(e)}, {point_weight, 0.0f, 0.0f}});
  }

  for (size_t i = 0; i < n_query_pts; ++i) {
    query_points.push_back({{d(e), d(e)}, {point_weight, 0.0f, 0.0f}});
  }

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double, std::milli> elapsed;

  // cpu
  std::vector<DevicePointType> cpu_refs(n_ref_pts), cpu_query(n_query_pts);
  for (int i = 0; i < n_ref_pts; ++i) { cpu_refs[i] = ref_points[i]; }
  for (int i = 0; i < n_query_pts; ++i) { cpu_query[i] = query_points[i]; }

  KernelType kernel; kernel.set_bandwidth(1.0);
  std::vector<FloatType> cpu_results(n_query_pts);

  start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < n_query_pts; ++i) {
    FloatType sum = bbrcit::ConstantTraits<FloatType>::zero();
    for (size_t j = 0; j < n_ref_pts; ++j) {
      sum += cpu_refs[j].w() * 
             kernel.unnormalized_eval(cpu_refs[j], cpu_query[i]);
    }
    sum *= kernel.normalization();
    cpu_results[i] = sum;
  }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  std::cout << "cpu result (i = " << n_query_pts-1 << "): " << cpu_results[n_query_pts-1] << std::endl;
  std::cout << "cpu time: " << elapsed.count() << " ms" << std::endl;
  std::cout << std::endl;

  // gpu
  bbrcit::CudaDirectKde<2,FloatType, KernelType> cuda_kde(
      ref_points, query_points);

  cuda_kde.kernel().set_bandwidth(1.0);
  std::vector<FloatType> gpu_results;

  start = std::chrono::high_resolution_clock::now();
  cuda_kde.eval(0, n_ref_pts-1, 0, n_query_pts-1, gpu_results, 256);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;

  std::cout << "gpu result (i = " << n_query_pts-1 << "): " << gpu_results[n_query_pts-1] << std::endl;
  std::cout << "gpu time: " << elapsed.count() << " ms" << std::endl;
  std::cout << std::endl;

  std::cout << std::setprecision(12) << std::fixed;

  bool is_cpu_gpu_consistent = true;
  for (size_t i = 0; i < n_query_pts; ++i) {
    if (!bbrcit::approximately_equal(cpu_results[i], gpu_results[i],
                                     static_cast<FloatType>(1e-6f),
                                     static_cast<FloatType>(1e-8f))) {
      std::cout << "cpu vs gpu results disagree (i = " << i << "): ";
      std::cout << cpu_results[i] << " " << gpu_results[i] << std::endl;
      is_cpu_gpu_consistent = false;
    }
  }

  std::cout << "correctness test: ";
  if (is_cpu_gpu_consistent) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
  }

  return 0;
}
