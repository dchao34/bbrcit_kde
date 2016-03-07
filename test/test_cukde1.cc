#include <iostream>
#include <vector>

#include <DecoratedPoint.h>
#include <KdeTraits.h>
#include <Attributes/KdeAttributes.h>
#include <Kernels/EpanechnikovKernel.h>
#include <CudaDirectKde.h>
#include <FloatUtils.h>

namespace {
  using FloatType = float;
  using HostPointType = 
    bbrcit::DecoratedPoint<2, bbrcit::KdeAttributes<FloatType>, FloatType>;
  using DevicePointType = bbrcit::DevicePointTraits<2,FloatType>::Type;
  using KernelType = bbrcit::EpanechnikovKernel<2, FloatType>;
}

int main() {

  std::vector<HostPointType> ref1, query1;
  ref1.push_back({{1.0, 1.0}, {1.0,1.0,1.0}});
  query1.push_back({{2.0, 2.0}, {1.0,1.0,1.0}});
  query1.push_back({{2.0, 2.0}, {1.0,1.0,1.0}});
  bbrcit::CudaDirectKde<2,FloatType> cukde1(ref1, query1);

  std::vector<HostPointType> ref2, query2;
  ref2.push_back({{3.0, 3.0}, {1.0,1.0,1.0}});
  ref2.push_back({{4.0, 4.0}, {1.0,1.0,1.0}});
  ref2.push_back({{5.0, 5.0}, {1.0,1.0,1.0}});
  query2.push_back({{6.0, 6.0}, {1.0,1.0,1.0}});
  query2.push_back({{7.0, 7.0}, {1.0,1.0,1.0}});
  query2.push_back({{8.0, 8.0}, {1.0,1.0,1.0}});
  query2.push_back({{9.0, 9.0}, {1.0,1.0,1.0}});
  bbrcit::CudaDirectKde<2,FloatType> cukde2(ref2, query2);

  // swap
  std::cout << " + swap() (1): " << std::endl;
  swap(cukde1, cukde2);
  std::cout << "\t(" << cukde1.reference_size() << ", ";
  std::cout << cukde1.query_size() << "), ";
  std::cout << "(" << cukde2.reference_size() << ", ";
  std::cout << cukde2.query_size() << ") ";
  std::cout << "\t(c.f. (3, 4), (1, 2) )" << std::endl;
  swap(cukde1, cukde2);
  std::cout << "\t(" << cukde1.reference_size() << ", ";
  std::cout << cukde1.query_size() << "), ";
  std::cout << "(" << cukde2.reference_size() << ", ";
  std::cout << cukde2.query_size() << ") ";
  std::cout << "\t(c.f. (1, 2), (3, 4) )" << std::endl;

  std::cout << std::endl; 

  // copy-control

  std::cout << " + copy-constructor (1): " << std::endl;

  bbrcit::CudaDirectKde<2,FloatType> cukde3 = cukde1;
  std::cout << "\t(" << cukde3.reference_size() << ", ";
  std::cout << cukde3.query_size() << ") ";
  std::cout << "(c.f. (1, 2))" << std::endl;

  std::cout << " + move-constructor (1): " << std::endl;

  bbrcit::CudaDirectKde<2,FloatType> cukde4(std::move(cukde3));
  std::cout << "\t(" << cukde4.reference_size() << ", ";
  std::cout << cukde4.query_size() << ") ";
  std::cout << "(c.f. (1, 2))" << std::endl;

  std::cout << std::endl; 


  return 0;
}
