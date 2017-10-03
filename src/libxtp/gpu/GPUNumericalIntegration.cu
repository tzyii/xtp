#include <algorithm>
#include <string>
#include <votca/xtp/aoshell.h>
#include <votca/xtp/gridbox.h>
#include <thrust/host_vector.h>
#include <votca/xtp/gpu/GPUNumericalIntegration.h>
#include <votca/xtp/gpu/GPUCTypes.cuh>
#include "kernels/IntegrationKernels.cuh"

namespace votca { namespace xtp { namespace gpu {

GPUNumericalIntegration::GPUNumericalIntegration(const AOBasis& aob,
                                                 const std::vector<GridBox>& gbs){
    _device=0;
    CUDA_API_CALL(cudaSetDevice(_device), true); // add some code here later to intelligently
                                                 // select a gpu
    
    // Later we need to also dissallow this gpu from begin used
    // by other processes (optionally?)...
    _gpuAOB = GPUAOBasis(aob);
    _gpuGridBox = GPUGridBox(gbs);
}

}}}
