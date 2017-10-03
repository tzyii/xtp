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

    size_t baseNumFuncs = _gpuAOB.sConts.size()*_gpuGridBox.h_gridPoints.x.size();

    //allocate funcvalues
    d_sFuncVals = gpu_vector(1 * baseNumFuncs, 0); 
    d_pFuncVals = gpu_vector(3 * baseNumFuncs, 0); 
    d_dFuncVals = gpu_vector(5 * baseNumFuncs, 0); 
    d_fFuncVals = gpu_vector(7 * baseNumFuncs, 0); 
    d_gFuncVals = gpu_vector(9 * baseNumFuncs, 0);

    _gpuFuncVals.sFuncVals.array = thrust::raw_pointer_cast(&d_sFuncVals[0]);
    _gpuFuncVals.sFuncVals.arraySize = d_sFuncVals.size();

    _gpuFuncVals.pFuncVals.array = thrust::raw_pointer_cast(&d_pFuncVals[0]);
    _gpuFuncVals.pFuncVals.arraySize = d_pFuncVals.size();

    _gpuFuncVals.dFuncVals.array = thrust::raw_pointer_cast(&d_dFuncVals[0]);
    _gpuFuncVals.dFuncVals.arraySize = d_dFuncVals.size();

    _gpuFuncVals.fFuncVals.array = thrust::raw_pointer_cast(&d_fFuncVals[0]);
    _gpuFuncVals.fFuncVals.arraySize = d_fFuncVals.size();

    _gpuFuncVals.gFuncVals.array = thrust::raw_pointer_cast(&d_gFuncVals[0]);
    _gpuFuncVals.gFuncVals.arraySize = d_gFuncVals.size();
    
}

}

}}}
