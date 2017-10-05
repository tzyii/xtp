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

void GPUNumericalIntegration::EvaluateFuncs(){
    // run the kernels baby

    // maybe this should be configurable??
    dim3 threadsperblock(512, 512); // this can change depending on hardware.
                                    // might need to be careful with this...

    // threadsperblock.x is the number of functions that are evaluated
    // simultaneously. threadsperblock.y is the number of grid points
    // at which each function is evaluated.

    // s
    dim3 numblocks(1 + _gpuAOB.sConts.size()/threadsperblock.x,
                   1 + _gpuGridBox.h_gridPoints.x.size()/threadsperblock.y);

    EvalSFuncs<<<numblocks, threadsperblock>>>(_gpuAOB.GetRawGPUArrs(),
                                               _gpuGridBox.GetRawGPUArrs(),
                                               _gpuFuncVals.sFuncVals);

    // p
    numblocks = dim3(1 + _gpuAOB.pConts.size()/threadsperblock.x,
                     1 + _gpuGridBox.h_gridPoints.x.size()/threadsperblock.y);

    EvalPFuncs<<<numblocks, threadsperblock>>>(_gpuAOB.GetRawGPUArrs(),
                                               _gpuGridBox.GetRawGPUArrs(),
                                               _gpuFuncVals.pFuncVals);

    // d
    numblocks = dim3(1 + _gpuAOB.dConts.size()/threadsperblock.x,
                     1 + _gpuGridBox.h_gridPoints.x.size()/threadsperblock.y);

    EvalDFuncs<<<numblocks, threadsperblock>>>(_gpuAOB.GetRawGPUArrs(),
                                               _gpuGridBox.GetRawGPUArrs(),
                                               _gpuFuncVals.dFuncVals);

    // f
    numblocks = dim3(1 + _gpuAOB.fConts.size()/threadsperblock.x,
                     1 + _gpuGridBox.h_gridPoints.x.size()/threadsperblock.y);

    EvalFFuncs<<<numblocks, threadsperblock>>>(_gpuAOB.GetRawGPUArrs(),
                                               _gpuGridBox.GetRawGPUArrs(),
                                               _gpuFuncVals.fFuncVals);


    // g
    numblocks = dim3(1 + _gpuAOB.gConts.size()/threadsperblock.x,
                     1 + _gpuGridBox.h_gridPoints.x.size()/threadsperblock.y);

    EvalGFuncs<<<numblocks, threadsperblock>>>(_gpuAOB.GetRawGPUArrs(),
                                               _gpuGridBox.GetRawGPUArrs(),
                                               _gpuFuncVals.gFuncVals);
}

}}}
