#ifndef GPU_AO_BASIS_H
#define GPU_AO_BASIS_H
#include <vector>
#include <votca/xtp/gpu/GPUTools.cuh>
#include <votca/xtp/aoshell.h>
#include <votca/xtp/aobasis.h>

namespace votca { namespace xtp { namespace gpu {
class GPUAOBasis{
public:

    GPUAOBasis();
    ~GPUAOBasis();

    thrust_vector sConts;
    thrust_vector pConts;
    thrust_vector dConts;
    thrust_vector fConts;
    thrust_vector gConts;
    thrust_vector alphas;
    thrust_vector powFactors;

    GPUAOBasis(const AOBasis& aob);

    gpuAOArrs GetRawGPUArrs();

private:
    gpuAOArrs rawGpuArrs;

    size_t padding = 16;

    gpu_vector d_sConts;
    gpu_vector d_pConts;
    gpu_vector d_dConts;
    gpu_vector d_fConts;
    gpu_vector d_gConts;
    gpu_vector d_alphas;
    gpu_vector d_powFactors;

    gpu_vector d_expoFactors;

    gpu_vector d_sFuncVals;
    gpu_vector d_pFuncVals;
    gpu_vector d_dFuncVals;
    gpu_vector d_fFuncVals;
    gpu_vector d_gFuncVals;
};

}}}
#endif // GPU_AO_BASIS_H
