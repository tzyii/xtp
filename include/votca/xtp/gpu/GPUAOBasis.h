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

    thrust_R3N_vectors sPos;
    thrust_R3N_vectors pPos;
    thrust_R3N_vectors dPos;
    thrust_R3N_vectors fPos;
    thrust_R3N_vectors gPos;

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

    device_R3N_vectors d_sPos;
    device_R3N_vectors d_pPos;
    device_R3N_vectors d_dPos;
    device_R3N_vectors d_fPos;
    device_R3N_vectors d_gPos;

    gpu_vector d_expoFactors;


};

}}}
#endif // GPU_AO_BASIS_H
