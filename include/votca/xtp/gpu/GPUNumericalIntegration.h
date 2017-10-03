#ifndef GPU_NUMERICAL_INTEGRATION
#define GPU_NUMERICAL_INTEGRATION

#include <cuda.h>
#include <vector>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/gridbox.h>
#include <votca/xtp/gpu/GPUAOBasis.h>
#include <votca/xtp/gpu/GPUGridBox.h>

namespace votca { namespace xtp { namespace gpu {

   class GPUNumericalIntegration{
   public:
       GPUNumericalIntegration(){}
       ~GPUNumericalIntegration(){}
       GPUNumericalIntegration(const AOBasis& aob, const std::vector<GridBox>& gbs);
   private:
       int _device;
       GPUAOBasis _gpuAOB;
       GPUGridBox _gpuGridBox;

       gpu_vector d_sFuncVals;
       gpu_vector d_pFuncVals;
       gpu_vector d_dFuncVals;
       gpu_vector d_fFuncVals;
       gpu_vector d_gFuncVals;

       GPUFuncVals _gpuFuncVals;
   };

}}}


#endif // GPU_NUMERICAL_INTEGRATION_CUH
