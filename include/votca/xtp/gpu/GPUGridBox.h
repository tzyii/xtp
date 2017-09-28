#ifndef GPU_GRID_BOX_CUH
#define GPU_GRID_BOX_CUH

#include <votca/xtp/gpu/GPUTools.cuh>
#include <votca/xtp/gpu/GPUCTypes.cuh>
#include <votca/xtp/gridbox.h>

namespace votca { namespace xtp { namespace gpu {

class GPUGridBox{
public:
    thrust_R3N_vectors h_gridPoints;
    thrust_vector h_weights;

    GPUGridBox();
    ~GPUGridBox();

    GPUGridBox(const std::vector<GridBox>& gridBoxes);

private:
    R3Nptrs gridPoints;
    device_R3N_vectors d_gridPoints;
    thrust_vector d_weights;
};


}}}
#endif // GPU_GRID_BOX_CUH
