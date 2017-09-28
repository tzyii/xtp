#include <votca/xtp/gpu/GPUGridBox.h>

namespace votca { namespace xtp { namespace gpu {

GPUGridBox::GPUGridBox(){}
GPUGridBox::~GPUGridBox(){}

GPUGridBox::GPUGridBox(const std::vector<GridBox>& gridBoxes){
        
    for (auto gridBox: gridBoxes){
        for (auto p: gridBox.getGridPoints()){
            h_gridPoints.x.push_back(p[0]); 
            h_gridPoints.y.push_back(p[1]); 
            h_gridPoints.z.push_back(p[2]); 
        }
            
        for (const double w: gridBox.getGridWeights()){
            h_weights.push_back(w);
        }
    }

    d_gridPoints = device_R3N_vectors(h_gridPoints);
    d_weights = gpu_vector(h_weights);
}
}}}
