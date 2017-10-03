#include <votca/xtp/gpu/GPUCTypes.cuh>
#include "VecFuncs.cuh"
namespace votca { namespace xtp { namespace gpu {

__global__ void EvalSFuncs(const gpuAOArrs aob, const R3Nptrs gridPoints,
                           gpuArr _sFuncVal){

    size_t fi = blockDim.x*blockIdx.x + threadIdx.x;
    size_t gi = blockDim.y*blockIdx.y + threadIdx.y;

    if (fi < aob.alphas.arraySize & gi < gridPoints.n){
        double3 center = R3NIndex(gridPoints, gi) - R3NIndex(aob.sFuncPos, fi);
        double powFactor = aob.powFactors.array[fi];
        double alpha = aob.alphas.array[fi];
        double distsq = dot(center, center);
        
        double expoFactor = powFactor * exp(-alpha*distsq); 
            
        _sFuncVal.array[fi*gridPoints.n + gi + 0] = expoFactor * aob.sConts.array[fi];
    }
    
}

__global__ void EvalPFuncs(const gpuAOArrs aob, const R3Nptrs gridPoints,
                           gpuArr _pFuncVal){
    
    size_t fi = blockDim.x*blockIdx.x + threadIdx.x;
    size_t gi = blockDim.y*blockIdx.y + threadIdx.y;

    if (fi < aob.alphas.arraySize & gi < gridPoints.n){
        double3 center = R3NIndex(gridPoints, gi) - R3NIndex(aob.sFuncPos, fi);
        double powFactor = aob.powFactors.array[fi];
        double alpha = aob.alphas.array[fi];
        double distsq = dot(center, center);
        
        double expoFactor = powFactor * exp(-alpha*distsq);
        
        double factor = 2*sqrt(aob.alphas.array[fi])*aob.pConts.array[fi];
        
        _pFuncVal.array[fi*gridPoints.n + gi + 0] = factor * center.z * expoFactor; 
        _pFuncVal.array[fi*gridPoints.n + gi + 1] = factor * center.y * expoFactor; 
        _pFuncVal.array[fi*gridPoints.n + gi + 2] = factor * center.x * expoFactor; 
    }
}

__global__ void EvalDFuncs(const gpuAOArrs aob, const R3Nptrs gridPoints,
                           gpuArr _dFuncVal){
    
    size_t fi = blockDim.x*blockIdx.x + threadIdx.x;
    size_t gi = blockDim.y*blockIdx.y + threadIdx.y;

    if (fi < aob.alphas.arraySize & gi < gridPoints.n){
        double3 center = R3NIndex(gridPoints, gi) - R3NIndex(aob.sFuncPos, fi);
        double powFactor = aob.powFactors.array[fi];
        double alpha = aob.alphas.array[fi];
        double distsq = dot(center, center);
        
        double expoFactor = powFactor * exp(-alpha*distsq); 

        double f = 2.*aob.alphas.array[fi]*aob.dConts.array[fi];
        
        _dFuncVal.array[fi*gridPoints.n + gi + 0] =
            f*rsqrt(3.) * (3.*center.z*center.z - dot(center, center)) * expoFactor;

        _dFuncVal.array[fi*gridPoints.n + gi + 1] =
            2.*f * (center.y * center.z) * expoFactor;

        _dFuncVal.array[fi*gridPoints.n + gi + 2] =
            2.*f * (center.x * center.z) * expoFactor;

        _dFuncVal.array[fi*gridPoints.n + gi + 3] =
            2.*f * (center.x * center.y) * expoFactor;

        _dFuncVal.array[fi*gridPoints.n + gi + 4] =
            f * (center.x*center.x - center.y*center.y) * expoFactor; 
    }
}

__global__ void EvalFFuncs(const gpuAOArrs aob, const R3Nptrs gridPoints,
                           gpuArr _fFuncVal){
    
    size_t fi = blockDim.x*blockIdx.x + threadIdx.x;
    size_t gi = blockDim.y*blockIdx.y + threadIdx.y;

    if (fi < aob.alphas.arraySize & gi < gridPoints.n){
        double3 center = R3NIndex(gridPoints, gi) - R3NIndex(aob.sFuncPos, fi);
        double powFactor = aob.powFactors.array[fi];
        double alpha = aob.alphas.array[fi];
        double distsq = dot(center, center);
        
        double expoFactor = powFactor * exp(-alpha*distsq); 
        
        double f = 2.*sqrt(alpha*alpha*alpha)*aob.fConts.array[fi];
        double f1 = f*2*rsqrt(15.);
        double f2 = f*sqrt(2.)*rsqrt(5.);
        double f3 = f*sqrt(2.)*rsqrt(3.);

        double cx2 = center.x*center.x;
        double cy2 = center.y*center.y;
        double cz2 = center.z*center.z;

        _fFuncVal.array[fi*gridPoints.n + gi + 0] =
            f1 * center.z * (5.*cz2 - 3.*distsq) * expoFactor; 

        _fFuncVal.array[fi*gridPoints.n + gi + 1] = 
            f2 * center.y * (5.*cz2 - distsq) * expoFactor;

        _fFuncVal.array[fi*gridPoints.n + gi + 2] =
            f2 * center.x * (5.*cz2 - distsq) * expoFactor;

        _fFuncVal.array[fi*gridPoints.n + gi + 3] =
            4.*f * center.x * center.y * center.z * expoFactor;

        _fFuncVal.array[fi*gridPoints.n + gi + 4] =
            2.*f * center.z * (cx2 - cy2) * expoFactor;
        
        _fFuncVal.array[fi*gridPoints.n + gi + 5] =
            f3 * center.y * (3.*cx2 - cy2) * expoFactor;

        _fFuncVal.array[fi*gridPoints.n + gi + 6] =
            f3 * center.x * (cx2 - 3.*cy2) * expoFactor; 

    }
}

__global__ void EvalGFuncs(const gpuAOArrs aob, const R3Nptrs gridPoints,
                           gpuArr _gFuncVal){
    
    size_t fi = blockDim.x*blockIdx.x + threadIdx.x;
    size_t gi = blockDim.y*blockIdx.y + threadIdx.y;

    if (fi < aob.alphas.arraySize & gi < gridPoints.n){
        double3 center = R3NIndex(gridPoints, gi) - R3NIndex(aob.sFuncPos, fi);
        double powFactor = aob.powFactors.array[fi];
        double alpha = aob.alphas.array[fi];
        double distsq = dot(center, center);
        
        double expoFactor = powFactor * exp(-alpha*distsq); 

        double f = 2.*rsqrt(3.)*alpha*alpha*aob.gConts.array[fi];
        double f1 = f*rsqrt(35.);
        double f2 = f*4.*rsqrt(14.);
        double f3 = f*2*rsqrt(7.);
        double f4 = f*2.*sqrt(2.);

        double cx2 = center.x*center.x;
        double cxcy = center.x*center.y;
        double cxcz = center.x*center.z;
        double cy2 = center.y*center.y;
        double cycz = center.y*center.z;
        double cz2 = center.z*center.z;

        _gFuncVal.array[fi*gridPoints.n + gi + 0] =
            f1 * (35.*cz2*cz2 - 30.*cz2*distsq + 3.*distsq*distsq)*expoFactor;

        _gFuncVal.array[fi*gridPoints.n + gi + 1] =
            f2 * cycz * (7.*cz2 - 3.*distsq) * expoFactor;
        
        _gFuncVal.array[fi*gridPoints.n + gi + 2] =
            f2 * cxcz * (7.*cz2 - 3.*distsq) * expoFactor; 

        _gFuncVal.array[fi*gridPoints.n + gi + 3] =
            2.*f3 * cxcy * (7.*cz2 - distsq) * expoFactor;

        _gFuncVal.array[fi*gridPoints.n + gi + 4] =
            f3 * (cx2 - cy2) * (7.*cz2 - distsq) * expoFactor; 
        
        _gFuncVal.array[fi*gridPoints.n + gi + 5] =
            f4 * cycz * (3.*cx2 - cy2) * expoFactor;

        _gFuncVal.array[fi*gridPoints.n + gi + 6] =
            f4 * cxcz * (cx2 - cy2) * expoFactor;
        
        _gFuncVal.array[fi*gridPoints.n + gi + 7] =
            4.*f * cxcy * (cx2 - cy2) * expoFactor;

        _gFuncVal.array[fi*gridPoints.n + gi + 8] =
            f * (cx2*cx2 - 6*cx2*cy2 + cy2*cy2)*expoFactor; 
    }
}

}}}
