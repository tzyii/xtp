#ifndef GPU_CTYPES_CUH
#define GPU_CTYPES_CUH

namespace votca { namespace xtp { namespace gpu {
typedef struct{
    double* array;
    size_t arraySize;
} gpuArr;

typedef struct{
    double* x;
    double* y;
    double* z;

    size_t n;
} R3Nptrs;

typedef struct{
    gpuArr sConts; 
    gpuArr pConts; 
    gpuArr dConts; 
    gpuArr fConts; 
    gpuArr gConts;

    gpuArr alphas;
    gpuArr powFactors;
    gpuArr expoFactors;

    R3Nptrs sFuncPos; 
    R3Nptrs pFuncPos; 
    R3Nptrs dFuncPos; 
    R3Nptrs fFuncPos; 
    R3Nptrs gFuncPos; 
} gpuAOArrs;

typedef struct{
    gpuArr sFuncVals; 
    gpuArr pFuncVals; 
    gpuArr dFuncVals; 
    gpuArr fFuncVals; 
    gpuArr gFuncVals;
} GPUFuncVals;

}}}

#endif // GPU_CTYPES_CUH
