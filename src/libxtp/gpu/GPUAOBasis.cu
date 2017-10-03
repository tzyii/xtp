#include <votca/xtp/gpu/GPUAOBasis.h>

namespace votca { namespace xtp { namespace gpu {

GPUAOBasis::GPUAOBasis(){}
GPUAOBasis::~GPUAOBasis(){}

GPUAOBasis::GPUAOBasis(const AOBasis& aob){
    for (AOBasis::AOShellIterator row = aob.firstShell(); row != aob.lastShell(); row++){
        const AOShell* shell = aob.getShell(row); 
        const std::string shell_type = shell->getType();
        vec shellPos = shell->getPos();
        
        for (AOShell::GaussianIterator itr = shell->firstGaussian(); itr != shell->lastGaussian(); ++itr){
            const std::vector<double>& contractions = (*itr)->getContraction();
            alphas.push_back((*itr)->getDecay());
            powFactors.push_back((*itr)->getPowfactor());            
            
            for (const char& c: shell_type){
                switch (::toupper(c)){
                case 'S':
                    sConts.push_back(contractions[0]);
                    sPos.x.push_back(shellPos.getX());
                    sPos.y.push_back(shellPos.getY());
                    sPos.z.push_back(shellPos.getZ());
                    
                case 'P':
                    pConts.push_back(contractions[1]);
                    pPos.x.push_back(shellPos.getX());
                    pPos.y.push_back(shellPos.getY());
                    pPos.z.push_back(shellPos.getZ());
                    
                case 'D':
                    dConts.push_back(contractions[2]);
                    dPos.x.push_back(shellPos.getX());
                    dPos.y.push_back(shellPos.getY());
                    dPos.z.push_back(shellPos.getZ());
                    
                case 'F':
                    fConts.push_back(contractions[3]);
                    fPos.x.push_back(shellPos.getX());
                    fPos.y.push_back(shellPos.getY());
                    fPos.z.push_back(shellPos.getZ());
                    
                case 'G':
                    gConts.push_back(contractions[4]);
                    gPos.x.push_back(shellPos.getX());
                    gPos.y.push_back(shellPos.getY());
                    gPos.z.push_back(shellPos.getZ());
                    
                case 'H':
                    std::cerr << "H functions not implemented at the moment!" << std::endl;
                    exit(EXIT_FAILURE);
                default:
                    std::cerr << "Shell type " << c << "unknown" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    
    Pad(sConts, 16);
    Pad(pConts, 16);
    Pad(dConts, 16);
    Pad(fConts, 16);
    Pad(gConts, 16);
    Pad(alphas, 16);
    Pad(powFactors, 16);

    d_sConts = gpu_vector(sConts);
    d_pConts = gpu_vector(pConts);
    d_dConts = gpu_vector(dConts);
    d_fConts = gpu_vector(fConts);
    d_gConts = gpu_vector(gConts);
    d_alphas = gpu_vector(alphas);
        
    d_powFactors = gpu_vector(powFactors);
        
    d_expoFactors = gpu_vector(alphas.size(), 0);

    // build the pointer struct that will be passed to the gpu kernels. 
    rawGpuArrs.sConts.array = thrust::raw_pointer_cast(&d_sConts[0]);
    rawGpuArrs.sConts.arraySize = d_sConts.size();

    rawGpuArrs.pConts.array = thrust::raw_pointer_cast(&d_pConts[0]);
    rawGpuArrs.pConts.arraySize = d_pConts.size();
        
    rawGpuArrs.dConts.array = thrust::raw_pointer_cast(&d_dConts[0]);
    rawGpuArrs.dConts.arraySize = d_dConts.size();

    rawGpuArrs.fConts.array = thrust::raw_pointer_cast(&d_fConts[0]);
    rawGpuArrs.fConts.arraySize = d_fConts.size();

    rawGpuArrs.gConts.array = thrust::raw_pointer_cast(&d_gConts[0]);
    rawGpuArrs.gConts.arraySize = d_gConts.size();

    rawGpuArrs.alphas.array = thrust::raw_pointer_cast(&d_alphas[0]);
    rawGpuArrs.alphas.arraySize = d_alphas.size();

    rawGpuArrs.powFactors.array = thrust::raw_pointer_cast(&d_powFactors[0]);
    rawGpuArrs.powFactors.arraySize = d_powFactors.size();

    rawGpuArrs.expoFactors.array = thrust::raw_pointer_cast(&d_expoFactors[0]);
    rawGpuArrs.expoFactors.arraySize = d_expoFactors.size();

    Pad(sPos, 16);
    Pad(pPos, 16);
    Pad(dPos, 16);
    Pad(fPos, 16);
    Pad(gPos, 16);

    d_sPos = device_R3N_vectors(sPos); 
    d_pPos = device_R3N_vectors(pPos); 
    d_dPos = device_R3N_vectors(dPos); 
    d_fPos = device_R3N_vectors(fPos); 
    d_gPos = device_R3N_vectors(gPos); 

    rawGpuArrs.sFuncPos.x = thrust::raw_pointer_cast(&d_sPos.x[0]);
    rawGpuArrs.sFuncPos.y = thrust::raw_pointer_cast(&d_sPos.y[0]);
    rawGpuArrs.sFuncPos.z = thrust::raw_pointer_cast(&d_sPos.z[0]);
    rawGpuArrs.sFuncPos.n = d_sPos.x.size();

    rawGpuArrs.pFuncPos.x = thrust::raw_pointer_cast(&d_pPos.x[0]);
    rawGpuArrs.pFuncPos.y = thrust::raw_pointer_cast(&d_pPos.y[0]);
    rawGpuArrs.pFuncPos.z = thrust::raw_pointer_cast(&d_pPos.z[0]);
    rawGpuArrs.pFuncPos.n = d_pPos.x.size();

    rawGpuArrs.dFuncPos.x = thrust::raw_pointer_cast(&d_dPos.x[0]);
    rawGpuArrs.dFuncPos.y = thrust::raw_pointer_cast(&d_dPos.y[0]);
    rawGpuArrs.dFuncPos.z = thrust::raw_pointer_cast(&d_dPos.z[0]);
    rawGpuArrs.dFuncPos.n = d_dPos.x.size();

    rawGpuArrs.fFuncPos.x = thrust::raw_pointer_cast(&d_fPos.x[0]);
    rawGpuArrs.fFuncPos.y = thrust::raw_pointer_cast(&d_fPos.y[0]);
    rawGpuArrs.fFuncPos.z = thrust::raw_pointer_cast(&d_fPos.z[0]);
    rawGpuArrs.fFuncPos.n = d_fPos.x.size();

    rawGpuArrs.gFuncPos.x = thrust::raw_pointer_cast(&d_gPos.x[0]);
    rawGpuArrs.gFuncPos.y = thrust::raw_pointer_cast(&d_gPos.y[0]);
    rawGpuArrs.gFuncPos.z = thrust::raw_pointer_cast(&d_gPos.z[0]);
    rawGpuArrs.gFuncPos.n = d_gPos.x.size();

}
    

gpuAOArrs GPUAOBasis::GetRawGPUArrs(){
    return rawGpuArrs;
}

}}}
