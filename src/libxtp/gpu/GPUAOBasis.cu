#include <votca/xtp/gpu/GPUAOBasis.h>

namespace votca { namespace xtp { namespace gpu {

GPUAOBasis::GPUAOBasis(){}
GPUAOBasis::~GPUAOBasis(){}

GPUAOBasis::GPUAOBasis(const AOBasis& aob){
    for (AOBasis::AOShellIterator row = aob.firstShell(); row != aob.lastShell(); row++){
        const AOShell* shell = aob.getShell(row); 
        const std::string shell_type = shell->getType();
        
        for (AOShell::GaussianIterator itr = shell->firstGaussian(); itr != shell->lastGaussian(); ++itr){
            const std::vector<double>& contractions = (*itr)->getContraction();
            alphas.push_back((*itr)->getDecay());
            powFactors.push_back((*itr)->getPowfactor());
            
            for (const char& c: shell_type){
                switch (::toupper(c)){
                case 'S':
                    sConts.push_back(contractions[0]); 
                case 'P':
                    pConts.push_back(contractions[1]);
                case 'D':
                    dConts.push_back(contractions[2]);
                case 'F':
                    fConts.push_back(contractions[3]);
                case 'G':
                    gConts.push_back(contractions[4]);
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
}
    

gpuAOArrs GPUAOBasis::GetRawGPUArrs(){
    return rawGpuArrs;
}

}}}
