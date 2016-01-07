/* 
 *            Copyright 2009-2015 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/ctp/votca_ctp_config.h>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/algorithm/string.hpp>
#include <votca/ctp/logger.h>
#include <votca/ctp/gdma_internal.h>



using namespace votca::tools;

namespace votca {
    namespace ctp {
        namespace ub = boost::numeric::ublas;

int GDMA_internal::Run_gdma( Orbitals *orbitals){
    const double itol=18; 
    const double tol=2.30258*itol;
    const double bigexp=4.0; //in dma this is bigexp_default
    std::vector< QMAtom* > atoms=orbitals->QMAtoms();
    sites=PolarSeg(-1,atoms);
    
    BasisSet bs;
    bs.LoadBasisSet(orbitals->getDFTbasis());
    AOBasis dftbasis;
    dftbasis.AOBasisFill(&bs, atoms);
    AOBasis basis;    
    basis.ReorderMOs(orbitals.MOCoefficients(), orbitals.getQMpackage(), "votca" );  
    
    // extend to all density matrices
    ub::matrix<double> dmat=orbitals->DensityMatrixGroundState();
    
    //looping over pairs of shells
    
    vector < QMAtom* > :: iterator atom1;
    vector < QMAtom* > :: iterator atom2;
    Element::ShellIterator eits1;
    Element::ShellIterator eits2;
    Shell::GaussianIterator sits1;
    Shell::GaussianIterator sits2;
    
    // startshell is necessary to identify shell with dmat entry
    int startshell1=0;
    //loop first atom
    for (atom1 = atoms.begin(); atom1 < atoms.end(); ++atom1) {
        vec r1=atom1->getPosition();

        Element* element1 = bs.getElement((*atom1)->type);
        for (eits1 = element1->firstShell(); eits1 != element1->lastShell(); eits1++) {
            Shell* shell1 = (*eits1);
            int s1_size=shell1->getnumofFunc();
            int l1=shell1->getLmax();
            int startshell2=0;
            // loop second atom
            for (atom2 = atoms.begin(); atom2 < atoms.end(); ++atom2) {
                vec r2=atom2->getPosition();
                vec r=r1-r2
                double dist=abs(r);
                            
                Element* element2 = bs.getElement((*atom2)->type);
                for (eits2 = element1->firstShell(); eits2 != element1->lastShell(); eits2++) {
                     Shell* shell2 = (*eits2);
                     int s2_size=shell2->getnumofFunc();
                     ub::matrix<double> dmat_shellpair=ub::zero_matrix<double>(s1_size,s2_size);
                     int l2=shell2->getLmax();
                     for(sits1=shell1->firstGaussian();sits1<shell1->lastGaussian();++sits1){
                         double alpha=sits1->decay;
                         for(sits2=shell2->firstGaussian();sits2<shell2->lastGaussian();++sits2){
                            double beta=sits2->decay;  
                            double ab=alpha+beta;
                            double dum=alpha*beta*dist*dist/ab;
                            if (dum>tol) continue;
                            if (ab>bigexp){
                                double p=alpha/ab;
                                vec posa=p*r;
                                vec posb=(1-p)*r;
                                vec posp=(alpha*r1+beta*r2)/ab;
                                double t=1/sqrt(ab);
                                int n=l1+l2+1; //max norm for gauss hermite quad
                                
                                        
                            }
                         }
                     }                         
                     
                     
                     
                     
                     
                     
                     startshell2+=s2_size;
                }
            }
            startshell1+=s1_size;
        }
    }
}

void GDMA_internal::MoveQ(vec pos){
    
}



    }
}
