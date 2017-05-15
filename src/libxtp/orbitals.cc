/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#include "votca/xtp/orbitals.h"
#include "votca/tools/globals.h"
#include "votca/xtp/elements.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <iterator>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>



namespace votca { namespace xtp {

template <class VEC>
class vector_less
{
        private:
        typedef typename VEC::size_type  size_type;
        typedef typename VEC::value_type value_type;
        vector_less();

        const VEC& data;
        public:
        vector_less(const VEC& vec) : data(vec) { }
        bool operator() (const size_type& left, const size_type& right) const
        {
                return std::less<value_type> () ( data(left), data(right) );
        }
};   
    
Orbitals::Orbitals() { 
    
    _basis_set_size = 0;
    _occupied_levels = 0;
    _unoccupied_levels = 0;
    _number_of_electrons = 0;
    _self_energy = 0.0;
    _qm_energy = 0.0;
    _couplingsA = 0;
    _couplingsB = 0;
    _ECP="";

    //_has_atoms = false;

    
    // GW-BSE
    _qpmin = 0;
    _qpmax = 0;
    _qptotal = 0;
    
    _rpamin = 0;
    _rpamax = 0;
    
    _ScaHFX = 0.0;
    
    _bse_cmin = 0;
    _bse_cmax = 0;
    _bse_vmin = 0;
    _bse_vmax = 0;
    _bse_vtotal = 0;
    _bse_ctotal = 0;
    _bse_size = 0;
    _bse_nmax = 0;
  
};   

Orbitals::~Orbitals() { 
    _mo_energies.clear();
    _mo_coefficients.clear();
    _overlap.clear();
    
    std::vector< ctp::QMAtom* >::iterator it;
    for ( it = _atoms.begin(); it != _atoms.end(); ++it ) delete *it;
    
};   

/*
const int    &Orbitals::getBasisSetSize() const { 
    if ( _has_basis_set_size ) {
        return _basis_set_size; 
    } else {
        throw std::runtime_error(" Basis set size is unknown. Parse a log file first. " );
    }
}
*/

//void          Orbitals::setBasisSetSize( const int &basis_set_size ){
//    _has_basis_set_size = true; 
//    _basis_set_size = basis_set_size; 
// }

/*
int     Orbitals::getNumberOfLevels() {
    if ( _has_occupied_levels && _has_unoccupied_levels ) {
        return  _occupied_levels + _unoccupied_levels; 
    } else {
        throw std::runtime_error(" Number of levels is unknown. Parse a log file first. " );
    }
}
*/

void   Orbitals::setNumberOfLevels( const int &occupied_levels,const int &unoccupied_levels  ){
    // _has_occupied_levels = true; 
    // _has_unoccupied_levels = true; 
    _occupied_levels = occupied_levels; 
    _unoccupied_levels = unoccupied_levels; 
}
/*
const int     &Orbitals::getNumberOfElectrons() const {
    if ( _has_number_of_electrons ) {
        return  _number_of_electrons; 
    } else {
        throw std::runtime_error(" Number of electrons is unknown. Parse a log file first. " );
    }
}
*/
//void          Orbitals::setNumberOfElectrons( const int &electrons ){
//    _has_number_of_electrons = true; 
//    _number_of_electrons = electrons; 
//}

/**
 * 
 * @param _energy_difference [ev] Two levels are degenerate if their energy is smaller than this value
 * @return A map with key as a level and a vector which is a list of close lying orbitals
 */
bool Orbitals::CheckDegeneracy( double _energy_difference ) {
    
    ub::vector<double>::iterator it1 = _mo_energies.begin();
    bool _degenerate = false;
    
    if ( tools::globals::verbose ) cout << endl <<"... ... Checking level degeneracy " << endl;
    
    _level_degeneracy.clear();
            
    while ( it1 !=_mo_energies.end() ) {
        
        // in all containers counters start with 0; real life - with 1
        int _level1 = std::distance(_mo_energies.begin(), it1) + 1;
        
        // add the level itself - it is easier to loo over all levels later
        _level_degeneracy[_level1].push_back( _level1 );        
        
        ub::vector<double>::iterator it2 = it1;
        it2++;
        
        while (  it2 !=_mo_energies.end() ) {
            //cout << _level1 << ":" << *it1 << ":" << *it2 << endl;
            double energy1 = *it1;
            double energy2 = *it2;
            
            // in all containers counters start with 0; real life - with 1
            int _level2 = std::distance(_mo_energies.begin(), it2) + 1;
            
            if ( std::abs(energy1 - energy2)*tools::conv::hrt2ev < _energy_difference ) {
                _level_degeneracy[_level1].push_back( _level2 );
                _level_degeneracy[_level2].push_back( _level1 );
                _degenerate = true;
            }
            it2++;
        }
        it1++;
    }

    if ( tools::globals::verbose ){ 

        if ( _degenerate ) {
            cout << "... ... Some levels are degenerate" << endl; 
            for (std::map<int, std::vector<int> >::iterator it = _level_degeneracy.begin();
                    it != _level_degeneracy.end();
                    ++it) {
                // output only degenerate levels
                if ( (it->second).size() > 1 ) {
                    std::cout << "... ... level  "<< it->first << " : ";
                    for (vector<int>::iterator lev = (it->second).begin(); lev != (it->second).end(); lev++)
                            cout << *lev << " ";
                    cout << endl;
                }
            }
        } else {
            cout << "... ... No degeneracy found" << endl;  
        }

        cout << "... ... Done checking level degeneracy" << endl;   
    
    }
    
    // _has_level_degeneracy = true;
    return _degenerate;
    
}    

std::vector<int>* Orbitals::getDegeneracy( int level, double _energy_difference ) {
    if ( ! hasDegeneracy() ) {
        
        CheckDegeneracy( _energy_difference );       
        /* 

        int _ld = _level_degeneracy.at(level).size();       
        if ( _ld > 1 ) {
                cout << "....level " << level << " degeneracy is: " <<  _ld << endl;
        } else {
                cout << "....level " << level << " is not degenerate" << endl;
        }
        */
        
    }        
    //cout << "Getting degeneracy for level " << level << endl;
    return &_level_degeneracy.at(level);
}

void Orbitals::SortEnergies(  std::vector<int>* index ) {
    if ( tools::globals::verbose )  cout << "... ... Sorting energies" << endl ;
    //cout << _mo_energies << endl;
    //cout << "MO Energies size" << _mo_energies.size() << endl ;
    //exit(0);
    index->resize( _mo_energies.size() );
    int i = 0;
            for ( vector< int > ::iterator soi = index->begin(); soi != index->end(); ++soi ) {
                index->at(i) = i;
                i++;
                //cout << *soi << " ";
            } 
    
    //cout << _mo_energies;
    
    std::stable_sort(index->begin(), index->end(), vector_less< ub::vector<double> >( _mo_energies ));

}

/// Writes a PDB file
void Orbitals::WritePDB( FILE *out, string tag ) {
    // out.setf(ios::fixed);
    string tag_extended="HEADER ! GENERATED BY VOTCA::XTP::"+tag+"\n";
    fprintf(out,"%s",tag_extended.c_str() );
    vector < ctp::QMAtom* > :: iterator atom;
    int id = 0;
    
    //cout << "Natoms " << _atoms.size() << endl;
    
    for (atom = _atoms.begin(); atom < _atoms.end(); ++atom){
         id++;      
         string resname = ( (*atom)->from_environment ) ? "MM" : "QM";
         int resnr = 1;
         
         //cout << id << " " << (*atom)->type << " " << endl;
         
         fprintf(out, "ATOM  %5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s  %8.3f\n",
                 id,                    // Atom serial number           %5d 
                 (*atom)->type.c_str(), // Atom name                    %4s
                 " ",                   // alternate location indicator.%1s
                 resname.c_str(),       // Residue name.                %3s
                 "A",                   // Chain identifier             %1s
                 resnr,                 // Residue sequence number      %4d
                 " ",                   // Insertion of residues.       %1s
                 (*atom)->x,            // X in Angstroms               %8.3f
                 (*atom)->y,            // Y in Angstroms               %8.3f
                 (*atom)->z,            // Z in Angstroms               %8.3f
                 1.0,                   // Occupancy                    %6.2f
                 0.0,                   // Temperature factor           %6.2f
                 " ",                   // Segment identifier           %4s
                 (*atom)->type.c_str(), // Element symbol               %2s
                 (*atom)->charge        // Charge on the atom.          %2s
                 );
    }  
}

// reduces the number of virtual orbitals to factor*number_of_occupied_orbitals
void Orbitals::Trim( int factor ) {
    
    if ( hasMOCoefficients() ) {
        _mo_coefficients.resize ( factor * _occupied_levels, _basis_set_size, true);
        _unoccupied_levels = ( factor -1 ) * _occupied_levels;        
    }

    if ( hasMOEnergies() ) {
        //cout << "\nBefore: " << _mo_energies.size();
        _mo_energies.resize(  factor * _occupied_levels, true );
        _unoccupied_levels = ( factor - 1) * _occupied_levels;
        //cout << " and After: " << _mo_energies.size() << endl;   
    }
}


// reduces the number of orbitals to [HOMO-degH:LUMO+degL]
void Orbitals::Trim( int degH, int degL ) {
    
    
    
    if ( hasMOCoefficients() ) {
        
        //cout << " COEF Before " << _mo_coefficients.size1() << endl;   
        
        _mo_coefficients = ub::project( _mo_coefficients , ub::range(_occupied_levels -degH , _occupied_levels + degL), ub::range ( 0, _basis_set_size ) ) ;
        
        //cout << " and After: " << _mo_coefficients.size1() << endl;   
        
    }

    if ( hasMOEnergies() ) {
        
        //cout << " EN Before " << _mo_energies.size() << endl;   
        
        ub::vector<double> _temp(degH+degL); 
        
        for ( int i = 0; i< degH+degL ; i++ ) {
            
            _temp(i) = _mo_energies(_occupied_levels-degH + i);
            
        }
        
        
        _mo_energies = _temp;
        //cout << " and After: " << _mo_energies.size() << endl;   
    }

    _occupied_levels = degH;
    _unoccupied_levels = degL;        



}











bool Orbitals::Load(string file_name) {
    try {
        std::ifstream ifs( file_name.c_str() );
        boost::archive::binary_iarchive ia( ifs );
        ia >> *this;
        ifs.close();
    } catch(std::exception &err) {
        std::cerr << "Could not load orbitals from " << file_name << flush; 
        std::cerr << "An error occurred:\n" << err.what() << endl;
        return false;
    } 
    return true;
}






 // Determine ground state density matrix
 ub::matrix<double> Orbitals::DensityMatrixGroundState(const ub::matrix<double>& _MOs ) {   
     // first fill Density matrix, if required
    //  if ( _dmatGS.size1() != _basis_set_size ) {
        ub::matrix<double> dmatGS = ub::zero_matrix<double>(_basis_set_size, _basis_set_size);
        #pragma omp parallel for
        for ( int _i=0; _i < _basis_set_size; _i++ ){
            for ( int _j=0; _j < _basis_set_size; _j++ ){
                for ( int _level=0; _level < _occupied_levels ; _level++ ){
                 
                    dmatGS(_i,_j) += 2.0 * _MOs( _level , _i ) * _MOs( _level , _j );
                 
                }
            }
         }
     //}    
     // return     
     return dmatGS;  
 }
 
 
 ub::matrix<double>  Orbitals::TransitionDensityMatrix(const ub::matrix<double>& _MOs ,const ub::matrix<real_gwbse>& _BSECoefs, int state){
    ub::matrix<double> dmatTS=ub::zero_matrix<double>(_basis_set_size);
    // The Transition dipole is sqrt2 bigger because of the spin, the excited state is a linear combination of 2 slater determinants, where either alpha or beta spin electron is excited
    double sqrt2=sqrt(2.0);
    /*Trying to implement D_{alpha,beta}= sqrt2*sum_{i}^{occ}sum_{j}^{virt}{BSEcoef(i,j)*MOcoef(alpha,i)*MOcoef(beta,j)} */
    // c stands for conduction band and thus virtual orbitals
    // v stand for valence band and thus occupied orbitals


     
    if ( _bse_size == 0 ) {
        _bse_vtotal = _bse_vmax - _bse_vmin +1 ;
        _bse_ctotal = _bse_cmax - _bse_cmin +1 ;
        _bse_size   = _bse_vtotal * _bse_ctotal;
          // indexing info BSE vector index to occupied/virtual orbital
        for ( unsigned _v = 0; _v < _bse_vtotal; _v++ ){
            for ( unsigned _c = 0; _c < _bse_ctotal ; _c++){
                  _index2v.push_back( _bse_vmin + _v );
                  _index2c.push_back( _bse_cmin + _c );
            }
        }
    }
    
    #pragma omp parallel for
    for(unsigned a=0;a<dmatTS.size1();a++){
        for(unsigned b=0;b<dmatTS.size2();b++){
            for(unsigned i=0;i<_bse_size;i++){
                int occ=_index2v[i];
                int virt=_index2c[i];
                dmatTS(a,b)+=sqrt2*_BSECoefs(i,state)*_MOs( occ , a ) * _MOs( virt , b ); //check factor 2??
            }
        }     
    }

    return dmatTS;
}
 
 
 
 // Excited state density matrix
std::vector<ub::matrix<double> > Orbitals::DensityMatrixExcitedState(const ub::matrix<double>& _MOs,const ub::matrix<real_gwbse>& _BSECoefs, int state ){
     
     
     /****** 
      * 
      *    Density matrix for GW-BSE based excitations
      * 
      *    - electron contribution
      *      D_ab = \sum{vc} \sum{c'} A_{vc}A_{vc'} mo_a(c)mo_b(c')
      * 
      *    - hole contribution 
      *      D_ab = \sum{vc} \sum{v'} A_{vc}A_{v'c} mo_a(v)mo_b(v')
      * 
      * 
      *   more efficient:
      * 
      *   - electron contribution
      *      D_ab = \sum{c} \sum{c'} mo_a(c)mo_b(c') [ \sum{v} A_{vc}A_{vc'} ]
      *           = \sum{c} \sum{c'} mo_a(c)mo_b(c') A_{cc'} 
      *    
      *   - hole contribution
      *      D_ab = \sum{v} \sum{v'} mo_a(v)mo_b(v') [ \sum{c} A_{vc}A_{v'c} ]
      *           = \sum{v} \sum{v'} mo_a(v)mo_b(v') A_{vv'} 
      *  
      */
     std::vector<ub::matrix<double> > dmatEX;        
    dmatEX.resize(2);
    dmatEX[0] = ub::zero_matrix<double>(_basis_set_size, _basis_set_size);
    dmatEX[1] = ub::zero_matrix<double>(_basis_set_size, _basis_set_size);

    int _vmin = this->_bse_vmin;
    int _vmax = this->_bse_vmax;
    int _cmin = this->_bse_cmin;
    int _cmax = this->_bse_cmax;
     
    if ( _bse_size == 0 ) {
        _bse_vtotal = _bse_vmax - _bse_vmin +1 ;
        _bse_ctotal = _bse_cmax - _bse_cmin +1 ;
        _bse_size   = _bse_vtotal * _bse_ctotal;
          // indexing info BSE vector index to occupied/virtual orbital
          for ( unsigned _v = 0; _v < _bse_vtotal; _v++ ){
              for ( unsigned _c = 0; _c < _bse_ctotal ; _c++){
                  _index2v.push_back( _bse_vmin + _v );
                  _index2c.push_back( _bse_cmin + _c );
               }
           }
     }
     
     // electron assist matrix A_{cc'}
     ub::matrix<real_gwbse> _Acc = ub::zero_matrix<real_gwbse>( this->_bse_ctotal , this->_bse_ctotal );
     ub::matrix<real_gwbse> _Avv = ub::zero_matrix<real_gwbse>( this->_bse_vtotal , this->_bse_vtotal );
  
     for ( unsigned _idx1 = 0 ; _idx1 < _bse_size ; _idx1++) {
         
         int _v = this->_index2v[_idx1];
         int _c = this->_index2c[_idx1];

         // electron assist matrix A_{cc'}
         #pragma omp parallel for
         for ( int _c2=_cmin; _c2 <= _cmax; _c2++ ){
             int _idx2 = (_cmax-_cmin+1)*(_v-_vmin)+(_c2-_cmin);
             
             _Acc(_c - _cmin ,_c2 - _cmin) +=   _BSECoefs(_idx1,state) * _BSECoefs(_idx2,state) ;
         }
         
         // hole assist matrix A_{vv'}
         #pragma omp parallel for
         for ( int _v2=_vmin; _v2 <= _vmax; _v2++ ){
                int _idx2 = (_cmax-_cmin+1)*(_v2-_vmin)+(_c-_cmin);
                
                _Avv(_v - _vmin ,_v2 - _vmin ) +=   _BSECoefs(_idx1,state) *_BSECoefs(_idx2,state) ;
    
            }
         
     }
     
     // hole part as matrix products
     // get slice of MOs of occs only
     ub::matrix<double> _occlevels = ub::project(_MOs, ub::range(_vmin, _vmax + 1), ub::range(0, _basis_set_size));
     ub::matrix<double> _temp = ub::prod( _Avv, _occlevels );
     dmatEX[0] = ub::prod(ub::trans(_occlevels), _temp);
     
     
     // electron part as matrix products
     // get slice of MOs of virts only
     ub::matrix<double> _virtlevels = ub::project(_MOs, ub::range(_cmin, _cmax + 1), ub::range(0, _basis_set_size));
     _temp = ub::prod( _Acc, _virtlevels );
     dmatEX[1] = ub::prod(ub::trans(_virtlevels), _temp);
     
     return dmatEX;
 }
 

 ub::vector<double> Orbitals::MullikenPopulation( const ub::matrix<double>& _densitymatrix, const ub::matrix<double>& _overlapmatrix, int _frag) {
     
     ub::vector<double> fragmentCharges=ub::vector<double>(2,0.0);
     
     
     ub::matrix<double> _prodmat = ub::prod( _densitymatrix, _overlapmatrix );
         
     for ( int _i = 0 ; _i < _frag; _i++){
        fragmentCharges(0) += _prodmat(_i,_i);
     }
     for ( unsigned _i = _frag ; _i < _overlapmatrix.size1(); _i++){
       fragmentCharges(1) += _prodmat(_i,_i);
     }
           
   return fragmentCharges; 
 }


 ub::vector<double> Orbitals::FragmentNuclearCharges(int _frag){
     Elements _elements;
     
     // go through atoms and count
    vector < ctp::QMAtom* > :: iterator atom;
    int id = 0;
    
    if ( _frag < 0 ) {
        throw runtime_error("Orbitals::FragmentNuclearCharges Fragment index is smaller than zero");
    }
    
    ub::vector<double> fragmentNuclearCharges=ub::vector<double>(2,0.0);
    
    for (atom = _atoms.begin(); atom < _atoms.end(); ++atom){
         id++;      
         // get element type and determine its nuclear charge
         double crg = _elements.getNucCrgECP((*atom)->type);
         // add to either fragment
         if ( id <= _frag ) {
             fragmentNuclearCharges(0) += crg;
         } else {
             fragmentNuclearCharges(1) += crg;
         }     
    }

    return fragmentNuclearCharges;
 }
 
  /** 
     * \brief Guess for a dimer based on monomer orbitals
     * 
     * Given two monomer orbitals (A and B) constructs a guess for dimer
     * orbitals: | A 0 | and energies: [EA, EB]
     *           | 0 B |
     */
void Orbitals::PrepareGuess( Orbitals* _orbitalsA, Orbitals* _orbitalsB, Orbitals* _orbitalsAB) {
       
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
       
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    int _electronsA = _orbitalsA->getNumberOfElectrons();
    int _electronsB = _orbitalsB->getNumberOfElectrons();
    
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
    
    ub::matrix<double>& _mo_coefficients = _orbitalsAB->MOCoefficients();    
    
    
    // AxB = | A 0 |  //   A = [EA, EB]  //
    //       | 0 B |  //                 //
    _mo_coefficients.resize( _levelsA + _levelsB, _basisA + _basisB  );
    _orbitalsAB->setBasisSetSize( _basisA + _basisB );
    _orbitalsAB->setNumberOfLevels( _electronsA - _electronsB , 
                                    _levelsA + _levelsB - _electronsA - _electronsB );
    _orbitalsAB->setNumberOfElectrons( _electronsA + _electronsB );
    
    ub::project( _mo_coefficients, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _mo_coefficients, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _mo_coefficients, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = _orbitalsA->MOCoefficients();
    ub::project( _mo_coefficients, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = _orbitalsB->MOCoefficients();   

    ub::vector<double>& _energies = _orbitalsAB->MOEnergies();
    _energies.resize( _levelsA + _levelsB );
     
    ub::project( _energies, ub::range (0, _levelsA ) ) = _orbitalsA->MOEnergies();
    ub::project( _energies, ub::range (_levelsA, _levelsA + _levelsB ) ) = _orbitalsB->MOEnergies();
    
    return;
}   
 
 

 
}}
