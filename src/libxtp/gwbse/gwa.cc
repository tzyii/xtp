/*
 *            Copyright 2009-2016 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>

#include <votca/xtp/gwbse.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>

#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/constants.h>

//#include "mathimf.h"

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        // +++++++++++++++++++++++++++++ //
        // MBPT MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //

        
        

        void GWBSE::FullQPHamiltonian(){
            
            // constructing full QP Hamiltonian, storage in vxc
            _vxc = -_vxc + _sigma_x + _sigma_c;
            // diagonal elements are given by _qp_energies
            for (unsigned _m = 0; _m < _vxc.size1(); _m++ ){
              _vxc( _m,_m ) = _qp_energies( _m + _qpmin );
            }

            
            // sigma matrices can be freed
            _sigma_x.resize(0);
            _sigma_c.resize(0);
            
            
            if ( _do_qp_diag ){
            /* diagonalize, since _vxc will be needed in BSE, and GSL
             * destroys the input array, we need to make a local copy first
             */
            
            // get eigenvalues and eigenvectors of this matrix
            ub::matrix<double> _temp = _vxc;
            _qp_diag_energies.resize(_temp.size1());
            _qp_diag_coefficients.resize(_temp.size1(), _temp.size1());
            linalg_eigenvalues(_temp, _qp_diag_energies, _qp_diag_coefficients);

            // TODO storage -> orbitals

            }
            
        }

        void GWBSE::sigma_c_setup(const TCMatrix& _Mmn) {

            // iterative refinement of qp energies
            int _max_iter = 15;
            unsigned _levelsum = _Mmn[0].size2(); // total number of bands
            unsigned _gwsize = _Mmn[0].size1(); // size of the GW basis
            const double pi = boost::math::constants::pi<double>();

            ub::vector<double>& dftenergies=_orbitals->MOEnergies();
            // initial _qp_energies are dft energies
            _qp_energies = _orbitals->MOEnergies(); // RANGES!
            ub::vector<double>_qp_old=_qp_energies;
            double _DFTgap =dftenergies(_homo + 1) - dftenergies(_homo);
            bool energies_converged=false;
            _sigma_c.resize(_qptotal);
	    // only diagonal elements except for in final iteration
            for (int _i_iter = 0; _i_iter < _max_iter - 1; _i_iter++) {
                // loop over all GW levels
                #pragma omp parallel for
                for (unsigned _gw_level = 0; _gw_level < _qptotal; _gw_level++) {
                    
                    _sigma_c(_gw_level, _gw_level) = 0;
                    double qpmin = _qp_energies(_gw_level + _qpmin);
                    const ub::matrix<real_gwbse>& Mmn = _Mmn[ _gw_level + _qpmin ];

                    // loop over all functions in GW basis
                    for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
                        
                        double ppm_freq = _ppm_freq(_i_gw);
                        double fac = _ppm_weight(_i_gw) * ppm_freq;
                        // loop over all bands
                        for (unsigned _i = 0; _i < _levelsum; _i++) {

                            double occ = 1.0;
                            if (_i > _homo) occ = -1.0; // sign for empty levels

                            // energy denominator
                            double _denom = qpmin - _qp_energies(_i) + occ*ppm_freq;

                            double _stab = 1.0;
                            if (std::abs(_denom) < 0.25) {
                                _stab = 0.5 * (1.0 - std::cos(4.0 * pi * std::abs(_denom)));
                            }

                            double _factor =0.5* fac * _stab / _denom; //Hartree

                            // sigma_c diagonal elements
                            _sigma_c(_gw_level, _gw_level) += _factor * Mmn(_i_gw, _i) * Mmn(_i_gw, _i);

                        }// bands

                    }// GW functions

                    // update _qp_energies
                    _qp_energies(_gw_level + _qpmin) = dftenergies(_gw_level + _qpmin) + _sigma_x(_gw_level, _gw_level) + _sigma_c(_gw_level, _gw_level) - _vxc(_gw_level, _gw_level);


                }// all bands
                //cout << " end of qp refinement step (diagonal) " << _i_iter << "\n" << endl;
                _qp_old = _qp_old - _qp_energies;
                energies_converged = true;
                for (unsigned l = 0; l < _qp_old.size(); l++) {
                    if (std::abs(_qp_old(l)) > _qp_limit) {
                        energies_converged = false;
                        break;
                    }
                }
                if (energies_converged) {
                    LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Converged after " << _i_iter << " qp_energy iterations." << flush;
                    break;
                } else {
                    _qp_old = _qp_energies;
                }

            } // iterations

             double _QPgap = _qp_energies( _homo +1 ) - _qp_energies( _homo  );
             double _shift_new = _QPgap - _DFTgap;
             
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << (format(" New shift [Ryd] : %1$+1.6f ") % _shift_new ).str() << flush;
            //cout << " shift new " << _shift_new << endl;
            if (std::abs((_shift_new - _shift)) > _shift_limit) {
                _shift = _shift_new;
            } else {
                _shift_converged = true;
            }


	    if ( ! _iterate_shift  ) _shift_converged = true;


            // only if _shift is converged
            if (_shift_converged) {
                // in final step, also calc offdiagonal elements
                // initialize sigma_c to zero at the beginning
                
                
             //this is not the fastest algorithm but faster ones throw igwbse off, so this is good enough.    
                
            #pragma omp parallel for
            for (unsigned _gw_level = 0; _gw_level < _qptotal; _gw_level++) {
                double qpmin_e=_qp_energies(_gw_level + _qpmin);

                const ub::matrix<real_gwbse>& Mmn = _Mmn[ _gw_level + _qpmin ];
                for (unsigned _m = 0; _m < _gw_level; _m++) {
                    _sigma_c(_gw_level, _m) = 0;
                    const ub::matrix<real_gwbse>& Mmn2 = _Mmn[_m + _qpmin];

                    // loop over all functions in GW basis
                    for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
                        double ppm_freq = _ppm_freq(_i_gw);
                        double fac = _ppm_weight(_i_gw) * ppm_freq;
                        // loop over all screening levels
                        for (unsigned _i = 0; _i < _levelsum; _i++) {

                            double occ = 1.0;
                            if (_i > _homo) occ = -1.0; // sign for empty levels

                            // energy denominator
                            double _denom = qpmin_e - _qp_energies(_i) + occ * ppm_freq;

                            double _stab = 1.0;
                            if (std::abs(_denom) < 0.25) {
                                _stab = 0.5 * (1.0 - std::cos(4.0 * pi * std::abs(_denom)));
                            }

                            double _factor = 0.5*fac * Mmn(_i_gw, _i) * _stab / _denom; //Hartree

                            _sigma_c(_gw_level, _m) += _factor * Mmn2(_i_gw, _i);


                        }// screening levels 
                    }// GW functions 
                }// GW row 
                _qp_energies(_gw_level + _qpmin) = dftenergies(_gw_level + _qpmin) + _sigma_x(_gw_level, _gw_level) + _sigma_c(_gw_level, _gw_level) - _vxc(_gw_level, _gw_level);
            } // GW col 
        } 
            
        return;
        } // sigma_c_setup


        void GWBSE::sigma_x_setup(const TCMatrix& _Mmn){
        
            // initialize sigma_x
            _sigma_x.resize(_qptotal);
            int _size  = _Mmn[0].size1();

            // band 1 loop over all GW levels
            #pragma omp parallel for
            for ( unsigned _m1 = 0 ; _m1 < _qptotal ; _m1++ ){
                
                const ub::matrix<real_gwbse>& M1mn =  _Mmn[ _m1 + _qpmin ];
                
                // band 2 loop over all GW levels
                //for ( int _m2 = _qpmin ; _m2 <= _qpmax ; _m2++ ){
                for ( unsigned _m2 = 0 ; _m2 <= _m1 ; _m2++ ){
                    _sigma_x( _m1, _m2 )=0;
                    const ub::matrix<real_gwbse>& M2mn =  _Mmn[ _m2 + _qpmin ];
                    
                    // loop over all basis functions
                    for ( int _i_gw = 0 ; _i_gw < _size ; _i_gw++ ){
                        // loop over all occupied bands used in screening
                        for ( unsigned _i_occ = 0 ; _i_occ <= _homo ; _i_occ++ ){
                            _sigma_x( _m1, _m2 ) -= M1mn( _i_gw , _i_occ ) * M2mn( _i_gw , _i_occ );
                        } // occupied bands
                    } // gwbasis functions
                } // level 2
            } // level 1

	    // factor for hybrid DFT
	    _sigma_x = ( 1.0 - _ScaHFX ) * _sigma_x;
            return;
        }


        void GWBSE::sigma_prepare_threecenters(TCMatrix& _Mmn){
            #pragma omp parallel for
            for ( int _m_level = 0 ; _m_level < _Mmn.get_mtot(); _m_level++ ){
                // get Mmn for this _m_level
                // and multiply with _ppm_phi = eigenvectors of epsilon
                // POTENTIAL BUG
	      // casting _Mmn to double for efficint prod() overload
	      ub::matrix<double> _Mmn_double = _Mmn[_m_level];
              ub::matrix<real_gwbse> _temp = ub::prod(  _ppm_phi , _Mmn_double );
                _Mmn[ _m_level ] = _temp;
            }
            return;
        }        
        


    }
    
 
};
