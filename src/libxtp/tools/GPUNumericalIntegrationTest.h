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

#ifndef _VOTCA_XTP_GPUNUMERICAL_INTEGRATION_TEST_H
#define _VOTCA_XTP_GPUNUMERICAL_INTEGRATION_TEST_H

#include <stdio.h>
#include <votca/ctp/logger.h>
// Overload of uBLAS prod function with MKL/GSL implementations

namespace votca {
    namespace xtp {
        using namespace std;
        class GPUNumericalIntegrationTest : public ctp::QMTool {
        public:

            GPUNumericalIntegrationTest() {
            };

            ~GPUNumericalIntegrationTest() {
            };

            string Identify() {
                return "GPUNumericalIntegrationTest";
            }

            void Initialize(Property *options);
            bool Evaluate();

        private:
            string _orbfile;
            ctp::Logger _log;

        };

        void GPUNumericalIntegrationTest::Initialize(Property* options) {
            // update options with the VOTCASHARE defaults
            UpdateWithDefaults( options, "xtp" );


            // orbitals file or pure DFT output
            _orbfile = options->get(key + ".input").as<string> ();

            // get the path to the shared folders with xml files
            char *votca_share = getenv("VOTCASHARE");
            if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
            // string xmlFile = string(getenv("VOTCASHARE")) + string("/xtp/packages/") + _package + string("_idft_pair.xml");
            // load_property_from_xml( _package_options, xmlFile );

            // register all QM packages (Gaussian, TURBOMOLE, etc)
            // QMPackageFactory::RegisterAll();

        }

bool GPUNumericalIntegrationTest::Evaluate() {
    _log.setReportLevel( ctp::logDEBUG );
    _log.setMultithreading( true );

    _log.setPreface(ctp::logINFO,    "\n... ...");
    _log.setPreface(ctp::logERROR,   "\n... ...");
    _log.setPreface(ctp::logWARNING, "\n... ...");
    _log.setPreface(ctp::logDEBUG,   "\n... ...");


    CTP_LOG(ctp::logDEBUG, _log) << "Reading serialized QM data from " << _orbfile << flush;

    Orbitals _orbitals;

    CTP_LOG(ctp::logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
    _orbitals.Load(_orbfile);

    // load DFT basis set (element-wise information) from xml file
    BasisSet dftbs;
    dftbs.LoadBasisSet(_orbitals.getDFTbasis());
    CTP_LOG(ctp::logDEBUG, _log) << " Loaded DFT Basis Set " << _orbitals.getDFTbasis() << flush;

    // fill DFT AO basis by going through all atoms
    AOBasis dftbasis;
    dftbasis.AOBasisFill(&dftbs, _orbitals.QMAtoms());

    ub::matrix<double> DMATGS = _orbitals.DensityMatrixGroundState();


    return true;
        }
}}


#endif
