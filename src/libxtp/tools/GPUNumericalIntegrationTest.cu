#include <iostream>
#include <stdexcept>
#include <stdlib.h>

#include "GPUNumericalIntegrationTest.h"

namespace votca {
    namespace xtp {

        void GPUNumericalIntegrationTest::Initialize(Property* options) {
            string key = "options." + Identify();

            _orbfile = options->get(key + ".input").as<string> ();

            // get the path to the shared folders with xml files
            char *votca_share = getenv("VOTCASHARE");
            if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

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

            gpu::GPUNumericalIntegration GpuNI; 

            return true;
        }

}}
