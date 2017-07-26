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

#ifndef _VOTCA_XTP_QMANALYZE_H
#define _VOTCA_XTP_QMANALYZE_H

#include <stdio.h>
#include <boost/format.hpp>
#include <votca/tools/constants.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/qmtool.h>
#include <votca/xtp/orbitals.h>
// #include <votca/xtp/mbgft.h>
// #include <votca/xtp/qmpackagefactory.h>

namespace votca {
    namespace xtp {
        using namespace std;
        namespace ub = boost::numeric::ublas;

        class QMAnalyze : public ctp::QMTool {
        public:

            QMAnalyze() {
            };

            ~QMAnalyze() {
            };

            string Identify() {
                return "qmanalyze";
            }

            void Initialize(tools::Property *options);
            bool Evaluate();


        private:

            string _orbfile;

            bool _print_BSE_singlets;
            bool _print_BSE_triplets;
            bool _print_GW_energies;
            bool _print_QP_energies;
            bool _print_DFT_energies;


            ctp::Logger _log;

            void CheckContent(Orbitals& _orbitals);

        };


    }
}


#endif
