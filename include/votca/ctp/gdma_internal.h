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

#ifndef __VOTCA_CTP_GDMA_INTERNAL_H
#define	__VOTCA_CTP_GDMA_INTERNAL_H

#include <string>
#include <map>
#include <votca/tools/property.h>
#include <votca/ctp/apolarsite.h>
#include <fstream>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/polarseg.h>



namespace votca { namespace ctp {


using namespace votca::tools;


class GDMA_internal
{
public:   

  

    GDMA_internal() { };
   ~GDMA_internal() { };
   
   int Run_gdma( Orbitals *orbitals);

private:
    
    int state;
    
    bool do_transition;
    PolarSeg sites;
    
};
}}

#endif	/* __VOTCA_CTP_GDMA_H */

