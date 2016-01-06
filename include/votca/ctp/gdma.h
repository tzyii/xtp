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

#ifndef __VOTCA_CTP_GDMA_H
#define	__VOTCA_CTP_GDMA_H

#include <string>
#include <map>
#include <votca/tools/property.h>
#include <fstream>




namespace votca { namespace ctp {


using namespace std;
using namespace votca::tools;

/**
    \brief information about an element
 
    The Atom class stores atom id, name, type, mass, charge, residue number
    
*/
class GDMA
{
public:   

  

    GDMA() { };
   ~GDMA() { };

   void WriteInputFile( );
   void RunExternal();
   void ParseOutputFile();
   void Initialize( Property *options  ); 
   
   // functions to override Initialize
   void SetLimit( double rank  ) { _limit = rank; } ;
   void SetRunDir( string dir ) { _runFolder = dir; }
   void SetChkFile (string file ) { _chkFile = file; };
   void SetExecutable( string exec ){ _executable = exec;};
   void SetDensity( string density ) { _density = density;};
   void SetRadius( double radius ) { _radius = radius;};
   void SetSwitch( double sw ) { _switch = sw; };
   
   void setLog( Logger* pLog ) { _log = pLog; };
   vector< vector<double> > &GetMultipoles() { return _multipoles; };
   
   
private:

    
    vector< vector<double> > _multipoles;
    
    string _runFolder;
    string _chkFile; 
    string _executable;

    string  _density; 
    int     _limit; 
    double  _radius; 
    double  _switch; 
    string  _outFile; 
    Logger*                             _log;
    
 
};
}}

#endif	/* __VOTCA_CTP_GDMA_H */
