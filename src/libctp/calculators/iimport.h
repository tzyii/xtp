#ifndef IIMPORT_H
#define IIMPORT_H

#include <votca/ctp/qmcalculator.h>
#include <sys/stat.h>
#include <votca/ctp/logger.h>

namespace votca { namespace ctp {

class IImport : public QMCalculator
{
public:

    string      Identify() { return "iimport"; }

    void        Initialize(Topology *top, Property *options);
    bool        EvaluateFrame(Topology *top);
    void        XML2PairTI(QMPair *qmpair, string &xmlDirFile);
    void        List2PairsTI(Topology *top, string &ti_file);
    void        FromIDFT(Topology *top, string &_idft_jobs_file);
    void        FromIDFTWithSuperExchange(Topology *top, string &_idft_jobs_file);

private:

    bool        _importFromDirs;
    bool        _importFromList;
    bool        _importFromIDFT;

    string      _TI_tag;
    string      _TI_file;
    string      _idft_jobs_file;
};


void IImport::Initialize(Topology *top, Property *options) {
    
    _importFromDirs = false;
    _importFromList = false;
    _importFromIDFT = false;

    // _options already has default values, update them with the supplied options
    _options.CopyValues("", *options );
    
    string key = "options." + Identify();

    //if (options->exists(key+".TI_tag")) {
    _TI_tag = _options.get(key+".TI_tag").as< string >();
    if (_TI_tag != "") {
        _importFromDirs = true;
        cout << endl << "... ... Using TI XML tag '" << _TI_tag << "'" << flush;
    }

    //else if (options->exists(key+".TI_file")) {
    
    _TI_file = _options.get(key+".TI_file").as< string >();
    if (_TI_file != "") {
        _importFromList = true;
        cout << endl << "... ... Using TI table '" << _TI_file << "'" << flush;
    }
    
    _idft_jobs_file = _options.get(key+".idft_jobs_file").as< string >();
    if ( _idft_jobs_file != "" ) {
        _importFromIDFT = true;
        cout << endl << "... ... Using IDFT jobs file '" << _TI_file << "'" << flush;
    }

    
}

bool IImport::EvaluateFrame(Topology *top) {
  QMNBList &nblist = top->NBList();
  QMNBList ::iterator nit;
    
  // Import from file
  if (_importFromList) {
      this->List2PairsTI(top, _TI_file);
  }

  // Import from DiPro folders
  else if (_importFromDirs) {
    string PARDIR = "frame" + boost::lexical_cast<string>(top->getDatabaseId());
    string SUBDIR = PARDIR + "/pair_0000_0000";



    for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
        string ID1  = boost::lexical_cast<string>((*nit)->Seg1()->getId());
        string ID2  = boost::lexical_cast<string>((*nit)->Seg2()->getId());

        SUBDIR = PARDIR + "/pair_" + ID1 + "_" + ID2 + "/TI.xml";
        this->XML2PairTI(*nit, SUBDIR);
    }
  }

  if ( _importFromIDFT ) FromIDFT(top, _idft_jobs_file);
  
  if ( _importFromIDFT && nblist.getSuperExchangeTypes().size() != 0 ) FromIDFTWithSuperExchange(top, _idft_jobs_file);
  
}


void IImport::List2PairsTI(Topology *top, string &ti_file) {

    QMNBList &nblist = top->NBList();
    int pair_count = 0;

    std::string line;
    std::ifstream intt;
    intt.open(ti_file.c_str());

    if (intt.is_open()) {
        while (intt.good()) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " \t");
            toker.ToVector(split);

            if (split.size() == 0) { continue; }
            
            if (split.size() != 4 && split.size() != 6) {
                cout << endl << "... ... Invalid line: " << line << flush;
                continue;
            }

            int seg1id = boost::lexical_cast<int>(split[0]);
            int seg2id = boost::lexical_cast<int>(split[1]);

            Segment *seg1 = top->getSegment(seg1id);
            Segment *seg2 = top->getSegment(seg2id);

            QMPair *qmp = nblist.FindPair(seg1,seg2);
            if (qmp == NULL) {
                cout << endl
                     << "... ... ERROR: " << line
                     << flush;
                cout << endl
                     << "... ... Line is not compatible with neighborlist. "
                     << flush;
                throw std::runtime_error("Forgot to run -e neighborlist?");
            }            

            int e_h_1 = (split[2] == "e") ? -1 : +1;
            double j2_1 = boost::lexical_cast<double>(split[3]);

            if (split.size() == 6) {

                int e_h_2 = (split[4] == "e") ? -1 : +1;
                if (e_h_1 == e_h_2) {
                    cout << endl << "... ... Invalid line: " << line << flush;
                    continue;
                }
                double j2_2 = boost::lexical_cast<double>(split[5]);

                qmp->setJeff2(j2_2, e_h_2);
                qmp->setIsPathCarrier(true, e_h_2);
            }

            qmp->setJeff2(j2_1, e_h_1);
            qmp->setIsPathCarrier(true, e_h_1);
            ++pair_count;
            
        }
    }
    else {
        cout << endl
             << "ERROR: No such file" << ti_file << ". "
             << flush;
        throw std::runtime_error("Missing input file in IImport calculator.");
    }

    cout << endl
         << "... ... Set transfer integrals for " << pair_count << " pairs. "
         << flush;
}


void IImport::XML2PairTI(QMPair *qmpair, string &xmlDirFile) {

    printf("\n... ... Import TIs for pair %5d (ID1 %4d    ID2 %4d) ",
            qmpair->getId(), qmpair->Seg1()->getId(), qmpair->Seg2()->getId());

    string TRANSPORT = "electron_or_hole";
    int    STATE     = 0;
    
    std::string line;
    std::ifstream intt;
    intt.open(xmlDirFile.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " \t<>=\"");
            toker.ToVector(split);

            if (split.size() < 3) {
                continue;
            }

//            cout << split.size() << " " << flush;
//            vector<string> ::iterator sit;
//            for (sit = split.begin(); sit < split.end(); ++sit) {
//                cout << (*sit) << " " << flush;
//            }

            if      (split[0] == "transport") {
                assert(split[1] == "name");
                TRANSPORT = split[2];
            }
            else if (split[0] == _TI_tag) {
                assert(TRANSPORT == "electron" || TRANSPORT == "hole");
                STATE = (TRANSPORT == "electron") ? -1 : +1;
                if (_TI_tag == "J" || _TI_tag == "T_00") {
                    double J = boost::lexical_cast<double>(split[1]);
                    qmpair->setJeff2(J*J, STATE);
                    qmpair->setIsPathCarrier(1, STATE);
                    printf("\n... ... ... J2(State = %+1d) = %4.7e",
                            STATE, qmpair->getJeff2(STATE));
                }
                else if (_TI_tag == "J_sq_degen" || _TI_tag == "J_sq_boltz") {
                    double J2 = boost::lexical_cast<double>(split[1]);
                    qmpair->setJeff2(J2, STATE);
                    qmpair->setIsPathCarrier(1, STATE);
                    printf("\n... ... ... J2(State = %+1d) = %4.7e",
                            STATE, qmpair->getJeff2(STATE));
                }
            }
        } /* Exit loop over lines */
    }
    else { cout << endl
                << "... ... ... ERROR: No TI.xml in pair folder. Skip..."
                << flush;
    }

    intt.close();

}


void IImport::FromIDFT(Topology *top, string &_idft_jobs_file) {

    Property xml;

    QMNBList &nblist = top->NBList();   
    int _number_of_pairs = nblist.size();
    int _current_pairs = 0;
    int _incomplete_jobs = 0;
    
    Logger _log;
    _log.setReportLevel(logINFO);
    
    load_property_from_xml(xml, _idft_jobs_file);
    
    list<Property*> jobProps = xml.Select("jobs.job");
    list<Property*> ::iterator it;

    for (it = jobProps.begin(); it != jobProps.end(); ++it) {
 
        if ( (*it)->exists("output") && (*it)->exists("output.pair") ) {
            
            //cout << **it;
            
            Property poutput = (*it)->get("output.pair");
            
            int homoA = poutput.getAttribute<int>("homoA");
            int homoB = poutput.getAttribute<int>("homoB");
            
            int idA = poutput.getAttribute<int>("idA");
            int idB = poutput.getAttribute<int>("idB");
                       
            string typeA = poutput.getAttribute<string>("typeA");
            string typeB = poutput.getAttribute<string>("typeB");
            //cout << idA << ":" << idB << "\n"; 
            Segment *segA = top->getSegment(idA);
            Segment *segB = top->getSegment(idB);

            QMPair *qmp = nblist.FindPair(segA,segB);
            
            if (qmp == NULL) { // there is no pair in the neighbor list with this name
                //LOG(logINFO, _log) << "No pair " <<  idA << ":" << idB << " found in the neighbor list. Ignoring" << flush; 
            }   else {
                
                _current_pairs++;
                
                list<Property*> pOverlap = poutput.Select("overlap");
                list<Property*> ::iterator itOverlap;

                    
                    for (itOverlap = pOverlap.begin(); itOverlap != pOverlap.end(); ++itOverlap) {

                        double energyA = (*itOverlap)->getAttribute<double>("eA");
                        double energyB = (*itOverlap)->getAttribute<double>("eB");
                        double overlapAB = (*itOverlap)->getAttribute<double>("jAB");
                        int orbA = (*itOverlap)->getAttribute<double>("orbA");
                        int orbB = (*itOverlap)->getAttribute<double>("orbB");

                        if ( orbA == homoA && orbB == homoB ) {
                                qmp->setJeff2(overlapAB*overlapAB, 1);
                                qmp->setIsPathCarrier(true, 1);
                        }

                        if ( orbA == homoA+1 && orbB == homoB+1 ) {
                                qmp->setJeff2(overlapAB*overlapAB, -1);
                                qmp->setIsPathCarrier(true, -1);
                        }
                    }    
            }
            
        } else {
            _incomplete_jobs++;
            LOG(logINFO, _log) << "Job " << (*it)->get( "id" ).as<string>() << " is " << (*it)->get( "status" ).as<string>() << endl;
        }
    }
    
    LOG(logINFO, _log) << "Pairs [total:updated] " <<  _number_of_pairs << ":" << _current_pairs << " Incomplete jobs: " << _incomplete_jobs << flush; 
    cout << _log;
}

    
void IImport::FromIDFTWithSuperExchange(Topology *top, string &_idft_jobs_file) {

    Property xml;

    vector<Property*> records;
            
    QMNBList &nblist = top->NBList();
    int _number_of_pairs = nblist.size();
    int _current_pairs = 0;
    int _incomplete_jobs = 0;
    
    Logger _log;
    _log.setReportLevel(logINFO);
    
    //generate lists of bridges for superexchange pairs
    nblist.GenerateSuperExchange();

    // load the QC results in a vector indexed by the pair ID
    load_property_from_xml(xml, _idft_jobs_file);
    list<Property*> jobProps = xml.Select("jobs.job");
    
    records.resize( jobProps.size() + 1  );
    
    for (list<Property*> ::iterator  it = jobProps.begin(); it != jobProps.end(); ++it) {
 
        if ( (*it)->exists("output") && (*it)->exists("output.pair") ) {
            
            Property poutput = (*it)->get("output.pair");
            
            int idA = poutput.getAttribute<int>("idA");
            int idB = poutput.getAttribute<int>("idB");
                       
            Segment *segA = top->getSegment(idA);
            Segment *segB = top->getSegment(idB);

            QMPair *qmp = nblist.FindPair(segA,segB);
            
            if (qmp == NULL) { // there is no pair in the neighbor list with this name
                ;//LOG(logINFO, _log) << "No pair " <<  idA << ":" << idB << " found in the neighbor list. Ignoring" << flush; 
            }   else {
                LOG(logINFO, _log) << "Store in record: " <<  idA << ":" << idB << flush; 
                records[qmp->getId()] = & ((*it)->get("output.pair"));
            }
        }
    } // finished loading from the file


    // loop over all pairs in the neighborlist
    for (QMNBList::iterator ipair = top->NBList().begin(); ipair != top->NBList().end(); ++ipair) {
        
        QMPair *pair = *ipair;
        Segment* segmentA = pair->Seg1PbCopy();
        Segment* segmentB = pair->Seg2PbCopy();
        
        double Jeff2_homo = 0;
        double Jeff2_lumo = 0;
        
        cout << "Processing pair " << segmentA->getId() << ":" << segmentB->getId() << endl;
        
        QMPair::PairType _ptype = pair->getType();
        Property* pair_property = records[ pair->getId() ];
 
        int homoA = pair_property->getAttribute<int>("homoA");
        int homoB = pair_property->getAttribute<int>("homoB");
       
        // If a pair is of a direct type 
        if ( _ptype == QMPair::Hopping ||  _ptype == QMPair::SuperExchangeAndHopping ) {
            cout << "Pair is hopping" << endl;
            list<Property*> pOverlap = pair_property->Select("overlap");
 
            for (list<Property*> ::iterator itOverlap = pOverlap.begin(); itOverlap != pOverlap.end(); ++itOverlap) {

                double overlapAB = (*itOverlap)->getAttribute<double>("jAB");
                int orbA = (*itOverlap)->getAttribute<double>("orbA");
                int orbB = (*itOverlap)->getAttribute<double>("orbB");

                if ( orbA == homoA && orbB == homoB ) {
                    Jeff2_homo += overlapAB*overlapAB;
                }

                if ( orbA == homoA+1 && orbB == homoB+1 ) {
                    Jeff2_lumo += overlapAB*overlapAB;
                }
            }    
            
        }
        
        // if pair has bridges only
        if ( _ptype == QMPair::SuperExchange  ||  _ptype == QMPair::SuperExchangeAndHopping ) {
            
            list<Property*> pOverlap = pair_property->Select("overlap");
            
            // this is to select HOMO_A and HOMO_B 
            double overlapAB;
            int orbA;
            int orbB;
            double energyA;
            double energyB;
            
            for (list<Property*> ::iterator itOverlap = pOverlap.begin(); itOverlap != pOverlap.end(); ++itOverlap) {
                if ( orbA == homoA && orbB == homoB ) {  
                    overlapAB = (*itOverlap)->getAttribute<double>("jAB");
                    orbA = (*itOverlap)->getAttribute<double>("orbA");
                    orbB = (*itOverlap)->getAttribute<double>("orbB");
                    energyA = (*itOverlap)->getAttribute<double>("eA");
                    energyB = (*itOverlap)->getAttribute<double>("eB");
                }
            }
            
            
            
            // loop over the bridging segments
            for ( vector< Segment* >::iterator itBridge = pair->getBridgingSegments().begin() ; itBridge != pair->getBridgingSegments().end(); itBridge++ ) {

                Segment* Bridge = *itBridge;
                int IDBridge = Bridge->getId();

                // pairs from the bridge to the donor and acceptor
                QMPair* Bridge_A = nblist.FindPair( segmentA, Bridge );
                QMPair* Bridge_B = nblist.FindPair( segmentB, Bridge );

                Property* pBridge_A = records[ Bridge_A->getId() ];
                Property* pBridge_B = records[ Bridge_B->getId() ];

                list<Property*> pOverlapA = pBridge_A->Select("overlap");
                list<Property*> pOverlapB = pBridge_B->Select("overlap");

                // IDs of the Donor and Acceptor
                int IdA = segmentA->getId();
                int IdB = segmentB->getId();

                // IDs stored in the file
                int id1A = pBridge_A->getAttribute<int>("idA");
                int id2A = pBridge_A->getAttribute<int>("idB");

                int id1B = pBridge_B->getAttribute<int>("idA");
                int id2B = pBridge_B->getAttribute<int>("idB");

                // suffix for the donor and acceptor
                string suffixA = ( id1A == IDBridge ) ? "B" : "A"; // use "A" as a bridge 
                string suffixB = ( id1B == IDBridge ) ? "B" : "A"; // use "A" as a bridge 
                string suffixBridgeA = ( id1A == IDBridge ) ? "A" : "B";
                string suffixBridgeB = ( id1B == IDBridge ) ? "A" : "B";
                
                int homoBridgeA = pBridge_A->getAttribute<int>("orb" + suffixBridgeA );
                int homoBridgeB = pBridge_B->getAttribute<int>("orb" + suffixBridgeB );
                assert( homoBridgeA == homoBridgeB );
                int homoBridge = homoBridgeA;
               
                for (list<Property*> ::iterator itOverlapA = pOverlapA.begin(); itOverlapA != pOverlapA.end(); ++itOverlapA) {
                for (list<Property*> ::iterator itOverlapB = pOverlapB.begin(); itOverlapB != pOverlapB.end(); ++itOverlapB) {
                    
                    int orbDonor = (*itOverlapA)->getAttribute<int>( "orb" + suffixA );
                    int orbAcceptor = (*itOverlapB)->getAttribute<int>( "orb" + suffixB );
                    int orbBridgeA  = (*itOverlapA)->getAttribute<int>( "orb" + suffixBridgeA );
                    int orbBridgeB = (*itOverlapB)->getAttribute<int>( "orb" + suffixBridgeB );
                    
                    if (  orbDonor == homoA && orbAcceptor == homoB && orbBridgeA == orbBridgeB && orbBridgeA <= homoBridge) {
                        
                        double jDB = (*itOverlapA)->getAttribute<double>( "jAB" );
                        double jBA = (*itOverlapB)->getAttribute<double>( "jAB" );
                        double eA  = (*itOverlapA)->getAttribute<double>( "e" + suffixA );
                        double eB  = (*itOverlapB)->getAttribute<double>( "e" + suffixB );
                        
                        double eBridgeA  = (*itOverlapA)->getAttribute<double>( "e" + suffixBridgeA );
                        double eBridgeB  = (*itOverlapB)->getAttribute<double>( "e" + suffixBridgeB );
                        
                        //assert( eBridgeA - eBridgeB < 1e-50 );
                     
                        cout << homoA << " " << homoB << " " << (*itOverlapA)->getAttribute<int>( "orb" + suffixBridgeA )
                             << " JDB " << jDB 
                             << " JBA " << jBA << endl;
                        
                        // This in principle violates detailed balance. Any ideas?
                        Jeff2_homo += 0.5 * (jDB*jBA / (eA - eBridgeA) + jDB*jBA / (eB - eBridgeB));
                        
                                
                    }

                    if (  orbDonor == homoA+1 && orbAcceptor == homoB+1 && orbBridgeA == orbBridgeB && orbBridgeA > homoBridge) {
                        
                        double jDB = (*itOverlapA)->getAttribute<double>( "jAB" );
                        double jBA = (*itOverlapB)->getAttribute<double>( "jAB" );
                        double eA  = (*itOverlapA)->getAttribute<double>( "e" + suffixA );
                        double eB  = (*itOverlapB)->getAttribute<double>( "e" + suffixB );
                        
                        double eBridgeA  = (*itOverlapA)->getAttribute<double>( "e" + suffixBridgeA );
                        double eBridgeB  = (*itOverlapB)->getAttribute<double>( "e" + suffixBridgeB );
                        
                         // This in principle violates detailed balance. Any ideas?
                        Jeff2_lumo += 0.5 * (jDB*jBA / (eA - eBridgeA) + jDB*jBA / (eB - eBridgeB));
                        //jDB*jBA / (eB - eBridgeB);
                                
                    }
                    
                     
                }}
            } // end over bridges 
            
            
            
        } // end of if superexchange
         
        pair->setJeff2(Jeff2_homo, 1);
        pair->setIsPathCarrier(true, 1);
        
        pair->setJeff2(Jeff2_lumo, -1);
        pair->setIsPathCarrier(true, -1);
       
        break;
    }
                    
    LOG(logINFO, _log) << "Pairs [total:updated] " <<  _number_of_pairs << ":" << _current_pairs << " Incomplete jobs: " << _incomplete_jobs << flush; 
    cout << _log;
}
    
    
}}

#endif
