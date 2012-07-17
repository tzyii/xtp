#ifndef PAIRDUMP_H
#define PAIRDUMP_H

#include <votca/ctp/qmcalculator.h>
#include <sys/stat.h>


namespace votca { namespace ctp {

class PairDump : public QMCalculator
{
public:

    PairDump() {};
   ~PairDump() {};

    string  Identify() { return "PairDump"; }

    void    Initialize(Topology *top, Property *options);
    bool    EvaluateFrame(Topology *top);

private:

    string _outParent;
    string _outMonDir;
    string _outDimDir;

    string _outFormat;
    bool   _subFolder;
};


void PairDump::Initialize(Topology *top, Property *options) {

    string key = "options.pairdump";   
}


bool PairDump::EvaluateFrame(Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }

    FILE *out;

    _outParent = "frame" + boost::lexical_cast<string>(top->getDatabaseId());
    mkdir(_outParent.c_str(), 0755);

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {

        string ID   = boost::lexical_cast<string>((*sit)->getId());
        string DIR  = _outParent + "/mol_" + ID;
        string FILE = _outParent + "/mol_" + ID + "/mol_" + ID + ".xyz";

        mkdir(DIR.c_str(), 0755);        
        out = fopen(FILE.c_str(),"w");
        (*sit)->WriteXYZ(out);
        fclose(out);
    }

    QMNBList ::iterator pit;
    QMNBList &nblist = top->NBList();
    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {

        string ID1  = boost::lexical_cast<string>((*pit)->Seg1()->getId());
        string ID2  = boost::lexical_cast<string>((*pit)->Seg2()->getId());
        string DIR1 = _outParent + "/pair_" + ID1 + "_" + ID2;
        string DIR2 = DIR1 + "/dim";
        string FILE = DIR2 + "/pair_" + ID1 + "_" + ID2 + ".xyz";

        mkdir(DIR1.c_str(), 0755);
        mkdir(DIR2.c_str(), 0755);
        out = fopen(FILE.c_str(),"w");
        (*pit)->WriteXYZ(out);
        fclose(out);
    }
}









}}

#endif
