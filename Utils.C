#ifndef INCLUDED_WWWUTILS
#define INCLUDED_WWWUTILS

// C++
#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <utility>
#include <fstream>


// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace WWWUtils{
  void printP4s(const std::vector<LorentzVector> & vecs, TString type="Jets"){
    cout<<"Printing "<<vecs.size()<<" "<<type<<endl;
    for (int i = 0; i<(int) vecs.size(); i++){
      cout<<" p4 with pt: "<<vecs.at(i).pt()<<" eta: "<<vecs.at(i).eta()<<" phi: "<<vecs.at(i).phi()<<endl;
    }
  }
}

#endif