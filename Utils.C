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

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

namespace WWWUtils{
  void printJetP4s(const std::vector<LorentzVector> & vecs){
    cout<<"Printing "<<vecs.size()<<" Jets:"<<endl;
    for (int i = 0; i<(int) vecs.size(); i++){
      cout<<" p4 with pt: "<<vecs.at(i).pt()<<" eta: "<<vecs.at(i).eta()<<" phi: "<<vecs.at(i).phi()<<endl;
    }
  }

  string printLep(short ind, bool verbose = false){
    std::stringstream ret;
    LorentzVector a = phys.lep_p4().at(ind);
    
    ret<<"(pt, eta, phi) = ("<<a.pt()<<", "<<a.eta()<<", "<<a.phi()<<") pass (cutbased_veto, cutbase_tight, mvabased_tight) = ("<<phys.lep_pass_VVV_cutbased_veto().at(ind)<<", "<<phys.lep_pass_VVV_cutbased_tight().at(ind)<<", "<<phys.lep_pass_VVV_MVAbased_tight().at(ind)<<") ";

    if (verbose){
      ret<<" pass ( [cutbased] veto_noiso, veto_noiso_noip, fo, fo_noiso, tight_noiso, [mvabased] tight, tight_noiso, [baseline]) ("<<phys.lep_pass_VVV_cutbased_veto_noiso().at(ind)<<", "<<phys.lep_pass_VVV_cutbased_veto_noiso_noip().at(ind)<<", "<<phys.lep_pass_VVV_cutbased_fo().at(ind)<<", "<<phys.lep_pass_VVV_cutbased_fo_noiso().at(ind)<<", "<<phys.lep_pass_VVV_cutbased_tight_noiso().at(ind)<<", "<<phys.lep_pass_VVV_MVAbased_tight_noiso().at(ind)<<", "<<phys.lep_pass_VVV_baseline().at(ind)<<")";
    }

    return ret.str();
  }
}

#endif