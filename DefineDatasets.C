#include "TChain.h"
#include "TString.h"
#include "ConfigHelper.C"

void addToChain(TChain *ch, TString set, bool hadoop, bool skimmed) {

  cout<<"Hadoop: "<<hadoop<<endl;
  cout<<"Skimmed: "<<skimmed<<endl;

  TString base_path;
  if (hadoop){
    base_path="/hadoop/cms/store/user/bhashemi/AutoTwopler_babies/merged/ZMET/";
  }
  else {
    base_path="/nfs-7/userdata/bhashemi/WWW_babies/";
  }

  TString version="WWW_v0.1.11";

  TString dir;
  if (skimmed){
    dir=base_path+version+"/skim/";
  }
  else{
    dir=base_path+version+"/output/";
  }

  cout<<"base_path: "<<base_path<<endl;
  
  //=============
  // V+Jets
  //=============

  if (set == "DY"){
    cout<<"Adding DY sample"<<endl;
    
    ch->Add(dir+"dy_m50_mgmlm_ext1*");
    ch->Add(dir+"dy_m50_mgmlm_ht100_ext1*");
    ch->Add(dir+"dy_m50_mgmlm_ht200_ext1*");
    ch->Add(dir+"dy_m50_mgmlm_ht400_ext1*"); 
    ch->Add(dir+"dy_m50_mgmlm_ht600_nonext*");
    ch->Add(dir+"dy_m50_mgmlm_ht800_nonext*");
    ch->Add(dir+"dy_m50_mgmlm_ht1200_nonext*");
    ch->Add(dir+"dy_m50_mgmlm_ht2500_nonext*");
    ch->Add(dir+"dy_m1050_mgmlm*.root");
  }

  else if (set == "WJets"){
    cout<<"Adding WJets"<<endl; 
    
    ch->Add(dir+"wjets_incl_mgmlm_*.root");
    ch->Add(dir+"wjets_ht100_mgmlm_ext1_*.root");
    ch->Add(dir+"wjets_ht200_mgmlm_ext1_*.root");
    ch->Add(dir+"wjets_ht400_mgmlm_ext1_*.root");
    ch->Add(dir+"wjets_ht600_mgmlm_ext1_*.root");
    ch->Add(dir+"wjets_ht800_mgmlm_ext1_*.root");
    ch->Add(dir+"wjets_ht1200_mgmlm_nonext_*.root");
    ch->Add(dir+"wjets_ht2500_mgmlm_ext1_*.root");
  }
  
  //=============
  // TTBar
  //=============

  else if (set == "TTBar2l"){
    cout<<"Adding TTBar to 2lep (smaller stats sample)"<<endl; 
    
    ch->Add(dir+"ttbar_dilep_mgmlm_ext1*");
  }
  else if (set == "TTBar2l-highstats"){
    cout<<"Adding TTBar to 2lep (large stats sample)"<<endl; 
    
    ch->Add(dir+"ttbar_dilep_powheg*"); //larger stats
  }
  else if (set == "TTBar1l"){
    cout<<"Adding TTBar to 1 lep"<<endl; 
    
    ch->Add(dir+"ttbar_1ltop_mgmlm_ext1*");
    ch->Add(dir+"ttbar_1ltbr_mgmlm_ext1*");
  }

  //=============
  // t/tV
  //=============
  
  else if (set == "TW-inclusive"){
    cout<<"Adding TW with inclusive decays"<<endl; 
    
    ch->Add(dir+"sttw_antitop_powheg*");
    ch->Add(dir+"sttw_top_powheg*");
  }
  else if (set == "SingleTop"){
    cout<<"Adding Single Top t and s channel"<<endl;  
    
    ch->Add(dir+"stt_antitop_incdec_powheg*");
    ch->Add(dir+"stt_top_incdec_powheg*");
    ch->Add(dir+"sts_4f_leptonic_amcnlo*");
  }
  else if (set == "TW"){
    cout<<"Adding TW with no fully hadronic decay"<<endl;
    
    ch->Add(dir+"sttw_antitop_nofullhaddecay_powheg*");
    ch->Add(dir+"sttw_top_nofullhaddecay_powheg*"); 
  }
  else if (set == "TZq"){
    cout<<"Adding TZq"<<endl; 
    
    ch->Add(dir+"tzq_ll_amcnlo*");
  }
  else if (set == "TWZ"){
    cout<<"Adding TWZ"<<endl; 
    
    ch->Add(dir+"sttwll_madgraph*");
  }
  //cout<<"Entries: "<<ch_fs->GetEntries()<<endl;

  //=============
  // ttV
  //=============

  else if (set == "TTW"){
    cout<<"Adding TTW MC"<<endl; 
    
    ch->Add(dir+"ttw_ln_amcnlo*");
    ch->Add(dir+"ttw_qq_amcnlo*");
  }
  else if (set == "TTZ"){
    cout<<"Adding TTZ"<<endl; 
    
    //ch->Add(dir+"ttz_2l2n_amcnlo*");
    ch->Add(dir+"ttz_incl_mgmlm*");
  }
  else if (set == "TTH"){
    cout<<"Adding TTH"<<endl; 
    
    ch->Add(dir+"tth_bb_powheg*");
    ch->Add(dir+"tth_nonbb_powheg*");
  }
  else if (set == "TTG"){
    cout<<"Adding TTG"<<endl; 
    
    ch->Add(dir+"ttg_incl_amcnlo*");
  }
  else if (set == "TTW-inclusive"){
    cout<<"Adding TTW inclusive MC"<<endl; 
    
    ch->Add(dir+"ttw_incl_mgmlm*");
  } 

  //=============
  // VV
  //=============

  else if (set == "WW"){
    cout<<"Adding WW"<<endl; 
    
    ch->Add(dir+"ww_2l2nu_powheg*.root");
    ch->Add(dir+"ww_lnuqq_powheg*.root");
    ch->Add(dir+"wpwpjj_ewk-qcd_madgraph*");
    ch->Add(dir+"ww_2l2nu_dbl_scat*");
  }
  else if (set == "WZ"){
    cout<<"Adding WZ MC"<<endl; 
    
    ch->Add(dir+"wz_3lnu_powheg*");
    ch->Add(dir+"wz_1l3n_amcnlo*");
    ch->Add(dir+"wz_lnqq_amcnlo*");
  }
  else if (set == "ZZ"){
    cout<<"Adding ZZ"<<endl; 
    
    ch->Add(dir+"zz_2l2n_powheg*");
    ch->Add(dir+"zz_2l2q_powheg*");
    ch->Add(dir+"zz_2q2n_powheg*");
    ch->Add(dir+"zz_4l_powheg*");
  }
  else if (set == "VH"){
    cout<<"Adding VH"<<endl; 

    ch->Add(dir+"vh_nonbb_amcnlo*");
  }
  //=============
  // VVV
  //=============

  else if (set == "WWW"){
    cout<<"Adding WWW"<<endl; 
    
    ch->Add(dir+"www_incl_amcnlo*"); 
  }
  else if (set == "WWW-mia"){
    cout<<"Adding WWW (Mia's samples)"<<endl; 
    
    ch->Add(dir+"www_2l_ext1_mia*"); 
    ch->Add(dir+"www_2l_mia*"); 
  }
  else if (set == "WWZ"){
    cout<<"Adding WWZ"<<endl; 
    
    ch->Add(dir+"wwz_incl_amcnlo*"); 
  }
  else if (set == "WZZ"){
    cout<<"Adding WZZ"<<endl; 
    
    ch->Add(dir+"wzz_incl_amcnlo*"); 
  }
  else if (set == "ZZZ"){
    cout<<"Adding ZZZ"<<endl; 
    
    ch->Add(dir+"zzz_incl_amcnlo*"); 
  }

//====================================
// QCD
//====================================
  else if (set == "QCD"){
    cout<<"Adding QCD (non-ext samples)"<<endl; 

    dir="/hadoop/cms/store/user/bhashemi/AutoTwopler_babies/merged/VVV/WWW_v0.1.11/";

    ch->Add(dir+"qcdht100_nonext*");
    ch->Add(dir+"qcdht200_nonext*");
    ch->Add(dir+"qcdht300_nonext*");
    ch->Add(dir+"qcdht500_nonext*");
    ch->Add(dir+"qcdht700_nonext*");
    ch->Add(dir+"qcdht1000_nonext*");
    ch->Add(dir+"qcdht1500_nonext*");
    ch->Add(dir+"qcdht2000_nonext*");
  }

  else if (set == "QCD-ext"){
    cout<<"Adding QCD (ext samples)"<<endl; 

    ch->Add(dir+"qcdht100_nonext*");
    ch->Add(dir+"qcdht200_ext1*");
    ch->Add(dir+"qcdht300_ext1*");
    ch->Add(dir+"qcdht500_ext1*");
    ch->Add(dir+"qcdht700_ext1*");
    ch->Add(dir+"qcdht1000_ext1*");
    ch->Add(dir+"qcdht1500_ext1*");
    ch->Add(dir+"qcdht2000_ext1*");
  }

//====================================
// Leptonic Data
//====================================
  
  version="WWW_v0.1.11";
  if (skimmed){
    dir=base_path+version+"/skim/";
  }
  else{
    dir=base_path+version+"/output/";
  }

  if (set == "Data-EE"){
    cout<<"Adding EE Trigger Data"<<endl;

    ch->Add(dir+"data_Run2016B_03feb2017rereco_ee_v2*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_ee_v1*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_ee_v1*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_ee_v1*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_ee_v1*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_ee_v1*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_ee_v2*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_ee_v3*");
  }
  else if (set == "Data-EMu"){  
    cout<<"Adding EMu Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_03feb2017rereco_em_v2*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_em_v1*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_em_v1*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_em_v1*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_em_v1*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_em_v1*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_em_v2*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_em_v3*");
  }
  else if (set == "Data-MuMu"){
    cout<<"Adding MuMu Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_03feb2017rereco_mm_v2*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_mm_v1*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_mm_v1*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_mm_v1*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_mm_v1*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_mm_v1*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_mm_v2*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_mm_v3*");
  }
  else if (set == "Data-Mu"){
    cout<<"Adding Single Muon Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_03feb2017rereco_sm_v2*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_sm_v1*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_sm_v1*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_sm_v1*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_sm_v1*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_sm_v1*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_sm_v2*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_sm_v3*");
  }

//====================================
// End Building TChain
//====================================

  return;
} 


TChain * getTChain(TString data_set, bool hadoop=false, bool skimmed=true) {
  TChain *ch = new TChain("t");

  vector<TString> sets;

  if (data_set[0]=='['){
     sets = sParseVector(data_set);
  }
  else{
    TString token;
    Ssiz_t from=0;
    //cout<<"got vector in string form: "<<opt<<endl;
    while(data_set.Tokenize(token, from, "_")){
      token.ReplaceAll("_", "");
      sets.push_back(TString(token.Data()));
    }
  }

  cout<<"Datasets Incoming: "<<endl;
  cout<<"===================================="<<endl;
  for (std::vector<TString>::iterator i=sets.begin(); i != sets.end(); i++){
    addToChain(ch, *i, hadoop, skimmed);
  }
  cout<<"===================================="<<endl;
  
  return ch;
}