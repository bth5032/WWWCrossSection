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
    base_path="/nfs-7/userdata/ZMEToutput/output/ZMETbabies/";
  }

  TString version="WWW_v0.1";

  TString dir;
  if (skimmed){
    dir=base_path+version+"/skims/";
  }
  else{
    dir=base_path+version+"/output/";
  }

  cout<<"base_path: "<<base_path<<endl;

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
  else if (set == "TTBar"){
    cout<<"Adding TTBar (smaller stats sample)"<<endl; 
    
    ch->Add(dir+"ttbar_dilep_mgmlm_ext1_*");
  }
  else if (set == "TTBar-highstats"){
    cout<<"Adding TTBar (large stats sample)"<<endl; 
    
    ch->Add(dir+"ttbar_dilep_powheg*"); //larger stats
  }
  else if (set == "TW-inclusive"){
    cout<<"Adding TW with inclusive decays"<<endl; 
    
    ch->Add(dir+"sttw_antitop_powheg*");
    ch->Add(dir+"sttw_top_powheg*");
  }
  else if (set == "TW"){
    cout<<"Adding TW with no fully hadronic decay"<<endl; 
    
    ch->Add(dir+"sttw_antitop_nofullhaddecay_powheg*");
    ch->Add(dir+"sttw_top_nofullhaddecay_powheg*");
  }
  //cout<<"Entries: "<<ch_fs->GetEntries()<<endl;
  else if (set == "WW"){
    cout<<"Adding WW"<<endl; 
    
    ch->Add(dir+"ww_2l2nu_powheg*.root");
  }
  else if (set == "TTW"){
    cout<<"Adding TTW MC"<<endl; 
    
    ch->Add(dir+"ttw_ln_amcnlo*");
    ch->Add(dir+"ttw_qq_amcnlo*");
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
    ch->Add(dir+"zz_4l_powheg*");
  }
  else if (set == "WWW"){
    cout<<"Adding WWW"<<endl; 
    
    ch->Add(dir+"www_incl_amcnlo*"); 
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
  else if (set == "TZq"){
    cout<<"Adding TZq"<<endl; 
    
    ch->Add(dir+"tzq_ll_amcnlo*");
  }
  else if (set == "TWZ"){
    cout<<"Adding TWZ"<<endl; 
    
    ch->Add(dir+"sttwll_madgraph*");
  }
  else if (set == "TTZ"){
    cout<<"Adding TTZ"<<endl; 
    
    //ch->Add(dir+"ttz_2l2n_amcnlo*");
    ch->Add(dir+"ttz_incl_mgmlm*");
  }
  else if (set == "TTW-inclusive"){
    cout<<"Adding TTW inclusive MC"<<endl; 
    
    ch->Add(dir+"ttw_incl_mgmlm*");
  } 

//====================================
// Photon Data
//====================================

  //Single Photon Trigger
  else if (set == "GammaData-SinglePhoton"){
    cout<<"Adding GammaData-SinglePhoton"<<endl; 
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-05/data_2016B_23sep2016rereco_ph_v3*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-05/data_2016C_23sep2016rereco_ph_v1*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-05/data_2016D_23sep2016rereco_ph_v1*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-05/data_2016E_23sep2016rereco_ph_v1*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-05/data_2016F_23sep2016rereco_ph_v1*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-05/data_2016G_23sep2016rereco_ph_v1*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-05/data_2016H_Prompt_ph_v2*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-05/data_2016H_Prompt_ph_v3*");
    
    //EWK Subtraction
  }
  else if (set == "GammaData-EWKSub"){
    cout<<"Adding EWK Subtraction Samples"<<endl;       
    //============
    // W+Gamma+Jets
    //============
    //This is the Wjets sample, it is intended to have events with a prompt photon vetod
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wjets_incl_mgmlm*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wjets_ht100_mgmlm*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wjets_ht200_mgmlm*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wjets_ht400_mgmlm*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wjets_ht600_mgmlm*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wjets_ht800_mgmlm*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wjets_ht1200_mgmlm*");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wjets_ht2500_mgmlm*");

    //This is the W+Gamma+Jets, it is inteded to have events with non-prompt photons vetod
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-04/wgjets_incl_mgmlm*");
    //============
    // Gamma+Z->NuNu
    //============
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/znunugamma_ptg130_mgmlm*.root");
    //============
    // TTbar
    //============
    //1lep
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/ttbar_1ltbr_mgmlm*.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/ttbar_1ltop_mgmlm*.root");
    //dilep
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/ttbar_dilep_mgmlm*.root");
    //============
    // Single Top
    //============
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/sttw_antitop_nofullhaddecay_powheg.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/sttw_top_nofullhaddecay_powheg.root");
  }
  else if (set == "GammaData-GammaRealMET"){
    cout<<"Adding GammaData-GammaRealMET"<<endl; 
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/ttbar_1ltbr_mgmlm.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/ttbar_1ltop_mgmlm.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wgjets_incl_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht100_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht1200_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht200_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht2500_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht400_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht600_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht800_amcnlo.root");
  }
  else if (set == "GammaData-JustWjetRealMET"){
    cout<<"Adding GammaData-JustWjetRealMET"<<endl; 
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht100_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht1200_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht200_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht2500_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht400_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht600_amcnlo.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wjets_ht800_amcnlo.root");
  }
  else if (set == "GammaData-JustWGjetRealMET"){
    cout<<"Adding GammaData-JustWGjetRealMET"<<endl; 
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/wgjets_incl_amcnlo.root");
  }
  else if (set == "GammaData-JustTTBarRealMET"){
    cout<<"Adding GammaData-JustTTBarRealMET"<<endl; 
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/ttbar_1ltbr_mgmlm.root");
    ch->Add("/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-11-10/ttbar_1ltop_mgmlm.root");
  }

//====================================
// Leptonic Data
//====================================

  else if (set == "DileptonData-ee"){
    TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-13/output/";
    
    cout<<"Adding EE Trigger Data"<<endl;

    ch->Add(dir+"data_Run2016B_23sep2016rereco_ee_v3*");
    ch->Add(dir+"data_Run2016C_23sep2016rereco_ee_v1*");
    ch->Add(dir+"data_Run2016D_23sep2016rereco_ee_v1*");
    ch->Add(dir+"data_Run2016E_23sep2016rereco_ee_v1*");
    ch->Add(dir+"data_Run2016F_23sep2016rereco_ee_v1*");
    ch->Add(dir+"data_Run2016G_23sep2016rereco_ee_v1*");
    ch->Add(dir+"data_Run2016H_Prompt_ee_v2*");
    ch->Add(dir+"data_Run2016H_Prompt_ee_v3*");
  }
  else if (set == "DileptonData-em"){
    TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-13/output/";
    
    cout<<"Adding EMu Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_23sep2016rereco_em_v3*");
    ch->Add(dir+"data_Run2016C_23sep2016rereco_em_v1*");
    ch->Add(dir+"data_Run2016D_23sep2016rereco_em_v1*");
    ch->Add(dir+"data_Run2016E_23sep2016rereco_em_v1*");
    ch->Add(dir+"data_Run2016F_23sep2016rereco_em_v1*");
    ch->Add(dir+"data_Run2016G_23sep2016rereco_em_v1*");
    ch->Add(dir+"data_Run2016H_Prompt_em_v2*");
    ch->Add(dir+"data_Run2016H_Prompt_em_v3*");
  }
  else if (set == "DileptonData-mm"){
    TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-13/output/";
    
    cout<<"Adding MuMu Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_23sep2016rereco_mm_v3*");
    ch->Add(dir+"data_Run2016C_23sep2016rereco_mm_v1*");
    ch->Add(dir+"data_Run2016D_23sep2016rereco_mm_v1*");
    ch->Add(dir+"data_Run2016E_23sep2016rereco_mm_v1*");
    ch->Add(dir+"data_Run2016F_23sep2016rereco_mm_v1*");
    ch->Add(dir+"data_Run2016G_23sep2016rereco_mm_v1*");
    ch->Add(dir+"data_Run2016H_Prompt_mm_v2*");
    ch->Add(dir+"data_Run2016H_Prompt_mm_v3*");
  }
  else if (set == "SingleLeptonData-SingleMu"){
    TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-13/output/";
    
    cout<<"Adding Single Muon Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_23sep2016rereco_sm_v3*");
    ch->Add(dir+"data_Run2016C_23sep2016rereco_sm_v1*");
    ch->Add(dir+"data_Run2016D_23sep2016rereco_sm_v1*");
    ch->Add(dir+"data_Run2016E_23sep2016rereco_sm_v1*");
    ch->Add(dir+"data_Run2016F_23sep2016rereco_sm_v1*");
    ch->Add(dir+"data_Run2016G_23sep2016rereco_sm_v1*");
    ch->Add(dir+"data_Run2016H_Prompt_sm_v2*");
    ch->Add(dir+"data_Run2016H_Prompt_sm_v3*");
  }

//====================================
// Skimmed Leptonic Data
//====================================

  else if (set == "DileptonData-ee-Skimmed"){
    //TString dir="/home/users/cwelke/ZMetbabyskims/V08-22-05/";
    //TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-11/skim/";
    TString dir="/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-15/skims/";

    cout<<"Adding Skimmed EE Trigger Data"<<endl;
    ch->Add(dir+"data_Run2016B_03feb2017rereco_ee_v2_skim*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_ee_v1_skim*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_ee_v1_skim*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_ee_v1_skim*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_ee_v1_skim*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_ee_v1_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_ee_v2_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_ee_v3_skim*");
  }
  else if (set == "DileptonData-em-Skimmed"){
    //TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-11/skim/";
    TString dir="/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-15/skims/";
    
    cout<<"Adding Skimmed EMu Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_03feb2017rereco_em_v2_skim*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_em_v1_skim*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_em_v1_skim*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_em_v1_skim*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_em_v1_skim*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_em_v1_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_em_v2_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_em_v3_skim*");
  }
  else if (set == "DileptonData-mm-Skimmed"){
    //TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-11/skim/";
    TString dir="/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-15/skims/";

    cout<<"Adding Skimmed MuMu Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_03feb2017rereco_mm_v2_skim*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_mm_v1_skim*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_mm_v1_skim*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_mm_v1_skim*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_mm_v1_skim*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_mm_v1_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_mm_v2_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_mm_v3_skim*");
  }
  else if (set == "SingleLeptonData-SingleMu-Skimmed"){
    //TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-11/skim/";
    TString dir="/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-15/skims/";
    
    cout<<"Adding Skimmed Single Muon Trigger Data"<<endl;
    
    ch->Add(dir+"data_Run2016B_03feb2017rereco_sm_v2_skim*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_sm_v1_skim*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_sm_v1_skim*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_sm_v1_skim*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_sm_v1_skim*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_sm_v1_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_sm_v2_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_sm_v3_skim*");
  }

//====================================
// Skimmed Photon Data
//====================================

  //Single Photon Trigger
  else if (set == "GammaData-SinglePhoton-Skimmed"){
    //TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-11/skim/";
    TString dir="/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-15/skims/";

    cout<<"Adding GammaData-SinglePhoton-Skimmed"<<endl; 
    
    ch->Add(dir+"data_Run2016B_03feb2017rereco_ph_v2_skim*");
    ch->Add(dir+"data_Run2016C_03feb2017rereco_ph_v1_skim*");
    ch->Add(dir+"data_Run2016D_03feb2017rereco_ph_v1_skim*");
    ch->Add(dir+"data_Run2016E_03feb2017rereco_ph_v1_skim*");
    ch->Add(dir+"data_Run2016F_03feb2017rereco_ph_v1_skim*");
    ch->Add(dir+"data_Run2016G_03feb2017rereco_ph_v1_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_ph_v2_skim*");
    ch->Add(dir+"data_Run2016H_03feb2017rereco_ph_v3_skim*");
  }
  //EWK Subtraction
  else if (set == "GammaData-EWKSub-Skimmed"){
    //TString dir="/hadoop/cms/store/user/olivito/AutoTwopler_babies/merged/ZMET/V08-22-16/skim/";
    TString dir="/nfs-7/userdata/ZMEToutput/output/ZMETbabies/V08-22-16/skims/";
    
    cout<<"Adding Skimmed EWK Subtraction Samples"<<endl;       
    
    //============
    // W+Gamma+Jets
    //============
    //This is the Wjets sample, it is intended to have events with a prompt photon vetod
    ch->Add(dir+"wjets_incl_mgmlm*");
    ch->Add(dir+"wjets_ht100_mgmlm*");
    ch->Add(dir+"wjets_ht200_mgmlm*");
    ch->Add(dir+"wjets_ht400_mgmlm*");
    ch->Add(dir+"wjets_ht600_mgmlm*");
    ch->Add(dir+"wjets_ht800_mgmlm*");
    ch->Add(dir+"wjets_ht1200_mgmlm*");
    ch->Add(dir+"wjets_ht2500_mgmlm*");

    //This is the W+Gamma+Jets, it is inteded to have events with non-prompt photons vetod
    ch->Add(dir+"wgjets_incl_mgmlm*");
    ch->Add(dir+"wgjets_ptg40_mgmlm*");
    ch->Add(dir+"wgjets_ptg130_mgmlm*");
    //============
    // Gamma+Z->NuNu
    //============
    ch->Add(dir+"znunugamma_ptg40_mgmlm*");
    ch->Add(dir+"znunugamma_ptg130_mgmlm*");
    //============
    // TTbar
    //============
    //1lep
    ch->Add(dir+"ttbar_1ltbr_mgmlm_ext1*");
    ch->Add(dir+"ttbar_1ltop_mgmlm_ext1*");
    //dilep
    ch->Add(dir+"ttbar_dilep_mgmlm_ext1*");
    //============
    // Single Top
    //============
    ch->Add(dir+"sttw_antitop_nofullhaddecay_powheg*");
    ch->Add(dir+"sttw_top_nofullhaddecay_powheg*");
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