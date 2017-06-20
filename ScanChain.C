#include "ScanChain.h"

//=============================
// Variable Computation
//=============================

pair<int, int> getMostBlike(){
  /* returns two most B-like jet indicies */
  int first_index = 0;
  int second_index = 1;

  if (g_jets_csv.at(first_index) < g_jets_csv.at(second_index)){
    first_index = 1;
    second_index = 0;
  }

  for (int i = 2; i < (int) g_jets_p4.size(); i++){
    if (g_jets_csv.at(first_index) < g_jets_csv.at(i)){
      second_index = first_index;
      first_index = i;
    }
    else if (g_jets_csv.at(second_index) < g_jets_csv.at(i)){
      second_index = i;
    }
  }

  return make_pair(first_index, second_index);
}

pair<int,int> getClosestBPairToHiggsMass(){
  int first = 0;
  int second = 1;
  
  if (g_jets_p4.size()<2){
    cout<<"Going to throw error finding closest B pair: Less than two jets!"<<endl;
    std::stringstream message;
    message<<"Can not find closest B pair to Higgs for event: "<<phys.evt()<<" run: "<<phys.run()<<" lumi: "<<phys.lumi()<<" the event only has "<<g_jets_p4.size()<<" jets";
    throw std::underflow_error(message.str());
  }

  double dist = fabs((g_jets_p4.at(0)+g_jets_p4.at(1)).M() - 125.0);

  for (int i = 0; i < (int) g_jets_p4.size(); i++) {
    for (int j = i+1; j < (int) g_jets_p4.size(); j++) {
      if (fabs((g_jets_p4.at(i)+g_jets_p4.at(j)).M() - 125 ) < dist){
        first = i;
        second = j;
      }
    }
  }

  return make_pair(first,second);
}

double getMT2B(){
  /*Builds MT2b from two highest CSV jets*/

  pair<int, int> b_index = getMostBlike();

  double mt2_1=MT2(g_met, g_met_phi, g_jets_p4.at(b_index.first)+phys.lep_p4().at(g_lep_inds.at(0)), g_jets_p4.at(b_index.second)+phys.lep_p4().at(g_lep_inds.at(1)), 0, 0);
  double mt2_2=MT2(g_met, g_met_phi, g_jets_p4.at(b_index.first)+phys.lep_p4().at(g_lep_inds.at(1)), g_jets_p4.at(b_index.second)+phys.lep_p4().at(g_lep_inds.at(0)), 0, 0);

  if (mt2_1 > mt2_2){
    return mt2_2;
  }
  else{
    return mt2_1;
  }
}

double getMbb(){
  /*Builds Mbb from two highest CSV jets*/
  
  pair<int, int> b_index = getMostBlike();

  return (g_jets_p4.at(b_index.first)+g_jets_p4.at(b_index.second)).M();
}

double getMT2_JL(){
  /*Builds MT2 (jet lepton sum) from two most W-like jets*/

  pair<int, int> indicies = getMostWlikePair(g_jets_p4);

  double mt2_1=MT2(g_met, g_met_phi, g_jets_p4.at(indicies.first)+phys.lep_p4().at(g_lep_inds.at(0)), g_jets_p4.at(indicies.second)+phys.lep_p4().at(g_lep_inds.at(1)), 0, 0);
  double mt2_2=MT2(g_met, g_met_phi, g_jets_p4.at(indicies.first)+phys.lep_p4().at(g_lep_inds.at(1)), g_jets_p4.at(indicies.second)+phys.lep_p4().at(g_lep_inds.at(0)), 0, 0);

  if (mt2_1 > mt2_2){
    return mt2_2;
  }
  else{
    return mt2_1;
  }
}

double getMT2ForBjets(bool select_highest_csv/*=false*/){
  /*This function gets the MT2 built out of the two Bjets in an event, no guarentee is made about selecting the highest csv jets*/
  double mt2;
  //cout<<__LINE__<<endl;
  if (select_highest_csv){
    pair<int, int> b_index = getMostBlike();
    //cout<<__LINE__<<endl;
    //make sure first index points to the higher csv of the first two jets
    mt2=MT2(g_met, g_met_phi, g_jets_p4.at(b_index.first), g_jets_p4.at(b_index.second), 0, 0);
    //cout<<__LINE__<<endl;
  }
  else{
    // MT2( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2, MASS_INVISIBLE_PARTICLE, Bool Verbose)
    //cout<<__LINE__<<endl;
    mt2=MT2(g_met, g_met_phi, g_jets_medb_p4.at(0), g_jets_medb_p4.at(1), 0, 0);
  }
  //cout<<__LINE__<<endl;
  return mt2;
}

double getMT2HiggsZ(bool select_highest_closest_higgs_mass/*=false*/){

  double mt2; 

  if(select_highest_closest_higgs_mass){
    pair<int, int> jet_indexes = getClosestBPairToHiggsMass();
    // MT2( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2, MASS_INVISIBLE_PARTICLE, Bool Verbose)
    mt2=MT2(g_met, g_met_phi, g_jets_p4.at(jet_indexes.first)+g_jets_p4.at(jet_indexes.second), phys.lep_p4().at(g_lep_inds.at(0))+phys.lep_p4().at(g_lep_inds.at(1)), 0, 0);
  }
  else{
    // MT2( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2, MASS_INVISIBLE_PARTICLE, Bool Verbose)
    mt2=MT2(g_met, g_met_phi, g_jets_medb_p4.at(0)+g_jets_medb_p4.at(1), phys.lep_p4().at(g_lep_inds.at(0))+phys.lep_p4().at(g_lep_inds.at(1)), 0, 0);
  }

  return mt2;
}

double bosonPt(){
  // Returns boson Pt, determines whether sample is gjets or zjets first.
  if (conf->get("event_type") == "dilepton") {
    return (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).pt();
  }
  else{
    if (phys.evt_type() == 2 && phys.ngamma() > 0){
      return phys.gamma_pt().at(0);
    }
    else
      return 0;
  }
}

double getMTLepMET(short id/*=0*/){
  /* Builds the MT from the lepton at index id and the MET vector (assumes massless particles)*/
  return sqrt(g_met*phys.lep_p4().at(id).pt()*(1 - cos(g_met_phi - phys.lep_p4().at(id).phi())));

  /* Massive Case
    ET1sq = m_1^2 + pt1^2 
    ET2sq = m_2^2 + pt2^2
    MT^2 = m_1^2  + m_2^2 + 2(sqrt(ET1sq*ET2sq) - pt1*pt2*cos(phi1 - phi2)) 
  */
}

double getdRGammaLep(short id/*=0*/){
  /* Builds the delta R (sqrt(dPhi^2 + dEta^2)) between the lepton at index id and the leading photon*/
  double dPhi=acos(cos(phys.gamma_p4().at(0).phi() - phys.lep_p4().at(id).phi()));
  double dEta=phys.gamma_p4().at(0).eta() - phys.lep_p4().at(id).eta();

  return sqrt(pow(dPhi, 2) + pow(dEta, 2));
}

pair<int, int> getPairWithMass(const vector<LorentzVector> &vecs, double target_mass, bool close){
  /*Searches the vector of lorentz vectors for the pair with mass nearest (furthest) from the target mass if 'close' is true (false)*/
  
  short mult = (close) ? 1 : -1;
  if (vecs.size() < 2){
    cout<<"Throwing error, vector does not have at least 2 elements in getPairWithMass()"<<endl;
    std::stringstream message;
    message<<"Can not find pair for mass "<<target_mass<<" in evt: "<<phys.evt()<<" run: "<<phys.run()<<" lumi: "<<phys.lumi()<<". The vector given has "<<vecs.size()<<" entries, must have at least 2.";
    throw std::underflow_error(message.str());
  }

  int first_index = 0;
  int second_index = 1;

  double best_mass = mult*fabs(target_mass - (vecs.at(first_index)+vecs.at(second_index)).M());

  for (int i = 0; i < (int) vecs.size(); i++){
    for (int j = i+1; j < (int) vecs.size(); j++){
      if (mult*fabs(target_mass - (vecs.at(i)+vecs.at(j)).M()) < best_mass){
        first_index=i;
        second_index=j;
        best_mass = mult*fabs(target_mass - (vecs.at(first_index)+vecs.at(second_index)).M());
      }
    }
  }

  return make_pair(first_index, second_index);
}

pair<int,int> getMostZlikePair(const vector<LorentzVector> &vecs){
  /*Finds the pair with 2 vector mass closest to the mass of the Z (91 GeV) for the vector of LorentzVector objects given*/
  return getPairWithMass(vecs, 91, true);
}

pair<int,int> getLeastZlikePair(const vector<LorentzVector> &vecs){
  /*Finds the pair with 2 vector mass furthest from the mass of the Z (91 GeV) for the vector of LorentzVector objects given*/
  return getPairWithMass(vecs, 91, false);
}

pair<int,int> getMostWlikePair(const vector<LorentzVector> &vecs){
  /*Finds the pair with 2 vector mass closest to the mass of the W (80 GeV) for the vector of LorentzVector objects given*/
  return getPairWithMass(vecs, 80, true);
}

pair<int,int> getClosestJetsInEta(){
  /* Goes pairwise through jet collection and finds the pair of jets with the minium delta eta. */
  if (g_njets < 2){
    cout<<"Going to throw error finding closest dEta jet pair: Less than two jets!"<<endl;
    std::stringstream message;
    message<<"Can not find closest dEta jet pair to Higgs for event: "<<phys.evt()<<" run: "<<phys.run()<<" lumi: "<<phys.lumi()<<" the event only has "<<g_jets_p4.size()<<" jets";
    throw std::underflow_error(message.str());
  }

  int j1 = 0;
  int j2 = 1;
  double best_deta = fabs(g_jets_p4.at(j1).eta() - g_jets_p4.at(j2).eta());

  for (int i = 0; i < (int) g_jets_p4.size()-1; i++){
    for (int j = i+1; j < (int) g_jets_p4.size(); j++){
      if (fabs(g_jets_p4.at(i).eta() - g_jets_p4.at(j).eta()) < best_deta){
        j1 = i;
        j2 = j;
      }
    }
  }

  return make_pair(j1,j2);
}

int getNumSFOSPairs(){
  /*Loops through pairs of entries in the lep_pdgId vector and counts how many have opposite value*/
  int num_SFOS = 0;
  for (int i = 0; i < (int) g_lep_inds.size(); i++){
    for (int j = i+1; j < (int) g_lep_inds.size(); j++){
      //cout<<"looking at pair with ID ("<<phys.lep_pdgId().at(i)<<", "<<phys.lep_pdgId().at(j)<<")"<<endl;
      if (phys.lep_pdgId().at(g_lep_inds.at(i)) == -phys.lep_pdgId().at(g_lep_inds.at(j))) num_SFOS++;
    }
  }
  //cout<<"num SFOS pairs: "<<num_SFOS<<endl;
  return num_SFOS;
}

TString getLepFlavorString(){
  /*Returns a string which gives the exact flavor composition of the leptons in the event*/
  int num_e=0;
  int num_ebar=0;
  int num_m=0;
  int num_mbar=0;

  for (int i = 0; i<(int) g_lep_inds.size(); i++){
    //cout<<"lep"<<i<<": "<<phys.lep_pdgId().at(i)<<endl;
    if (phys.lep_pdgId().at(g_lep_inds.at(i)) == 11) num_e++;
    else if (phys.lep_pdgId().at(g_lep_inds.at(i)) == -11) num_ebar++;
    else if (phys.lep_pdgId().at(g_lep_inds.at(i)) == 13) num_m++;
    else if (phys.lep_pdgId().at(g_lep_inds.at(i)) == -13) num_mbar++;
  }

  TString flavor = "";
  for (int j = 0; j<num_e; j++){
    flavor += "E";
  }
  for (int j = 0; j<num_ebar; j++){
    flavor += "Ebar";
  }
  for (int j = 0; j<num_m; j++){
    flavor += "Mu";
  }
  for (int j = 0; j<num_mbar; j++){
    flavor += "Mubar";
  }
  //cout<<flavor<<endl;

  return flavor;
}

double DeltaR(const LorentzVector p1, const LorentzVector p2){
  /*Returns the DeltaR between objects p1 and p2.*/
  //cout<<__LINE__<<endl;
  return sqrt( (p1.eta() - p2.eta())*(p1.eta() - p2.eta())+(p1.phi() - p2.phi())*(p1.phi() - p2.phi()) );
}

bool isCleanLepFromW(int index){
  /* checks the gen record for a lepton with the index specified */
  LorentzVector lp4 = phys.lep_p4().at(index);
  int lpdgId = phys.lep_pdgId().at(index);

  for (int i = 0; i < (int) phys.genPart_p4().size(); i++){
    if (phys.genPart_pdgId().at(i) == lpdgId){
      //cout<<"Found lepton "<<index<<" in evt: "<<phys.evt()<<" Gen Record with (pt, eta, phi) = ("<<phys.genPart_p4().at(i).pt()<<", "<<phys.genPart_p4().at(i).eta()<<", "<<phys.genPart_p4().at(i).phi()<<"). Lep has ("<<phys.lep_p4().at(index).pt()<<", "<<phys.lep_p4().at(index).eta()<<", "<<phys.lep_p4().at(index).phi()<<")."<<endl;
      if (DeltaR(phys.genPart_p4().at(i), lp4) < 0.2){
        if (fabs(phys.genPart_p4().at(i).pt() - lp4.pt()) < 5){
          if (fabs(phys.genPart_motherId().at(i)) == 24){
            //cout<<"Was a match!"<<endl;
            //Found Lepton with:
            // DR within 0.2
            // Pt within 5 GeV
            // W mother
            return true;
          }
        }
      }
    }
  }

  //cout<<"No Match Found for evt: "<<phys.evt()<<endl;

  return false;
}

flavor_type getFlavorType(){
  /* Checks whether the two leading analysis leptons are EE, EMu, or MuMu */
  if ( (abs(phys.lep_pdgId().at(g_lep_inds.at(0))) == 11) && (abs(phys.lep_pdgId().at(g_lep_inds.at(1))) == 13) ) return EMu;
  //cout<<__LINE__<<endl;
  if ( (abs(phys.lep_pdgId().at(g_lep_inds.at(0))) == 13) && (abs(phys.lep_pdgId().at(g_lep_inds.at(1))) == 11) ) return EMu;
  //cout<<__LINE__<<endl;
  if ( (abs(phys.lep_pdgId().at(g_lep_inds.at(0))) == 11) && (abs(phys.lep_pdgId().at(g_lep_inds.at(1))) == 11) ) return EE;
  //cout<<__LINE__<<endl;
  if ( (abs(phys.lep_pdgId().at(g_lep_inds.at(0))) == 13) && (abs(phys.lep_pdgId().at(g_lep_inds.at(1))) == 13) ) return MuMu;
  //cout<<__LINE__<<endl;

  cout<<"Going to throw error finding flavor type in event"<<endl;
  std::stringstream message;
  message<<"Can not find flavor type for event: "<<phys.evt()<<" run: "<<phys.run()<<" lumi: "<<phys.lumi()<<" the event only has objects with pdgIds "<<phys.lep_pdgId().at(g_lep_inds.at(0))<<" and "<<phys.lep_pdgId().at(g_lep_inds.at(1))<<" jets";
  throw std::invalid_argument(message.str());
}

charge_type getChargeType(){
  /* Returns the charge type, SS or OS for 2 lep events, and SFOS0/1/2 for more than 2leps */
  if (g_lep_inds.size() == 2){
    int id1 = phys.lep_pdgId().at(g_lep_inds.at(0));
    int id2 = phys.lep_pdgId().at(g_lep_inds.at(1));
    if ( (id1/abs(id1)) == (id2/abs(id2)) )   return SS;
    else                                      return OS;
  }
  
  //cout<<__LINE__<<endl;

  int nsfof = getNumSFOSPairs();

  if (nsfof == 0)      return SFOS0;
  else if (nsfof == 1) return SFOS1;
  else if (nsfof == 2) return SFOS2;
  
  //cout<<__LINE__<<endl;

  return NA; //In case we have more than 3 leps which can happen without lep veto. 
}

//=============================
// Triggers
//=============================

bool passPhotonEmulatedTrigger() {
  if( phys.gamma_r9()            .at(0) < 0.92 ) return false;
  if( phys.gamma_hOverE()        .at(0) > 0.2  ) return false;
  if( phys.gamma_hollowtkiso03() .at(0) > 5    ) return false;
  if( fabs(phys.gamma_eta().at(0)) < 1.4 && phys.gamma_ecpfclusiso().at(0) > 4 + phys.gamma_pt().at(0) * 0.0053  ) return false;
  if( fabs(phys.gamma_eta().at(0)) < 1.4 && phys.gamma_hcpfclusiso().at(0) > 8 + phys.gamma_pt().at(0) * 0.014   ) return false;
  if( fabs(phys.gamma_eta().at(0)) > 1.6 && phys.gamma_ecpfclusiso().at(0) > 4 + phys.gamma_pt().at(0) * 0.0034  ) return false;
  if( fabs(phys.gamma_eta().at(0)) > 1.6 && phys.gamma_hcpfclusiso().at(0) > 8 + phys.gamma_pt().at(0) * 0.0139  ) return false;

  return true;
}

bool passPhotonTriggers(){
  if ( (! MCTriggerEmulation) && (! phys.isData()) ){
    return true;
  }
  else{
    if( ((phys.HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon165_R9Id90_HE10_IsoM() > 0) || (phys.HLT_Photon165_HE10_matchedtophoton() && phys.HLT_Photon165_HE10() > 0)) && phys.gamma_pt().at(0) > 180 ) return true;
    else if( phys.HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon120_R9Id90_HE10_IsoM() > 0 && phys.gamma_pt().at(0) > 135 && phys.gamma_pt().at(0) < 180 ) return true;
    else if( phys.HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon90_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 105 && phys.gamma_pt().at(0) < 135  ) return true;
    else if( phys.HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon75_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 85 && phys.gamma_pt().at(0) < 105   ) return true;
    else if( phys.HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon50_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 55 && phys.gamma_pt().at(0) < 85    ) return true;
    else if( phys.HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon36_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 40 && phys.gamma_pt().at(0) < 55    ) return true;
    else if( phys.HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon30_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 33 && phys.gamma_pt().at(0) < 40    ) return true;
    else if( phys.HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon22_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) < 33 ) return true;
    return false;
  }
}

bool passMuonTriggers(){
  if ( (! MCTriggerEmulation) && (! phys.isData()) ){
    return true;
  }
  else{
    //cout<<__LINE__<<endl;
    if ( conf->get("use_muon_DZ_triggers") == "true" ){
      //cout<<"Using DZ triggers"<<endl;
      //cout<<__LINE__<<endl;
      //if (printStats) { cout<<"HLT_DoubleMu: "<<phys.HLT_DoubleMu()<<" HLT_DoubleMu_tk: "<<phys.HLT_DoubleMu_tk()<<" "<<" HLT_DoubleMu_noiso: "<<phys.HLT_DoubleMu_noiso()<<" "; }
      return (phys.HLT_DoubleMu() || phys.HLT_DoubleMu_tk() || phys.HLT_DoubleMu_noiso());
    }
    else{
      //cout<<"Using Non DZ triggers"<<endl;
      //cout<<__LINE__<<endl;
      //if (printStats) { cout<<"HLT_DoubleMu_nonDZ: "<<phys.HLT_DoubleMu_nonDZ()<<" HLT_DoubleMu_tk_nonDZ: "<<phys.HLT_DoubleMu_tk_nonDZ()<<" "<<" HLT_DoubleMu_noiso: "<<phys.HLT_DoubleMu_noiso()<<" "; }
      if(conf->get("signal_region") == "LeonoraEvtLists"){
        return (phys.HLT_DoubleMu() || phys.HLT_DoubleMu_tk() /*|| phys.HLT_DoubleMu_dbltk()*/ || phys.HLT_DoubleMu_nonDZ() || phys.HLT_DoubleMu_tk_nonDZ() || phys.HLT_DoubleMu_noiso());
      }
      else{
        return (phys.HLT_DoubleMu() || phys.HLT_DoubleMu_tk() /*|| phys.HLT_DoubleMu_dbltk()*/ || phys.HLT_DoubleMu_nonDZ() || phys.HLT_DoubleMu_tk_nonDZ() || phys.HLT_DoubleMu_noiso());
      }
    } 
  }
}

bool passElectronTriggers(){
  if ( (! MCTriggerEmulation) && (! phys.isData()) ){
    return true;
  }
  else{
    //cout<<__LINE__<<endl;
    //if (printStats) { cout<<"HLT_DoubleEl_DZ_2: "<<phys.HLT_DoubleEl_DZ_2()<<" HLT_DoubleEl_noiso: "<<phys.HLT_DoubleEl_noiso()<<" HLT_DoubleEl_DZ(): "<<phys.HLT_DoubleEl_DZ()<<endl; }
    return (phys.HLT_DoubleEl_DZ_2() || phys.HLT_DoubleEl_noiso() || phys.HLT_DoubleEl_DZ() );
  }
}

bool passEMuTriggers(){
  if ( (! MCTriggerEmulation) && (! phys.isData()) ){
    return true;
  }
  else{
    return (phys.HLT_MuEG() || phys.HLT_MuEG_2() || phys.HLT_MuEG_noiso() || phys.HLT_MuEG_noiso_2() || phys.HLT_Mu8_EG23_DZ() || phys.HLT_Mu23_EG12_DZ() || phys.HLT_Mu23_EG8_DZ()); //updated for Moriond 2017
  }
}

bool passSingleMuTriggers(){
  if ( (! MCTriggerEmulation) && (! phys.isData()) ){
    return true;
  }
  else{
    return (phys.HLT_singleMu());
  }
}

bool passLeptonHLTs(){
  flavor_type ft = getFlavorType();

  if ( ft == MuMu ){ //Muon Event
    return passMuonTriggers();
  }
  else if ( ft == EE ){ //Electron Event
    return passElectronTriggers();
  }
  else if ( ft == EMu) { //Emu event
    return passEMuTriggers(); 
  }
  else{
    cout<<"Throwing error, got bad flavor type? This shouldn't be possible"<<endl;
    std::stringstream message;
    message<<"Flavor type for event: "<<phys.evt()<<" run: "<<phys.run()<<" lumi: "<<phys.lumi()<<" is not possible... Inside passLeptonHLTs()";
    throw std::underflow_error(message.str());
    return false;
  }
}

//=============================
// Has Good Event Functions
//=============================

bool hasGood3l(){
  //if (printStats) { cout<<"Number of Leptons: "<<phys.nlep()<<" "; }
  
  if ( ((int) g_lep_inds.size()) != stoi(conf->get("num_leptons"))){
      numEvents->Fill(10); 
      if (printFail) cout<<phys.evt()<<" :Failed num leptons cut. Needed "<<conf->get("num_leptons")<<" got "<<g_lep_inds.size()<<endl;
      return false; // require 2 leps
  }

  int lep_charge_sum = abs(phys.lep_pdgId().at(g_lep_inds.at(0))/abs(phys.lep_pdgId().at(g_lep_inds.at(0))) + phys.lep_pdgId().at(g_lep_inds.at(1))/abs(phys.lep_pdgId().at(g_lep_inds.at(1))) + phys.lep_pdgId().at(g_lep_inds.at(2))/abs(phys.lep_pdgId().at(g_lep_inds.at(2))));

  if( lep_charge_sum != 1 ){ 
    numEvents->Fill(70);
    if (printFail) cout<<phys.evt()<<" :Failed total lepton sum not 1"<<endl;
    return false; // no WWW production should have total charge greater than 1
  }

  //cout<<__LINE__<<endl;

  if ( conf->get("3lep_phi_MET_min") != "" ){
    if ( fabs(acos(cos((phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).phi() - g_met_phi))) <= stod(conf->get("3lep_phi_MET_min")) ){
      numEvents->Fill(31); 
      if (printFail) cout<<phys.evt()<<" :Failed phi of 3lep system within "<<stod(conf->get("3lep_phi_MET_min"))<<" of MET phi"<<endl;
      return false;
    }
  }

  if (conf->get("pt_3l_min") != ""){
    double pt_3l = ( phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2)) ).pt();
    if (pt_3l <= stod(conf->get("pt_3l_min"))){
      numEvents->Fill(27);
      if (printFail) cout<<phys.evt()<<" :Failed 3 lep system pt cut of "<<conf->get("pt_3l_min")<<" with value "<<pt_3l<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if( conf->get("num_SFOS_pairs") != "" ) {
    if ( stoi(conf->get("num_SFOS_pairs")) != getNumSFOSPairs()){
      numEvents->Fill(21); 
      if (printFail) cout<<phys.evt()<<" :Failed Num SFOS pairs cut"<<endl;
      return false; // require correct number of SFOS pairs
    }
  }

  //cout<<__LINE__<<endl;

  //if (printStats) { cout<<"evt_type: "<<phys.evt_type()<<" "; }*/

  if ( conf->get("Mll_SF_min") != "" ){
    if ( fabs(phys.lep_pdgId().at(g_lep_inds.at(0))) == fabs(phys.lep_pdgId().at(g_lep_inds.at(1))) ){
      if ( (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M() <= stod(conf->get("Mll_SF_min")) ){
        numEvents->Fill(24); 
        if (printFail) cout<<phys.evt()<<" :Failed Mll SF cut in 3 lep region for leps 1 and 2"<<endl;
        return false;
      }
    }

    //cout<<__LINE__<<endl;

    if ( fabs(phys.lep_pdgId().at(g_lep_inds.at(0))) == fabs(phys.lep_pdgId().at(g_lep_inds.at(2))) ){
      if ( (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(2))).M() <= stod(conf->get("Mll_SF_min")) ){
        numEvents->Fill(24); 
        if (printFail) cout<<phys.evt()<<" :Failed Mll SF cut in 3 lep region for leps 1 and 3"<<endl;
        return false;
      }
    }

    //cout<<__LINE__<<endl;

    if ( fabs(phys.lep_pdgId().at(g_lep_inds.at(1))) == fabs(phys.lep_pdgId().at(g_lep_inds.at(2))) ){
      if ( (phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).M() <= stod(conf->get("Mll_SF_min")) ){
        numEvents->Fill(24); 
        if (printFail) cout<<phys.evt()<<" :Failed Mll SF cut in 3 lep region for leps 2 and 3"<<endl;
        return false;
      }
    }

    //cout<<__LINE__<<endl;
  }

  if( ( conf->get("dilmass_Z_veto") == "true" ) ) {
    if ( conf->get("signal_region") == "3lep_0SFOS" ){
      //Cut has been removed
      /*if ( ( fabs(phys.lep_pdgId().at(g_lep_inds.at(0))) == fabs(phys.lep_pdgId().at(g_lep_inds.at(1))) ) && ( fabs(phys.lep_pdgId().at(g_lep_inds.at(0))) == 11) ){
        if ( fabs((phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M() - Z_MASS ) <= 15 ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 0SFOS pairs"<<endl;
          return false;
        }
      }

      //cout<<__LINE__<<endl;

      if ( ( fabs(phys.lep_pdgId().at(g_lep_inds.at(0))) == fabs(phys.lep_pdgId().at(g_lep_inds.at(2))) ) && ( fabs(phys.lep_pdgId().at(g_lep_inds.at(0))) == 11) ){
        if ( fabs((phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(2))).M() - Z_MASS ) <= 15 ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 0SFOS pairs"<<endl;
          return false;
        }
      }

      //cout<<__LINE__<<endl;

      if ( ( fabs(phys.lep_pdgId().at(g_lep_inds.at(1))) == fabs(phys.lep_pdgId().at(g_lep_inds.at(2))) ) && ( fabs(phys.lep_pdgId().at(g_lep_inds.at(1))) == 11) ){
        if ( fabs((phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).M() - Z_MASS ) <= 15 ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 0SFOS pairs"<<endl;
          return false;
        }
      }*/
    }

    else if ( conf->get("signal_region") == "3lep_1SFOS" ){
      if ( ( phys.lep_pdgId().at(g_lep_inds.at(0)) == -phys.lep_pdgId().at(g_lep_inds.at(1)) ) ){
        //veto if mass within [MZ-35,MZ+20]
        if ( ( (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M() >= 55 ) && ( (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M() <= 110 ) ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 1SFOS pairs"<<endl;
          return false;
        }
      }

      //cout<<__LINE__<<endl;

      if ( ( phys.lep_pdgId().at(g_lep_inds.at(0)) == -phys.lep_pdgId().at(g_lep_inds.at(2)) ) ){
        //veto if mass within [MZ-35,MZ+20]
        if ( ( (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(2))).M() >= 55 ) && ( (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(2))).M() <= 110 ) ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 1SFOS pairs"<<endl;
          return false;
        }
      }

      //cout<<__LINE__<<endl;

      if ( ( phys.lep_pdgId().at(g_lep_inds.at(1)) == -phys.lep_pdgId().at(g_lep_inds.at(2)) ) ){
        //veto if mass within [MZ-35,MZ+20]
        if ( ( (phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).M() >= 55 ) && ( (phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).M() <= 110 ) ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 1SFOS pairs"<<endl;
          return false;
        }
      }

      //cout<<__LINE__<<endl;
    }

    else if ( conf->get("signal_region") == "3lep_2SFOS" ){
      if ( ( phys.lep_pdgId().at(g_lep_inds.at(0)) == -phys.lep_pdgId().at(g_lep_inds.at(1)) ) ){
        //veto if mass within [MZ-35,MZ+20]
        if ( fabs((phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M() - Z_MASS ) <= 20 ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 2SFOS pairs"<<endl;
          return false;
        }
      }

      //cout<<__LINE__<<endl;

      if ( ( phys.lep_pdgId().at(g_lep_inds.at(0)) == -phys.lep_pdgId().at(g_lep_inds.at(2)) ) ){
        //veto if mass within [MZ-35,MZ+20]
        if ( fabs((phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(2))).M() - Z_MASS ) <= 20 ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 2SFOS pairs"<<endl;
          return false;
        }
      }

      //cout<<__LINE__<<endl;

      if ( ( phys.lep_pdgId().at(g_lep_inds.at(1)) == -phys.lep_pdgId().at(g_lep_inds.at(2)) ) ){
        //veto if mass within [MZ-35,MZ+20]
        if ( fabs((phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).M() - Z_MASS ) <= 20 ){
          numEvents->Fill(19); 
          if (printFail) cout<<phys.evt()<<" :Failed dilmass Z veto for 3 lepton region with 2SFOS pairs"<<endl;
          return false;
        }
      }

      //cout<<__LINE__<<endl;
    }
    else{
      cout<<"Throwing error, can not apply dilmass Z veto without a known signal region"<<endl;
      std::stringstream message;
      message<<"Don't know how to apply the dilepton mass Z veto with signal region: "<<conf->get("signal_region");
      throw std::invalid_argument(message.str());
    }
  }

  //cout<<__LINE__<<endl;

  if( ( conf->get("dilmass_ee_veto") == "true" ) ) {
    //check each pair of leps to see if they are electron type, then check for the mass window cut

    if ( (abs(phys.lep_pdgId().at(g_lep_inds.at(0))) == 11) && (abs(phys.lep_pdgId().at(g_lep_inds.at(1))) == 11) ){
      double dilmass = ( phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1)) ).M();
      if (fabs(dilmass - Z_MASS) <= 15 ){
        numEvents->Fill(25);
        if (printFail) cout<<phys.evt()<<" :Failed Mll ee cut in 3 lep region for leps 1 and 2 with dilmass"<<dilmass<<endl;
        return false;
      }
    }

    if ( (abs(phys.lep_pdgId().at(g_lep_inds.at(0))) == 11) && (abs(phys.lep_pdgId().at(g_lep_inds.at(2))) == 11) ){
      double dilmass = ( phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(2)) ).M();
      if (fabs(dilmass - Z_MASS) <= 15 ){
        numEvents->Fill(25);
        if (printFail) cout<<phys.evt()<<" :Failed Mll ee cut in 3 lep region for leps 1 and 3 with dilmass"<<dilmass<<endl;
        return false;
      }
    }

    if ( (abs(phys.lep_pdgId().at(g_lep_inds.at(1))) == 11) && (abs(phys.lep_pdgId().at(g_lep_inds.at(2))) == 11) ){
      double dilmass = ( phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2)) ).M();
      if (fabs(dilmass - Z_MASS) <= 15 ){
        numEvents->Fill(25);
        if (printFail) cout<<phys.evt()<<" :Failed Mll ee cut in 3 lep region for leps 2 and 3 with dilmass"<<dilmass<<endl;
        return false;
      }
    }
  }

  //cout<<__LINE__<<endl;

  return true;
}

bool hasGood2l(){
  //cout<<__LINE__<<endl;

  if ( ((int) g_lep_inds.size()) != stoi(conf->get("num_leptons"))){
      numEvents->Fill(10); 
      if (printFail) cout<<phys.evt()<<" :Failed num leptons cut. Needed "<<conf->get("num_leptons")<<" got "<<g_lep_inds.size()<<endl;
      return false; // require 2 leps
  }

  if (conf->get("dil_sign") == "same"){
    if( !( CT == SS ) ) {
      numEvents->Fill(21); 
      if (printFail) cout<<phys.evt()<<" :Failed same sign cut"<<endl;
      return false; // require same sign
    }
  }
  else if (conf->get("dil_sign") == "opposite"){
    if( !( CT == OS ) ) {
      numEvents->Fill(21); 
      if (printFail) cout<<phys.evt()<<" :Failed opposite sign cut"<<endl;
      return false; // require opposite sign
    }
  }

  if (conf->get("dil_flavor") == "emu"){ 
    if (! ( FT == EMu ) ){ //require explicit emu event
      numEvents->Fill(15); 
      if (printFail) cout<<phys.evt()<<" :Failed not explicit e/mu event"<<endl;
      return false; // require explicit opposite flavor event
    }
  }
  else if (conf->get("dil_flavor") == "same"){
    if( !( FT == EE || FT == MuMu ) ) {
        numEvents->Fill(15); 
        if (printFail) cout<<phys.evt()<<" :Failed explicit mu/mu or e/e cut"<<endl;
        return false; // require explicit same flavor event
    }
  }
  else if (conf->get("dil_flavor") == "ee"){
    if( !( FT == EE ) ) {
        numEvents->Fill(15); 
        if (printFail) cout<<phys.evt()<<" :Failed explicit e/e cut"<<endl;
        return false; // require explicit same flavor event
    }
  }
  else if (conf->get("dil_flavor") == "mumu"){
    if( !( FT == MuMu ) ) {
        numEvents->Fill(15); 
        if (printFail) cout<<phys.evt()<<" :Failed explicit mu/mu cut"<<endl;
        return false; // require explicit same flavor event
    }
  }
  //cout<<"pdgIDs: "<<phys.lep_pdgId().at(g_lep_inds.at(0))<<" "<<phys.lep_pdgId().at(g_lep_inds.at(1))<<endl;

  //cout<<__LINE__<<endl;

  if ((conf->get("dilmass_Z_veto") == "true") && (FT == EE) ){ //only apply for EE events
    //cout<<g_lep_inds.at(0)<<" "<<g_lep_inds.at(1)<<" "<<g_lep_inds.size()<<" "<<phys.nlep()<<endl;
    if ( ( (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M() > Z_VETO_WINDOW_LOW ) && ( (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M() < Z_VETO_WINDOW_HIGH ) ) {
      numEvents->Fill(22); 
      if (printFail) cout<<phys.evt()<<" :Failed lepton pair on Z mass cut with mass: "<<(phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M()<<endl;
      return false; // require opposite sign
    }
  }

  //cout<<__LINE__<<endl;

  if (g_njets >= 2){
    // MJJ for closest jets in eta
    pair<int,int> jets_dEta = getClosestJetsInEta();
    double dijet_mass = ( g_jets_p4.at(jets_dEta.first) + g_jets_p4.at(jets_dEta.second) ).M();
    if( (dijet_mass < W_JET_WINDOW_LOW) || (dijet_mass > W_JET_WINDOW_HIGH) ) {
      numEvents->Fill(61); 
      if (printFail) cout<<phys.evt()<<" :Failed jet pair near W mass cut with dijet mass: "<< dijet_mass <<endl;
      return false; // require at at least 2 Jets within the W mass window.
    }
    //Safety cuts against Wpm Wpm VBS
    double leading_dijet_mass = ( g_jets_p4.at(0) + g_jets_p4.at(1) ).M();
    if ( leading_dijet_mass >= 400 ){
      numEvents->Fill(58); 
      if (printFail) cout<<phys.evt()<<" :Failed leading dijet mass >= 400 GeV: "<< leading_dijet_mass <<endl;
      return false; // require at at least 2 Jets within the W mass window.
    }
  }
  else{
    numEvents->Fill(34); 
    if (printFail) cout<<phys.evt()<<" :Failed Not enough jets"<<endl;
    return false; // require at least 2 jets
  }

  //cout<<__LINE__<<endl;
  
  //if (printPass) cout<<phys.evt()<<": Passes good Z Cuts"<<endl;
  return true;
}

bool hasGoodEvent() {
  if (conf->get("num_leptons") == "3"){
    return hasGood3l();
  }
  else{
    return hasGood2l();
  }
}

//=============================
// Event Weight Assignment
//=============================

void readyReweightHists(){
  TString conf_name = conf->get("Name");

  //cout<<"FINDFIND Adding "<<conf->get("Name");

  cout<<"Reweighting with "<<TString(conf->get("histo_output_dir")+"ct_"+conf->get("rwt_var")+"_"+conf->get("signal_region")+"_rwt.root")<<endl;
  TString rwt_hist_name = "h_"+conf->get("rwt_var")+"_ratio";
  TFile *reweight_file = TFile::Open( TString(conf->get("histo_output_dir")+"ct_"+conf->get("rwt_var")+"_"+conf->get("signal_region")+"_rwt.root"), "READ");
  g_reweight_pairs.push_back(make_pair( (TH1D*) reweight_file->Get(rwt_hist_name)->Clone(TString("reweight_hist_")+conf->get("rwt_var")),conf->get("rwt_var")));
  g_reweight_pairs.back().first->SetDirectory(rootdir);
  reweight_file->Close();

  while (conf->get("weight_from") != "" ){
    conf->loadConfig(conf->get("weight_from"));
    //cout<<"FINDFIND Adding "<<conf->get("Name");
    cout<<"Reweighting with "<<TString(conf->get("histo_output_dir")+"ct_"+conf->get("rwt_var")+"_"+conf->get("signal_region")+"_rwt.root")<<endl;
    rwt_hist_name = "h_"+conf->get("rwt_var")+"_ratio";
    reweight_file = TFile::Open( TString(conf->get("histo_output_dir")+"ct_"+conf->get("rwt_var")+"_"+conf->get("signal_region")+"_rwt.root"), "READ");
    g_reweight_pairs.push_back(make_pair( (TH1D*) reweight_file->Get(rwt_hist_name)->Clone(TString("reweight_hist_")+conf->get("rwt_var")),conf->get("rwt_var")));
    g_reweight_pairs.back().first->SetDirectory(rootdir);
    reweight_file->Close();      
  }

  conf->loadConfig(conf_name.Data());
  cout<<"Reweight hists loaded, proceeding with conf "<<conf->get("Name")<<endl;
}

void readyVPTReweight(TString save_path){
  /* Adds the vpt reweighting histogram to the g_reweight_pairs vector */

  TString vpt_weight_path = save_path+conf->get("Name")+"_vpt_rwt.root";
  TString rwt_hist_name = "vpt_ratio";

  cout<<"Reweighting with "<<vpt_weight_path<<endl;

  TFile *reweight_file = TFile::Open(vpt_weight_path, "READ");
  
  //Add pair (vpt_weight, "vpt") to g_reweight_pairs
  g_reweight_pairs.push_back(make_pair( (TH1D*) reweight_file->Get(rwt_hist_name)->Clone(TString("vpt_reweight_hist")),"vpt"));
  g_reweight_pairs.back().first->SetDirectory(rootdir);
  
  reweight_file->Close();
}

double getEff(const double &pt, const double &eta){
  /* Returns the trigger efficiency from g_pt_eff */
  if (fabs(eta) < 1.4){
    return g_pt_eff_barrel->GetEfficiency(g_pt_eff_barrel->FindFixBin(pt));
  }
  else{
    return g_pt_eff_endcap->GetEfficiency(g_pt_eff_endcap->FindFixBin(pt)); 
  }
}

double getReweight(){
  double weight = 1;
  
  TH1D* rwt_hist;
  TString rwt_var;
  //cout<<"Size: "<<g_reweight_pairs.size()<<endl;
  for (int i=0; i<(int)g_reweight_pairs.size(); i++){
    rwt_hist = g_reweight_pairs.at(i).first;
    rwt_var = g_reweight_pairs.at(i).second;
    //cout<<rwt_var<<endl;

    if (rwt_var == "vpt"){
      //cout<<"Applying vpt reweight -- pt: "<<bosonPt()<<" weight: "<<rwt_hist->GetBinContent(rwt_hist->FindBin(bosonPt()))<<endl;
      weight *= rwt_hist->GetBinContent(rwt_hist->FindBin(bosonPt()));
    }
    else if (rwt_var == "ht_wide"){
      //cout<<"Adding HT weight: "<<rwt_hist->GetBinContent(rwt_hist->FindBin(ht))<<endl;
      weight *= rwt_hist->GetBinContent(rwt_hist->FindBin(g_ht)); 
    }
    else if (rwt_var == "muonPT_endcap"){
      //cout<<"Adding muon pt trigger efficiency weight in endcap: "<<rwt_hist->GetBinContent(rwt_hist->FindBin(phys.lep_pt().at(0)))<<endl;
      if (phys.lep_p4().at(g_lep_inds.at(0)).eta() > 1.6){
        weight *= rwt_hist->GetBinContent(rwt_hist->FindBin(phys.lep_pt().at(g_lep_inds.at(0)))); 
      }
    }
    else if (rwt_var == "muonPT_barrel"){
      //cout<<"Adding muon pt trigger efficiency weight in barrel: "<<rwt_hist->GetBinContent(rwt_hist->FindBin(phys.lep_pt().at(0)))<<endl;
      if (phys.lep_p4().at(g_lep_inds.at(0)).eta() < 1.4){
        weight *= rwt_hist->GetBinContent(rwt_hist->FindBin(phys.lep_pt().at(g_lep_inds.at(0)))); 
      }
    }
    else{
      std::stringstream message;
      message<<"Reweight varible is not a valid option, please choose vpt, or ht_wide, got: "<<rwt_var<<".";
      throw std::invalid_argument(message.str());
    }
  }

  return weight;
}

double scale1fbFix(){
  /*This method stores fixes to the evt_scale1fb in the event of file corruptions. It's basically just a lookup table*/

  if (TString(currentFile->GetTitle()).Contains("sttw_antitop_nofullhaddecay_powheg")){
    //cout<<"Scale 1fb fixed for "<<TString(currentFile->GetTitle())<<endl;
    return 1.03;
  }
  else{
    //cout<<"Scale 1fb is good for "<<TString(currentFile->GetTitle())<<endl;
    return 1;
  }
}

double getWeight(){
  /*Gets the proper weight for the sample. */
  double weight=1;
  //cout<<__LINE__<<endl;
  if ( ! ( phys.isData() ) ){
    weight *= phys.evt_scale1fb();
    
    //Weight to some other lumi
    if ( conf->get("scaleTofb") != "" ){
      weight *= stod(conf->get("scaleTofb"));
    }
    
    //cout<<__LINE__<<endl;

    if (conf->get("pileup_reweight") == "true"){
      weight*=g_pileup_hist->GetBinContent(g_pileup_hist->FindBin(phys.nTrueInt()));
    }
    
    //cout<<__LINE__<<endl;
    
    if (TString(conf->get("data_set")).Contains("GammaData-EWKSub")){
      weight *= -1; //EWK Subtraction
    }
  }
  //cout<<__LINE__<<endl;

  /*weight *= g_scale_factor;

  if ( conf->get("reweight") == "true" || conf->get("vpt_reweight") == "true") {
    weight *= getReweight();
  }

  //cout<<__LINE__<<endl;

  if (conf->get("rwt_photon_eff") == "true" ){
    weight *= getEff(phys.gamma_pt().at(0), phys.gamma_eta().at(0));
  }

  //cout<<__LINE__<<endl;

  if (conf->get("rwt_muon_eff") == "true"){
    weight *= getEff(phys.lep_pt().at(g_lep_inds.at(0)), phys.lep_eta().at(g_lep_inds.at(0)));
  }

  //cout<<__LINE__<<endl;

 if ((! phys.isData()) ){    
    if ( conf->get("event_type") == "dilepton" && (! MCTriggerEmulation)){
      //simulated trigger efficiencies
      if ( FT == EE   ) weight *= 0.969; //ee
      if ( FT == MuMu ) weight *= 0.980; //mumu
      if ( FT == EMu  ) weight *= 0.932; //emu
    }

    //cout<<__LINE__<<endl;

    for (int i = 0; i < phys.nlep(); i++){
      //cout<<__LINE__<<endl;
      weight*=phys.weightsf_lepid().at(g_lep_inds.at(i));
      //cout<<__LINE__<<endl;
      weight*=phys.weightsf_lepiso().at(g_lep_inds.at(i));
      //cout<<__LINE__<<endl;
      weight*=phys.weightsf_lepip().at(g_lep_inds.at(i));
      //cout<<__LINE__<<endl;
      weight*=phys.weightsf_lepreco().at(g_lep_inds.at(i));
      //cout<<__LINE__<<endl;
      weight*=phys.weightsf_lepconv().at(g_lep_inds.at(i));
      //cout<<__LINE__<<endl;
    }

    //cout<<__LINE__<<endl;

    if (conf->get("no_btag_sf") != "true"){
      //cout<<"Applying Btag Scale Factors"<<endl;
      weight *= phys.weight_btagsf();
    }

    if (conf->get("susy_mc") == "true"){
      double ISR_norm, btag_norm;
      if(conf->get("SUSY_Glu_LSP_scan") == "true"){
        ISR_norm=1./g_isr_norm->GetBinContent(g_isr_norm->GetXaxis()->FindBin(phys.mass_gluino()), g_isr_norm->GetYaxis()->FindBin(phys.mass_LSP()));
        //cout<<__LINE__<<endl;
        btag_norm=1./g_btagsf_norm->GetBinContent(g_btagsf_norm->GetXaxis()->FindBin(phys.mass_gluino()), g_btagsf_norm->GetYaxis()->FindBin(phys.mass_LSP()));
      }
      else if (conf->get("SUSY_chi_scan") == "true"){
        //JUST READ FROM BIN 1 IN THE LSP ROW
        ISR_norm=1./g_isr_norm->GetBinContent(g_isr_norm->GetXaxis()->FindBin(phys.mass_chi()), 1);
        //cout<<__LINE__<<endl;
        btag_norm=1./g_btagsf_norm->GetBinContent(g_btagsf_norm->GetXaxis()->FindBin(phys.mass_chi()), 1);
      }
      else{
        std::stringstream message;
        message<<"Can not get ISR or Btag normalization if SUSY_chi_scan or SUSY_Glu_LSP_scan are not set.";
        throw std::invalid_argument(message.str());
      }
      //cout<<__LINE__<<endl;
      weight *= phys.isr_weight(); //ISR scale factor
      weight *= ISR_norm;
      
      //cout<<__LINE__<<endl;

      weight *= btag_norm;

      //cout<<__LINE__<<endl;
      
      for (int i = 0; i < phys.nlep(); i++){
        weight *= phys.weightsf_lepid_FS().at(i); //Fast Sim Lepton ID
        weight *= phys.weightsf_lepiso_FS().at(i); //Fast Sim Lepton isolation
        weight *= phys.weightsf_lepip_FS().at(i); //Fast Sim Lepton impact parameter
      } 
      //cout<<__LINE__<<endl;
    }
  }
  //cout<<__LINE__<<endl;

  if (phys.isData() && phys.ngamma() > 0 && TString(currentFile->GetTitle()).Contains("data") && TString(currentFile->GetTitle()).Contains("_ph")){
    weight *= getPrescaleWeight();
  }


  //cout<<__LINE__<<endl;*/

  /*if (weight < 0 || weight > 0.3){
    cout<<"Odd Weight: "<<weight<<" "<<phys.evt()<<endl;
  }*/

  //weight *= scale1fbFix();
  //cout<<"weight: "<<weight<<" evt: "<<phys.evt()<<endl;

  return weight;
}

double getPrescaleWeight(){
  //cout<<__LINE__<<endl;
  //cout<<"Getting Prescale Weights"<<endl;
  if( ((phys.HLT_Photon165_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon165_R9Id90_HE10_IsoM() > 0) || (phys.HLT_Photon165_HE10_matchedtophoton() && phys.HLT_Photon165_HE10() > 0)) && phys.gamma_pt().at(0) > 180 ) return 1;
  else if( phys.HLT_Photon120_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon120_R9Id90_HE10_IsoM() > 0 && phys.gamma_pt().at(0) > 135 && phys.gamma_pt().at(0) < 180 ) return phys.HLT_Photon120_R9Id90_HE10_IsoM();
  else if( phys.HLT_Photon90_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon90_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 105 && phys.gamma_pt().at(0) < 135  ) return phys.HLT_Photon90_R9Id90_HE10_IsoM();
  else if( phys.HLT_Photon75_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon75_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 85 && phys.gamma_pt().at(0) < 105   ) return phys.HLT_Photon75_R9Id90_HE10_IsoM();
  else if( phys.HLT_Photon50_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon50_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 55 && phys.gamma_pt().at(0) < 85    ) return phys.HLT_Photon50_R9Id90_HE10_IsoM();
  else if( phys.HLT_Photon36_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon36_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 40 && phys.gamma_pt().at(0) < 55    ) return /*166*/ 168.309;
  else if( phys.HLT_Photon30_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon30_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 33 && phys.gamma_pt().at(0) < 40    ) return /*354*/ 379.058;
  else if( phys.HLT_Photon22_R9Id90_HE10_IsoM_matchedtophoton() && phys.HLT_Photon22_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) < 33 ) return /*1871*/ 1985.08;
  return 0;

  //if( (phys.HLT_Photon165_R9Id90_HE10_IsoM() > 0 || phys.HLT_Photon165_HE10() > 0) && phys.gamma_pt().at(0) > 180. ) return 1;
  //else if( phys.HLT_Photon120_R9Id90_HE10_IsoM() > 0 && phys.gamma_pt().at(0) > 135. ) return phys.HLT_Photon120_R9Id90_HE10_IsoM();
  //else if( phys.HLT_Photon90_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 105. ) return phys.HLT_Photon90_R9Id90_HE10_IsoM();
  //else if( phys.HLT_Photon75_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 85. ) return phys.HLT_Photon75_R9Id90_HE10_IsoM();
  //else if( phys.HLT_Photon50_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) > 55. ) return phys.HLT_Photon50_R9Id90_HE10_IsoM();
  //else if( phys.HLT_Photon36_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) < 55. && phys.gamma_pt().at(0) > 40. ) {
  //  return /*g_l1prescale_hist36->GetBinContent(g_l1prescale_hist36->FindBin(phys.nVert())) */ 166 /*134*/;
  //}
  //else if( phys.HLT_Photon30_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) < 40. && phys.gamma_pt().at(0) > 33. ){
  //  return /*g_l1prescale_hist30->GetBinContent(g_l1prescale_hist30->FindBin(phys.nVert())) */ 354 /*269*/;
  //}
  //else if( phys.HLT_Photon22_R9Id90_HE10_IsoM()  > 0 && phys.gamma_pt().at(0) < 33. ) {
  //  return /*g_l1prescale_hist22->GetBinContent(g_l1prescale_hist22->FindBin(phys.nVert())) */ 1871 /*1667*/;
  //}
  //cout<<__LINE__<<endl;
  //return 1; // should not get here
}

//=============================
// Cuts
//=============================

bool passSignalRegionCuts(){
  
  //Njets Min Cut
  if (conf->get("Njets_min") != ""){
    if (g_njets < stod(conf->get("Njets_min"))){
      numEvents->Fill(34);
      if (printFail) cout<<phys.evt()<<" :Failed min jets cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  //Njets Max Cut
  if (conf->get("Njets_max") != ""){
    if (g_njets > stod(conf->get("Njets_max"))){
      numEvents->Fill(35);
      if (printFail) cout<<phys.evt()<<" :Failed max jets cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;
  //if (printStats) { cout<<"NbjetsMed: "<<g_nBJetMedium<<" "; }

  //Num Bottom jets Min Cut
  if (conf->get("NBjets_min") != ""){
    if (g_nBJetMedium < stod(conf->get("NBjets_min"))){
      numEvents->Fill(36);
      if (printFail) cout<<phys.evt()<<" :Failed min bjet cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  //Num Bottom jets Max Cut
  if (conf->get("NBjets_max") != ""){
    if (g_nBJetMedium > stod(conf->get("NBjets_max"))){
      numEvents->Fill(37);
      if (printFail) cout<<phys.evt()<<" :Failed max bjet cut"<<endl;
      return false;
    }
  }

  //Num Bottom jets Max Cut
  if (conf->get("NBjets_loose_max") != ""){
    //cout<<"BJET evt: "<<phys.evt()<<" -- Checking for max loose bjets, have: "<<g_nBJetLoose<<endl;
    if (g_nBJetLoose > stod(conf->get("NBjets_loose_max"))){
      numEvents->Fill(37);
      if (printFail) cout<<phys.evt()<<" :Failed max bjet cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  //Num Bottom jets Min Cut
  if (conf->get("NBjets_loose_min") != ""){
    if (g_nBJetLoose < stod(conf->get("NBjets_loose_min"))){
      numEvents->Fill(36);
      if (printFail) cout<<phys.evt()<<" :Failed min bjet cut"<<endl;
      return false;
    }
  }


  //cout<<__LINE__<<endl;
  //if (printStats) { cout<<"g_dphi_metj1: "<<g_dphi_metj1<<" "; }
  //Leading Jet/MET Phi min
  if (conf->get("dPhi_MET_j1") != ""){
    if (g_dphi_metj1 < stod(conf->get("dPhi_MET_j1"))){
      numEvents->Fill(38);
      if (printFail) cout<<phys.evt()<<" :Failed dPhi MET with jet 1 cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;
  //if (printStats) { cout<<"g_dphi_metj2: "<<g_dphi_metj2<<" "; }
  //Trailing Jet/MET Phi min
  if (conf->get("dPhi_MET_j2") != ""){
    if (g_dphi_metj2 < stod(conf->get("dPhi_MET_j2"))){
      numEvents->Fill(39);
      if (printFail) cout<<phys.evt()<<" :Failed dPhi MET with jet 2 cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;
  //if (printStats) { cout<<"mt2b: "<<g_mt2b<<" "; }
  //MT2b min
  if (conf->get("MT2b_min") != "" && conf->get("event_type") != "photon"){
    if (g_mt2b < stod(conf->get("MT2b_min"))){
      numEvents->Fill(40);
      if (printFail) cout<<phys.evt()<<" :Failed MT2b cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  //MT2b min
  if (conf->get("MT2b_loose_min") != "" && conf->get("event_type") != "photon"){
    if (getMT2B() < stod(conf->get("MT2b_loose_min"))){
      numEvents->Fill(40);
      if (printFail) cout<<phys.evt()<<" :Failed MT2b cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;
  //if (printStats) { cout<<"mt2b: "<<g_mt2b<<" "; }
  //MT2 min
  if (conf->get("MT2_min") != ""){
    if (g_mt2 < stod(conf->get("MT2_min"))){
      numEvents->Fill(59);
      if (printFail) cout<<phys.evt()<<" :Failed MT2 cut"<<endl;
      return false;
    }
  }


  //cout<<__LINE__<<endl;
  //HT min
  if (conf->get("HT_min") != ""){
  //if (printStats) { cout<<"ht: "<<g_ht<<" "; }
    if (g_ht < stod(conf->get("HT_min"))){
      numEvents->Fill(41);
      if (printFail) cout<<phys.evt()<<" :Failed sum HT min cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;
  //DiBottom mass difference from Higgs Mass
  if (conf->get("mbb_mh_diff") != ""){
  //if (printStats) { cout<<"mbb_mh_diff: "<<fabs(mbb - 125)<<" "; }
    if (fabs(g_mbb - 125) < stod(conf->get("mbb_mh_diff"))){
      numEvents->Fill(42);
      if (printFail) cout<<phys.evt()<<" :Failed sum mbb_mh diff cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  //Wierd ATLAS SR cut
  if (conf->get("sum_HT_pt_pt") != ""){
    double pt;
    //cout<<__LINE__<<endl;
    if (phys.evt_type() == 2 && phys.ngamma() > 0){
      pt = phys.gamma_pt().at(0);
    }
    else{
      pt = phys.lep_pt().at(g_lep_inds.at(0)) + phys.lep_pt().at(g_lep_inds.at(1));
    }
    //cout<<__LINE__<<endl;
    //if (printStats) { cout<<"sum_HT_pt_pt: "<<fabs(ht + pt )<<" "; }
    if ( fabs(g_ht + pt ) < stod(conf->get("sum_HT_pt_pt") ) ){
      numEvents->Fill(43);
      if (printFail) cout<<phys.evt()<<" :Failed sum HT pt pt cut"<<endl;
      return false;
    }
  }
  
  //cout<<__LINE__<<endl;

  if (conf->get("lep1_pt_min") != "" && conf->get("event_type") != "photon" ){
    if ( phys.lep_pt().at(g_lep_inds.at(0)) < stod( conf->get("lep1_pt_min") )){
      numEvents->Fill(45);
      if (printFail) cout<<phys.evt()<<" :Failed lep1 min pt"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("lep2_pt_min") != "" && conf->get("event_type") != "photon" ){
    if ( phys.lep_pt().at(g_lep_inds.at(1)) < stod( conf->get("lep2_pt_min") )){
      numEvents->Fill(46);
      if (printFail) cout<<phys.evt()<<" :Failed lep2 min pt cut"<<endl;
      return false;
    }
  }

  if (conf->get("lep3_pt_min") != "" && conf->get("event_type") != "photon" ){
    if ( phys.lep_pt().at(g_lep_inds.at(2)) < stod( conf->get("lep3_pt_min") )){
      numEvents->Fill(46);
      if (printFail) cout<<phys.evt()<<" :Failed lep3 min pt cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("MET_min") != ""){
    if ( g_met < stod( conf->get("MET_min") )){
      numEvents->Fill(56);
      if (printFail) cout<<phys.evt()<<" :Failed MET min cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("MET_max") != ""){
    if ( g_met > stod( conf->get("MET_max") )){
      numEvents->Fill(57);
      if (printFail) cout<<phys.evt()<<" :Failed MET max cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("Mbb_max") != ""){
    if ( g_mbb > stod( conf->get("Mbb_max") )){
      numEvents->Fill(58);
      if (printFail) cout<<phys.evt()<<" :Failed Mbb max cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("Mbb_loose_max") != ""){
    if ( getMbb() > stod( conf->get("Mbb_loose_max") )){
      numEvents->Fill(58);
      if (printFail) cout<<phys.evt()<<" :Failed Mbb loose max cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("MT_LepMET_min") != ""){
    if ( getMTLepMET() < stod( conf->get("MT_LepMET_min") ) ){
      numEvents->Fill(63);
      if (printFail) cout<<phys.evt()<<" :Failed MT from Lepton and MET min cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("Mjj_dphi_max") != ""){
    if ( g_mjj_mindphi > stod( conf->get("Mjj_dphi_max") ) ){
      numEvents->Fill(68);
      if (printFail) cout<<phys.evt()<<" :Failed Mjj cut"<<endl;
      return false;
    }
  }

  /*if (conf->get("Mjj_Z_veto") == "true"){
    pair<int,int> jetind = getMostZlikePair(g_jets_p4);
    if( fabs((g_jets_p4.at(jetind.first) + g_jets_p4.at(jetind.second)).M() - Z_MASS) < Z_VETO_WINDOW ) {
      numEvents->Fill(9); 
      if (printFail) cout<<phys.evt()<<" :Failed jet pair on Z mass cut"<<endl;
      return false; // Mjj Z Veto
    }
  }*/

  //cout<<__LINE__<<endl;
  

  //cout<<__LINE__<<endl;

  if (conf->get("flavor") == "dimuon"){
    if ( FT != MuMu ){
      numEvents->Fill(73);
      if (printFail) cout<<phys.evt()<<" :Failed dimuon cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("flavor") == "dielectron"){
    if ( FT != EE ){
      numEvents->Fill(73);
      if (printFail) cout<<phys.evt()<<" :Failed dielectron cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("force_true_bjets") != ""){
    pair<int, int> b_index = getMostBlike();

    if ( fabs(phys.jets_mcFlavour().at(b_index.first)) != 5 ){
      numEvents->Fill(67);
      if (printFail) cout<<phys.evt()<<" :Failed truth level bjet cut"<<endl;
      return false;
    }

    if ( fabs(phys.jets_mcFlavour().at(b_index.second)) != 5 ){
      numEvents->Fill(67);
      if (printFail) cout<<phys.evt()<<" :Failed truth level bjet cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("force_fake_bjets") != ""){
    pair<int, int> b_index = getMostBlike();

    if ( (fabs(phys.jets_mcFlavour().at(b_index.first)) == 5) && (fabs(phys.jets_mcFlavour().at(b_index.second)) == 5)){
      numEvents->Fill(67);
      if (printFail) cout<<phys.evt()<<" :Failed truth level bjet cut"<<endl;
      return false;
    }
  }

  //cout<<__LINE__<<endl;

  if (conf->get("3chargeagree") == "true"){
    for (int i = 0; i < stoi(conf->get("num_leptons")); i++){
      if ( ! phys.lep_3ch_agree().at(g_lep_inds.at(i)) ){
        numEvents->Fill(65);
        if (printFail) cout<<phys.evt()<<" :Failed Three Charge Agree cut at lepton "<<i<<endl;
        return false;
      }
    }
  }

  if (conf->get("dEta_jj_max") != ""){
      if ( fabs(g_jets_p4.at(0).eta() - g_jets_p4.at(1).eta()) >  stod(conf->get("dEta_jj_max"))){
        numEvents->Fill(74);
        if (printFail) cout<<phys.evt()<<" :Failed Delta Eta JJ Max cut at "<<stod(conf->get("dEta_jj_max"))<<endl;
        return false;
      }
  }

  if (conf->get("Mll_min") != ""){
      if ( (phys.lep_p4().at(g_lep_inds.at(0))+phys.lep_p4().at(g_lep_inds.at(1))).M() <  stod(conf->get("Mll_min"))){
        numEvents->Fill(48);
        if (printFail) cout<<phys.evt()<<" :Failed Mll min cut at "<<stod(conf->get("Mll_min"))<<endl;
        return false;
      }
  }

  if (conf->get("reliso04_max") != ""){
    for (int i = 0; i < stoi(conf->get("num_leptons")); i++){
      if ( phys.lep_relIso04EA().at(g_lep_inds.at(i)) > stod(conf->get("reliso04_max")) ){
        numEvents->Fill(66);
        if (printFail) cout<<phys.evt()<<" :Failed reliso04 max: "<<stod(conf->get("reliso04_max"))<<" at lepton "<<i<<endl;
        return false;
      }
    }
  }

  if (conf->get("reliso03_max") != ""){
    for (int i = 0; i < stoi(conf->get("num_leptons")); i++){
      if ( phys.lep_relIso03EA().at(g_lep_inds.at(i)) > stod(conf->get("reliso03_max")) ){
        numEvents->Fill(66);
        if (printFail) cout<<phys.evt()<<" :Failed reliso03 max: "<<stod(conf->get("reliso03_max"))<<" at lepton "<<i<<endl;
        return false;
      }
    }
  }

  if (conf->get("ip3d_max") != ""){
    for (int i = 0; i < stoi(conf->get("num_leptons")); i++){
      if ( fabs(phys.lep_ip3d().at(g_lep_inds.at(i))) > stod(conf->get("ip3d_max")) ){
        numEvents->Fill(66);
        if (printFail) cout<<phys.evt()<<" :Failed ip3d_max max: "<<stod(conf->get("ip3d_max"))<<" at lepton "<<i<<endl;
        return false;
      }
    }
  }

  //cout<<__LINE__<<endl;

  if ( (conf->get("jet1_pt_min") != "" ) && (g_njets > 0) ){
    if ( g_jets_p4.at(0).pt() < stod( conf->get("jet1_pt_min") )){
      numEvents->Fill(49);
      if (printFail) cout<<phys.evt()<<" :Failed jet1 min pt"<<endl;
      return false;
    }
  }

  if (conf->get("use_veto_leps") == "true"){
    if (phys.nlep_VVV_cutbased_veto() != stoi(conf->get("num_leptons"))){ //so long as we use cutbased tight ID this is okay, otherwise we need more work here to ensure we veto 'extra' leps that pass this selection (i.e. there could be some leps that pass MVA based but not cut based, so would need to do overlap checking.)
      numEvents->Fill(54);
      if (printFail) cout<<phys.evt()<<" :Failed veto lepton selections"<<endl;
      return false;
    }
  }


  if (conf->get("isotrack_veto") == "true" ){
    if (conf->get("use_veto_leps") == "true"){
      if ( phys.nisoTrack_mt2_cleaned_VVV_cutbased_veto() != 0){
        numEvents->Fill(62);
        if (printFail) cout<<phys.evt()<<" :Failed isotrack veto"<<endl;
        return false;
      }
    }
    else{
      if ( phys.nisoTrack_mt2_cleaned_VVV_cutbased_tight() != 0){
        numEvents->Fill(62);
        if (printFail) cout<<phys.evt()<<" :Failed isotrack veto"<<endl;
        return false;
      }
    }
  }

  //cout<<__LINE__<<endl;
  //if (printPass) cout<<phys.evt()<<": Passes Signal Region Cuts"<<endl;
  return true;
}

bool isDuplicate(){
  //cout<<__LINE__<<endl;
  if( phys.isData() ) {
    DorkyEventIdentifier id(phys.run(), phys.evt(), phys.lumi());
    //cout<<__LINE__<<endl;
    if (is_duplicate(id) ){
      ++nDuplicates;
      //cout<<__LINE__<<endl;
      if (printFail) cout<<phys.evt()<<" :Is a duplicate"<<endl;
      return true;
    }
    //cout<<__LINE__<<endl;
  }
  //cout<<__LINE__<<endl;
  //if (printPass) cout<<phys.evt()<<": Passes not a duplicate"<<endl;
  return false;
}

bool passMETFilters(){
  //updated for Moriond 2017
  if (!phys.Flag_EcalDeadCellTriggerPrimitiveFilter()      ) { 
    numEvents->Fill(5);
    if (printFail) cout<<phys.evt()<<" :Failed EcalDeadCellTriggerPrimativeFilter cut"<<endl;
    return false;
  }
  if (!phys.Flag_HBHENoiseFilter                   ()      ){ 
    numEvents->Fill(2);
    if (printFail) cout<<phys.evt()<<" :Failed HBHENoiseFilter cut"<<endl;
    return false;
  }
  if (!phys.Flag_HBHEIsoNoiseFilter                ()      ){ 
    numEvents->Fill(3);
    if (printFail) cout<<phys.evt()<<" :Failed HBHEIsoNoiseFilter cut"<<endl;
    return false;
  } 
  if (!phys.Flag_goodVertices                      ()      ) { 
    numEvents->Fill(6);
      if (printFail) cout<<phys.evt()<<" :Failed goodVerticies cut"<<endl;
    return false;
  }
  if(conf->get("signal_region") != "LeonoraEvtLists"){
    if (phys.nJet200MuFrac50DphiMet() > 0){
      numEvents->Fill(70);
      if (printFail) cout<<phys.evt()<<" :Failed nJet200MuFrac50DphiMet cut"<<endl;
      return false;
    }
    if ((g_met / phys.met_calo_pt()) > 5){
      numEvents->Fill(71);
      if (printFail) cout<<phys.evt()<<" :Failed T1MET/CaloMET cut"<<endl;
      return false;
    }
  }
  if (conf->get("susy_mc") != "true"){
    if (!phys.Flag_globalTightHalo2016            ()      ){ 
      numEvents->Fill(4);
      if (printFail) cout<<phys.evt()<<" :Failed CSCTightHalo2015Filter cut"<<endl;
      return false;
    }
  }
  if ( phys.isData() ) {
    if (!phys.Flag_eeBadScFilter                     ()      ) { 
      numEvents->Fill(7);
      if (printFail) cout<<phys.evt()<<" :Failed eeBadScFilter cut"<<endl;
      return false;
    }
  }
  
  //if (printPass) cout<<phys.evt()<<": Passes MET Filters"<<endl;
  return true;
}

bool passBaseCut(){
  //if (printStats) { cout<<"goodrun : "<<goodrun(phys.evt(), phys.lumi())<<" "; }
  //if (printStats) { cout<<"njets : "<<g_njets<<" "; }
  
  //bool pass=true;

  if (phys.isData()){
    if (! (goodrun(phys.run(), phys.lumi()))){ 
      numEvents->Fill(8);
      if (printFail) cout<<phys.evt()<<" :Failed golden JSON cut"<<endl;
      return false; //golden json
    }
  }

  if(conf->get("num_leptons") != ""){
    if (conf->get("fakerate_study") == "true"){ 
      //When doing FR study, there are two ways to pass, either 2 tight leps, or 1 tight lep and 2 fo (which implies veto) leps.
      if ((phys.nlep_VVV_cutbased_tight() != stod(conf->get("num_leptons")) ) && ( !(phys.nlep_VVV_cutbased_tight() == (stod(conf->get("num_leptons"))-1) && phys.nlep_VVV_cutbased_fo() == stod(conf->get("num_leptons"))))){ 
        numEvents->Fill(10);
        if (printFail) cout<<phys.evt()<<" :Failed 2 lepton cut"<<endl;
        return false; // require at least 2 good leptons
      } 
    }
    //if (printStats) { cout<<"Number of Leptons: "<<phys.nlep()<<" "; }
    }
    else{
      if( phys.nlep_VVV_cutbased_tight() != stod(conf->get("num_leptons")) ){ 
        numEvents->Fill(10);
        if (printFail) cout<<phys.evt()<<" :Failed 2 lepton cut"<<endl;
        return false; // require at least 2 good leptons
      }
  }

  if(conf->get("tau_veto") == "true"){
    if (phys.nTaus20() >= 1){
      numEvents->Fill(75);
      if (printFail) cout<<phys.evt()<<" :Failed Tau veto"<<endl;
      return false;
    }
  }

  //if (printPass) cout<<phys.evt()<<": Passes Base Cuts"<<endl;
  //return pass;
  return true;
}

bool passFileSelections(){
  /* Method which holds all the file specific selections, for instance cutting out the
  events with genht > 100 in the DY inclusive samples */


  //Zjets Monte Carlo samples
  if ( (! phys.isData()) && TString(conf->get("data_set")).Contains("DY")){
    //cout<<"Zjets MC event"<<endl;
    if( TString(currentFile->GetTitle()).Contains("dy_m50_mgmlm") && (! TString(currentFile->GetTitle()).Contains("_ht")) ){
      //cout<<"File: "<<currentFile->GetTitle()<<" with gen_ht: "<<phys.gen_ht()<<endl;
      if( phys.gen_ht() > 100 ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(44);
        return false;
      }
      if(phys.evt_scale1fb() > 30){
        numEvents->Fill(60);
        return false;
      }
      if( fabs(phys.gen_ht() - g_ht) > 300 ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(69);
        return false;
      }
    }
    if( TString(currentFile->GetTitle()).Contains("dy_m50_mgmlm_ht100")){
      if( fabs(phys.gen_ht() - g_ht) > 300 ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(69);
        return false;
      }
    }
  }

  //Wjets Monte Carlo samples
  if ( (! phys.isData()) && TString(conf->get("data_set")).Contains("WJets")){
    //cout<<"Wjets MC event"<<endl;
    if( TString(currentFile->GetTitle()).Contains("wjets_incl_mgmlm_") && (! TString(currentFile->GetTitle()).Contains("_ht")) ){
      //cout<<"File: "<<currentFile->GetTitle()<<" with gen_ht: "<<phys.gen_ht()<<endl;
      if( phys.gen_ht() > 100 ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(44);
        return false;
      }
    }
  }

  //Photon MC samples
  if ( (! phys.isData()) && TString(conf->get("data_set")).Contains("GammaMC")){
    //cout<<"Photon MC event"<<endl;
    if( TString(currentFile->GetTitle()).Contains("gjetsht40") ||  TString(currentFile->GetTitle()).Contains("gjetsht100") ){
      if( fabs(phys.gen_ht() - g_ht) > 300 ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(69);
        return false;
      }
    }
  }

  //WJets cocktail for inclusive photon sample and Electroweak Subtraction
  if ( TString(conf->get("data_set")).Contains("GammaMC-WGamma") || TString(conf->get("data_set")).Contains("GammaData-EWKSub")){
    
    //Inclusive GenHT Cut
    if( TString(currentFile->GetTitle()).Contains("wjets_incl_mgmlm") ){
      //cout<<"File: "<<currentFile->GetTitle()<<" with gen_ht: "<<phys.gen_ht()<<endl;
      if( phys.gen_ht() > 100 ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(44);
        return false;
      }
    }

    //Remove overlap between WGammaJets and WJets
    if( TString(currentFile->GetTitle()).Contains("wjets") ){ //WJets
      if( phys.ngamma() > 0 && phys.gamma_genIsPromptFinalState().at(0) == 1 ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(64);
        return false;
      }
    }
    else if ( TString(currentFile->GetTitle()).Contains("wgjets_") ){ //WGammaJets
      if( phys.ngamma() > 0 && phys.gamma_genIsPromptFinalState().at(0) != 1 ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(64);
        return false;
      }
    }

    //remove the events with less than 40 gen pt
    if ( TString(currentFile->GetTitle()).Contains("wgjets_incl_mgmlm") && phys.ngamma() > 0){
      //First get best photon match
      int bestMatch = -1;
      float bestDR = 0.1;
      double eta = phys.gamma_eta().at(0);
      double phi = phys.gamma_phi().at(0);
      for(unsigned int iGen = 0; iGen < phys.genPart_pdgId().size(); iGen++){
        if ((phys.genPart_pdgId().at(iGen) != 22) /*&& (fabs(phys.genPart_pdgId().at(iGen)) != 11)*/ ) continue; // accept gen photons
        if ((fabs(phys.genPart_motherId().at(iGen)) > 24) && (fabs (phys.genPart_motherId().at(iGen)) != 2212) ) continue; // don't want stuff from pions etc 
        if (phys.genPart_status().at(iGen) != 1  ) continue; 
        if (fabs(eta - phys.genPart_eta().at(iGen)) > 0.1 ) continue;
        float thisDR = sqrt(pow(eta-phys.genPart_eta().at(iGen), 2) + pow(phi-phys.genPart_phi().at(iGen), 2)); //DeltaR( cms3.genps_p4() .at(iGen).eta(), eta, cms3.genps_p4().at(iGen).phi(), phi);
        if (thisDR < bestDR) {
          bestMatch = iGen;
          bestDR=thisDR;
        }
      }
      if (bestMatch < 0){
        numEvents->Fill(72);
        return false;
      }
      else if (phys.genPart_pt().at(bestMatch) > 40){
        numEvents->Fill(72);
        return false; 
      }
    } 
  }

  if ( TString(conf->get("data_set")).Contains("FSMC-TTBar-NoPromptGamma") ){
    
    //Remove prompt photons from TTBar
    if( TString(currentFile->GetTitle()).Contains("ttbar") ){
      if( phys.ngamma() > 0 && (phys.gamma_genIsPromptFinalState().at(0) == 1 && phys.gamma_mcMatchId().at(0) == 22)) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(64);
        num_events_veto_ttbar++;
        return false;
      }
    }   
  }

  if ( TString(conf->get("data_set")).Contains("FSMC-TTGamma-NoNonPromptGamma") ){
    
    //Remove Non-prompt from TTGamma
    if ( TString(currentFile->GetTitle()).Contains("ttgamma_incl_amcnlo") ){
      if( phys.ngamma() > 0 && (phys.gamma_genIsPromptFinalState().at(0) != 1 || phys.gamma_mcMatchId().at(0) != 22 ) ) {
        //cout<<"skipped"<<endl;
        numEvents->Fill(64);
        num_events_veto_ttgamma++;
        return false;
      }
    }   
  }

  return true;
}

//=============================
// Setup
//=============================

void setLepIndexes(){
  /* Loops through lepton objects and adds indexed to g_lep_inds if they pass the tight selection (or fakable object when we are doing fake rate study)*/
  bool FRS = (conf->get("fakerate_study") == "true") ? true : false;

  for (short i = 0; i < (short) phys.lep_p4().size(); i++){
    if(phys.lep_pass_VVV_cutbased_tight().at(i))                g_lep_inds.push_back(i);
    else if (FRS && phys.lep_pass_VVV_cutbased_fo().at(i))      g_lep_inds.push_back(i); 
  }
}

bool isOverlapJet(const LorentzVector &jet_p4){
  /*Takes in a p4 for a jet and determines whether that jet is less than MAX_DR_JET_LEP_OVERLAP for all leptons*/ 
  //cout<<__LINE__<<endl;
  for (int i = 0; i< (int)g_lep_inds.size(); i++){
    if (DeltaR(jet_p4, phys.lep_p4().at(g_lep_inds.at(i))) > MAX_DR_JET_LEP_OVERLAP){
      //cout<<"Failed JET at (eta, phi, pt) = ("<<jet_p4.eta()<<", "<<jet_p4.phi()<<", "<<jet_p4.pt()<<") lep at ("<<phys.lep_p4().at(i).eta()<<", "<<phys.lep_p4().at(i).phi()<<", "<<phys.lep_p4().at(i).pt()<<") DR="<<DeltaR(jet_p4, phys.lep_p4().at(i))<<endl;
      return false;
    }
  }
  //cout<<__LINE__<<endl;

  return true;
}

void writeCleanedJets(const vector<LorentzVector> &vecs, const vector<float> &csvs){
  /*This function writes only elements of the given vector to g_jets_p4 variable, it can be passed the up and down variations as well.*/
  int n_pass = 0;
  //cout<<__LINE__<<endl;
  for (int i = 0; i < (int) vecs.size(); i++){
    if (vecs.at(i).pt() < JET_PT_MIN) continue;
    if (fabs(vecs.at(i).eta()) > JET_ETA_MAX) continue;
    if ( !isOverlapJet(vecs.at(i)) ) { 
      g_jets_p4.push_back(vecs.at(i)); 
      g_jets_csv.push_back(csvs.at(i));
      n_pass++;
    }
  }
  //cout<<__LINE__<<endl;

  g_njets = (int) g_jets_p4.size();
  //cout<<"g_njets: "<<g_njets<<endl;
}

void writeCleanedBJets(const vector<LorentzVector> &vecs, const vector<float> &csvs){
  /*This function writes only elements of the given vector to g_jets_medb_p4 variable, it can be passed the up and down variations as well. Distingushed from the regular jet function because we need to keep track of the CSV values as well.*/
  g_nBJetLoose = 0;
  //cout<<__LINE__<<endl;
  for (int i = 0; i < (int) vecs.size(); i++){
    if (vecs.at(i).pt() < BJET_PT_MIN) continue;
    if (fabs(vecs.at(i).eta()) > BJET_ETA_MAX) continue;
    if ( !isOverlapJet(vecs.at(i)) ) { 
      if (csvs.at(i) > BJET_CSV_MED) {
        g_jets_medb_p4.push_back(vecs.at(i)); 
        g_nBJetLoose++;
      }
      else if (csvs.at(i) > BJET_CSV_LOOSE){
        g_nBJetLoose++; 
      }

    }
  }
  //cout<<__LINE__<<endl;
  g_nBJetMedium = (int) g_jets_medb_p4.size();
  //cout<<"g_nBJetMedium: "<<g_nBJetMedium<<endl;
  //cout<<"g_jets_csv_size: "<<g_jets_csv.size()<<endl;
}

void setupGlobals(){
  Z_VETO_WINDOW_LOW = (conf->get("z_veto_window_low") == "") ? 80 : stod(conf->get("z_veto_window_low"));
  Z_VETO_WINDOW_HIGH = (conf->get("z_veto_window_high") == "") ? 100 : stod(conf->get("z_veto_window_high"));
  W_JET_WINDOW_LOW = (conf->get("w_jet_window_low") == "") ? 60 : stod(conf->get("w_jet_window_low"));
  W_JET_WINDOW_HIGH = (conf->get("w_jet_window_high") == "") ? 100 : stod(conf->get("w_jet_window_high"));
  MAX_DR_JET_LEP_OVERLAP = (conf->get("max_dr_jet_lep") == "") ? 0.4 : stod(conf->get("max_dr_jet_lep"));
  JET_PT_MIN = (conf->get("jet_pt_min") == "") ? 20 : stod(conf->get("jet_pt_min"));
  JET_ETA_MAX = (conf->get("jet_eta_max") == "") ? 2.5 : stod(conf->get("jet_eta_max"));
  BJET_PT_MIN = (conf->get("bjet_pt_min") == "") ? 20 : stod(conf->get("bjet_pt_min"));
  BJET_ETA_MAX = (conf->get("bjet_eta_max") == "") ? 2.4 : stod(conf->get("bjet_eta_max"));
  
  g_jets_p4.clear();
  g_jets_medb_p4.clear();
  g_jets_csv.clear();
  g_lep_inds.clear();

  if ( conf->get("uncertainty_mode") == "JES_up" ){
    g_dphi_metj1 = phys.dphi_metj1_up();
    g_dphi_metj2 = phys.dphi_metj2_up();
    g_mbb = phys.mbb_csv_up();
    g_mjj_mindphi = phys.mjj_mindphi_up();
    g_nBJetMedium = phys.nBJetMedium_up();
    g_met = phys.met_T1CHS_miniAOD_CORE_up_pt();
    g_met_phi = phys.met_T1CHS_miniAOD_CORE_up_phi();
    g_mt2 = phys.mt2_up();
    g_mt2b = phys.mt2b_up();
    g_ht = phys.ht_up();

    /*g_jets_p4 = phys.jets_up_p4();
    g_njets = phys.njets_up();*/
    /*g_nBJetMedium = phys.nBJetMedium_up();
    g_nBJetLoose = phys.nBJetLoose_up();
    g_jets_medb_p4 = phys.jets_medb_up_p4();
    g_jets_csv = phys.jets_up_csv();*/

    writeCleanedJets(phys.jets_up_p4(), phys.jets_csv()); //g_jets_p4, g_njets
    writeCleanedBJets(phys.jets_medb_up_p4(), phys.jets_up_csv());
  }
  else if (conf->get("uncertainty_mode") == "JES_dn"){
    g_dphi_metj1 = phys.dphi_metj1_dn();
    g_dphi_metj2 = phys.dphi_metj2_dn();
    g_mbb = phys.mbb_csv_dn();
    g_mjj_mindphi = phys.mjj_mindphi_dn();
    g_nBJetMedium = phys.nBJetMedium_dn();
    g_met = phys.met_T1CHS_miniAOD_CORE_dn_pt();
    g_met_phi = phys.met_T1CHS_miniAOD_CORE_dn_phi();
    g_mt2 = phys.mt2_dn();
    g_mt2b = phys.mt2b_dn();
    g_ht = phys.ht_dn();

    /*g_jets_p4 = phys.jets_dn_p4();
    g_njets = phys.njets_dn();*/
    /*g_nBJetMedium = phys.nBJetMedium_dn();
    g_nBJetLoose = phys.nBJetLoose_dn();
    g_jets_medb_p4 = phys.jets_medb_dn_p4();
    g_jets_csv = phys.jets_dn_csv();*/

    writeCleanedJets(phys.jets_dn_p4(), phys.jets_csv()); //g_jets_p4, g_njets
    writeCleanedBJets(phys.jets_medb_dn_p4(), phys.jets_dn_csv());
  }
  else if (conf->get("uncertainty_mode") == "GenMet"){
    g_dphi_metj1 = phys.dphi_genmetj1();
    g_dphi_metj2 = phys.dphi_genmetj2();
    g_met = phys.met_genPt();
    g_met_phi = phys.met_genPhi();
    g_mt2 = phys.mt2_genmet();
    g_mt2b = phys.mt2b_genmet();
    g_mjj_mindphi = phys.mjj_mindphi();
    g_mbb = phys.mbb_csv();
    g_ht = phys.ht();

    /*g_jets_p4 = phys.jets_p4();
    g_njets = phys.njets();*/
    /*g_nBJetMedium = phys.nBJetMedium();
    g_nBJetLoose = phys.nBJetLoose();
    g_jets_medb_p4 = phys.jets_medb_p4();
    g_jets_csv = phys.jets_csv();*/

    writeCleanedJets(phys.jets_p4(), phys.jets_csv()); //g_jets_p4, g_njets
    writeCleanedBJets(phys.jets_medb_p4(), phys.jets_csv());
  }
  else{
    g_dphi_metj1 = phys.dphi_metj1();
    g_dphi_metj2 = phys.dphi_metj2();
    g_mbb = phys.mbb_csv();
    g_mjj_mindphi = phys.mjj_mindphi();
    g_met = phys.met_T1CHS_miniAOD_CORE_pt();
    g_met_phi = phys.met_T1CHS_miniAOD_CORE_phi();
    g_mt2 = phys.mt2();
    g_mt2b = phys.mt2b();
    g_ht = phys.ht();
    
    /*g_jets_p4 = phys.jets_p4();
    g_njets = phys.njets();*/
    /*g_nBJetMedium = phys.nBJetMedium();
    g_nBJetLoose = phys.nBJetLoose();
    g_jets_medb_p4 = phys.jets_medb_p4();
    g_jets_csv = phys.jets_csv();*/

    writeCleanedJets(phys.jets_p4(), phys.jets_csv()); //g_jets_p4, g_njets
    writeCleanedBJets(phys.jets_p4(), phys.jets_csv()); //g_jets_medb_p4, g_jets_csv, g_nBJetMedium
  }

  //cout<<__LINE__<<endl;
  
  setLepIndexes();
  if (g_lep_inds.size() >= 2){
    FT = getFlavorType();
    CT = getChargeType();
  }

  //cout<<__LINE__<<endl;
}

int ScanChain( TChain* chain, ConfigParser *configuration, bool fast/* = true*/, int nEvents/* = -1*/) {
  /* Runs through baby files and makes histogram files. 
  
  Inputs:
  chain -- contains the files to make the histograms from, 
  configuration -- pointer to the configuration object
  */  
  //cout<<__LINE__<<endl;
  //Set Global Vars
  conf=configuration;
  //cout<<__LINE__<<endl;
  g_sample_name=conf->get("Name");

  if (conf->get("MCTriggerEmulation") != ""){
    if (conf->get("MCTriggerEmulation") == "true"){
      cout<<"Manually setting MC trigger emulation to true"<<endl;
      MCTriggerEmulation=true;
    }
    else if (conf->get("MCTriggerEmulation") == "false"){
      cout<<"Manually setting MC trigger emulation to false"<<endl;
      MCTriggerEmulation=false;
    }
  }

  TString savePath = getOutputDir(conf, "hist");
  ofstream files_log;
  files_log.open((savePath+TString(g_sample_name+"_files.log")).Data());
  //cout<<__LINE__<<endl;
  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  if (conf->get("scale_factor") != ""){
    g_scale_factor*=stod(conf->get("scale_factor"));
  }

//===========================================
// Define Histograms
//===========================================

  clear_list(); //Event duplicate removal clear list

  cout<<"Opening file "<<TString(savePath+conf->get("Name")+".root")<<endl;
  TFile * output = new TFile(TString(savePath+conf->get("Name")+".root"), "recreate");
  output->cd();

  numEvents = new TH1I("numEvents", "Number of events in "+g_sample_name, 80, 0, 80);
  numEvents->SetDirectory(rootdir);

  const int n_weight_log_bins = 54;
  const double weight_log_bins[n_weight_log_bins+1] = {-5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, -.09, -.08, -.07, -.06, -.05, -.04, -.03, -.02, -.01, 0, .01, .02, .03, .04, .05, .06, .07, .08, .09, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5};

  TH1D* weight_log = new TH1D("weight_log", "Event weights in "+g_sample_name, n_weight_log_bins , weight_log_bins);
  weight_log->SetDirectory(rootdir);

  TH1D* weight_log_flat = new TH1D("weight_log_flat", "Event weights in "+g_sample_name, 101 , 0, 1.01);
  weight_log_flat->SetDirectory(rootdir);

  //MET Histos
  TH1D *t1met = new TH1D("type1MET", "Type 1 MET for "+g_sample_name, 6000,0,6000);
  t1met->SetDirectory(rootdir);
  t1met->Sumw2();

  const int n_metbins_wide_std = 6;
  const double metbins_wide_std[n_metbins_wide_std+1] = {0, 50, 100, 150, 225, 300, 500};

  TH1D *t1met_widebin = new TH1D("type1MET_widebin", "Type 1 MET for "+g_sample_name, n_metbins_wide_std, metbins_wide_std);
  t1met_widebin->SetDirectory(rootdir);
  t1met_widebin->Sumw2();

  //MET Histos
  TH1D *nVert = new TH1D("nVert", "Number of verticies in "+g_sample_name, 150,0,150);
  nVert->SetDirectory(rootdir);
  nVert->Sumw2();

  TH1D *mt2 = new TH1D("mt2", "MT2 for "+g_sample_name, 500,0,500);
  mt2->SetDirectory(rootdir);
  mt2->Sumw2();

  TH1D *dilmass_zlike = new TH1D("dilmass_zlike", "Dilepton Mass for "+g_sample_name, 500,0,500);
  dilmass_zlike->SetDirectory(rootdir);
  dilmass_zlike->Sumw2();

  TH1D *dilmass_zless = new TH1D("dilmass_zless", "Dilepton Mass for "+g_sample_name, 500,0,500);
  dilmass_zless->SetDirectory(rootdir);
  dilmass_zless->Sumw2();

  TH1D *dilmass_sum = new TH1D("dilmass_sum", "Dilepton Mass for "+g_sample_name, 500,0,500);
  dilmass_sum->SetDirectory(rootdir);
  dilmass_sum->Sumw2();

  TH1D *trilep_mass = new TH1D("trilep_mass", "Three Lepton System Mass for "+g_sample_name, 500,0,500);
  trilep_mass->SetDirectory(rootdir);
  trilep_mass->Sumw2();

  TH1D *mt2b = new TH1D("mt2b", "MT2b for "+g_sample_name, 6000,0,6000);
  mt2b->SetDirectory(rootdir);
  mt2b->Sumw2();

  TH1D *nlep = new TH1D("nlep", "Number of Leptons for "+g_sample_name, 20,0,20);
  nlep->SetDirectory(rootdir);
  nlep->Sumw2();

  TH1D *nisotrack = new TH1D("nisotrack", "Number of Iso Track Leptons (MT2 style) for "+g_sample_name, 20,0,20);
  nisotrack->SetDirectory(rootdir);
  nisotrack->Sumw2();

  TH1D *dphi_jet1_met = new TH1D("dphi_jet1_met", "#Delta#Phi(jet_{1}, E^{miss}_{T}) for "+g_sample_name, 100,0,3.15);
  dphi_jet1_met->SetDirectory(rootdir);
  dphi_jet1_met->Sumw2();

  TH1D *dphi_jet2_met = new TH1D("dphi_jet2_met", "#Delta#Phi(jet_{2}, E^{miss}_{T}) for "+g_sample_name, 100,0,3.15);
  dphi_jet2_met->SetDirectory(rootdir);
  dphi_jet2_met->Sumw2();

  TH1D *ht = new TH1D("ht", "Scalar sum of hadronic pt (HT) for "+g_sample_name, 6000,0,6000);
  ht->SetDirectory(rootdir);
  ht->Sumw2();

  TH1D *ht_wide = new TH1D("ht_wide", "Scalar sum of hadronic pt (HT) for "+g_sample_name, 60,0,6000);
  ht_wide->SetDirectory(rootdir);
  ht_wide->Sumw2();

  TH1D *gen_ht = new TH1D("genht", "Scalar sum of generated hadronic pt (Gen HT) for "+g_sample_name, 6000,0,6000);
  gen_ht->SetDirectory(rootdir);
  gen_ht->Sumw2();

  TH1D *numMETFilters = new TH1D("numMETFilters", "Number of MET Filters passed for events in "+g_sample_name, 50,0,50);
  numMETFilters->SetDirectory(rootdir);
  numMETFilters->Sumw2();

  const int n_ptbins_std = 10;
  const double ptbins_std[n_ptbins_std+1] = {0, 22, 33, 40, 55, 85, 105, 135, 180, 250, 6000};

  const int n_ptbins_fine = 51;
  const double ptbins_fine[n_ptbins_fine+1] = {0, 22, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 300, 350, 400, 6001};

  TH1D *vpt;
  if (conf->get("signal_region") == "2j"){
    vpt = new TH1D("vpt", "Boson Pt for events in "+g_sample_name, n_ptbins_fine, ptbins_fine);
  }  
  else{
    vpt = new TH1D("vpt", "Boson Pt for events in "+g_sample_name, n_ptbins_std, ptbins_std);
  }

  vpt->SetDirectory(rootdir);
  vpt->Sumw2();

  TH1D *vpt_flat = new TH1D("vpt_flat", "Boson P_{T} for events in "+g_sample_name, 6000,0,6000);
  vpt_flat->SetDirectory(rootdir);
  vpt_flat->Sumw2();

  TH1D *njets = new TH1D("njets", "Number of jets for events in "+g_sample_name, 50,0,50);
  njets->SetDirectory(rootdir);
  njets->Sumw2();

  TH1D *nbtags_m = new TH1D("nbtags_m", "Number of \"medium\" B-tagged jets for events in "+g_sample_name, 50,0,50);
  nbtags_m->SetDirectory(rootdir);
  nbtags_m->Sumw2();

  TH1D *nbtags_l = new TH1D("nbtags_l", "Number of \"loose\" B-tagged jets for events in "+g_sample_name, 50,0,50);
  nbtags_l->SetDirectory(rootdir);
  nbtags_l->Sumw2();

  TH1D *nbtags_t = new TH1D("nbtags_t", "Number of \"tight\" B-tagged jets for events in "+g_sample_name, 50,0,50);
  nbtags_t->SetDirectory(rootdir);
  nbtags_t->Sumw2();

  TH1D* mjj_min_dphi = new TH1D("mjj_min_dphi", "M_{jj} (for minimum #Delta #Phi Jets) in "+g_sample_name, 6000,0,6000);
  mjj_min_dphi->SetDirectory(rootdir);
  mjj_min_dphi->Sumw2();

  TH1D* mjj_zlike = new TH1D("mjj_zlike", "M_{jj} for jets closest to Z Mass in "+g_sample_name, 6000,0,6000);
  mjj_zlike->SetDirectory(rootdir);
  mjj_zlike->Sumw2();

  TH1D* dEta_jj_zlike = new TH1D("dEta_jj_zlike", "#Delta#Eta_{jj} for jets closest to Z Mass in "+g_sample_name, 700,0,7);
  dEta_jj_zlike->SetDirectory(rootdir);
  dEta_jj_zlike->Sumw2();

  TH1D* mjj_wlike = new TH1D("mjj_wlike", "M_{jj} for jets closest to W Mass in "+g_sample_name, 6000,0,6000);
  mjj_wlike->SetDirectory(rootdir);
  mjj_wlike->Sumw2();

  TH1D* dEta_jj_wlike = new TH1D("dEta_jj_wlike", "#Delta#Eta_{jj} for jets closest to W Mass in "+g_sample_name, 700,0,7);
  dEta_jj_wlike->SetDirectory(rootdir);
  dEta_jj_wlike->Sumw2();

  TH1D* dPhi_3l_MET = new TH1D("dPhi_3l_MET", "#Delta#Phi between 3 lepton system and MET"+g_sample_name, 100,0,3.15);
  dPhi_3l_MET->SetDirectory(rootdir);
  dPhi_3l_MET->Sumw2();

  TH1D *lep1_pt = new TH1D("lep1_pt", "p_{T} for leading lepton "+g_sample_name, 500,0,500);
  lep1_pt->SetDirectory(rootdir);
  lep1_pt->Sumw2();

  TH1D *lep2_pt = new TH1D("lep2_pt", "p_{T} for subleading lepton "+g_sample_name, 500,0,500);
  lep2_pt->SetDirectory(rootdir);
  lep2_pt->Sumw2();

  TH1D *lep3_pt = new TH1D("lep3_pt", "p_{T} for trailing lepton "+g_sample_name, 500,0,500);
  lep3_pt->SetDirectory(rootdir);
  lep3_pt->Sumw2();

  TH1D *lep1_eta = new TH1D("lep1_eta", "#eta for leading lepton "+g_sample_name, 100,0,3);
  lep1_eta->SetDirectory(rootdir);
  lep1_eta->Sumw2();

  TH1D *lep2_eta = new TH1D("lep2_eta", "#eta for subleading lepton "+g_sample_name, 100,0,3);
  lep2_eta->SetDirectory(rootdir);
  lep2_eta->Sumw2();

  TH1D *lep3_eta = new TH1D("lep3_eta", "#eta for trailing lepton "+g_sample_name, 100,0,3);
  lep3_eta->SetDirectory(rootdir);
  lep3_eta->Sumw2();

  TH1D *jet1_pt = new TH1D("jet1_pt", "p_{T} for leading jet "+g_sample_name, 500,0,500);
  jet1_pt->SetDirectory(rootdir);
  jet1_pt->Sumw2();

  TH1D *jet2_pt = new TH1D("jet2_pt", "p_{T} for subleading jet "+g_sample_name, 500,0,500);
  jet2_pt->SetDirectory(rootdir);
  jet2_pt->Sumw2();

  TH1D *jet3_pt = new TH1D("jet3_pt", "p_{T} for trailing jet "+g_sample_name, 500,0,500);
  jet3_pt->SetDirectory(rootdir);
  jet3_pt->Sumw2();

  TH1D *mt2_jl = new TH1D("mt2_jl", "MT2 (Jet+lepton sum) for "+g_sample_name, 6000,0,6000);
  mt2_jl->SetDirectory(rootdir);
  mt2_jl->Sumw2();

  TH1D *lep_signs = new TH1D("lep_signs", "Product of lepton signs"+g_sample_name, 2,-1.01,1.01);
  lep_signs->SetDirectory(rootdir);
  lep_signs->Sumw2();

  TH1D *lep_flavor = new TH1D("lep_flavor", "dil_flavor"+g_sample_name, 0,0,0);
  lep_flavor->SetDirectory(rootdir);
  lep_flavor->Sumw2();

  TH1D *lep1_ip3d = new TH1D("lep1_ip3d", "Leading lepton ip3d for "+g_sample_name, 6000,0,6);
  lep1_ip3d->SetDirectory(rootdir);
  lep1_ip3d->Sumw2();

  TH1D *lep1_ip3derr = new TH1D("lep1_ip3derr", "Leading lepton ip3d error for "+g_sample_name, 6000,0,6);
  lep1_ip3derr->SetDirectory(rootdir);
  lep1_ip3derr->Sumw2();

  TH1D *lep1_sip3d = new TH1D("lep1_sip3d", "Leading lepton sip3d for "+g_sample_name, 6000,0,6);
  lep1_sip3d->SetDirectory(rootdir);
  lep1_sip3d->Sumw2();

  TH1D *lep2_ip3d = new TH1D("lep2_ip3d", "Subleading lepton ip3d for "+g_sample_name, 6000,0,6);
  lep2_ip3d->SetDirectory(rootdir);
  lep2_ip3d->Sumw2();

  TH1D *lep2_ip3derr = new TH1D("lep2_ip3derr", "Subleading lepton ip3d error for "+g_sample_name, 6000,0,6);
  lep2_ip3derr->SetDirectory(rootdir);
  lep2_ip3derr->Sumw2();

  TH1D *lep2_sip3d = new TH1D("lep2_sip3d", "Subleading lepton sip3d for "+g_sample_name, 6000,0,6);
  lep2_sip3d->SetDirectory(rootdir);
  lep2_sip3d->Sumw2();

  TH1D *lep3_ip3d = new TH1D("lep3_ip3d", "Trailing lepton ip3d for "+g_sample_name, 6000,0,6);
  lep1_ip3d->SetDirectory(rootdir);
  lep1_ip3d->Sumw2();

  TH1D *lep3_ip3derr = new TH1D("lep3_ip3derr", "Trailing lepton ip3d error for "+g_sample_name, 6000,0,6);
  lep3_ip3derr->SetDirectory(rootdir);
  lep3_ip3derr->Sumw2();

  TH1D *lep3_sip3d = new TH1D("lep3_sip3d", "Trailing lepton sip3d for "+g_sample_name, 6000,0,6);
  lep3_sip3d->SetDirectory(rootdir);
  lep3_sip3d->Sumw2();

  //-------------------------------------------
  //Truth level tagging W leps ip distributions

  TH1D *Wleps_ip3d = new TH1D("Wleps_ip3d", "Ip3d for leptons tagged as from a W"+g_sample_name, 6000,0,6);
  Wleps_ip3d->SetDirectory(rootdir);
  Wleps_ip3d->Sumw2();

  TH1D *otherleps_ip3d = new TH1D("otherleps_ip3d", "Ip3d for leptons not tagged as from a W"+g_sample_name, 6000,0,6);
  otherleps_ip3d->SetDirectory(rootdir);
  otherleps_ip3d->Sumw2();

  TH1D *Wleps_ip3derr = new TH1D("Wleps_ip3derr", "Ip3derr for leptons tagged as from a W"+g_sample_name, 6000,0,6);
  Wleps_ip3derr->SetDirectory(rootdir);
  Wleps_ip3derr->Sumw2();

  TH1D *Wleps_reliso04 = new TH1D("Wleps_reliso04", "RelIso04 for leptons tagged as from a W"+g_sample_name, 6000,0,6);
  Wleps_reliso04->SetDirectory(rootdir);
  Wleps_reliso04->Sumw2();

  TH1D *Wleps_ptRatio = new TH1D("Wleps_ptRatio", "Pt Ratio for leptons tagged as from a W"+g_sample_name, 6000,0,6);
  Wleps_ptRatio->SetDirectory(rootdir);
  Wleps_ptRatio->Sumw2();

  TH1D *otherleps_ip3derr = new TH1D("otherleps_ip3derr", "Ip3derr for leptons not tagged as from a W"+g_sample_name, 6000,0,6);
  otherleps_ip3derr->SetDirectory(rootdir);
  otherleps_ip3derr->Sumw2();

  TH1D *Wleps_sip3d = new TH1D("Wleps_sip3d", "Sip3d for leptons tagged as from a W"+g_sample_name, 6000,0,6);
  Wleps_sip3d->SetDirectory(rootdir);
  Wleps_sip3d->Sumw2();

  TH1D *otherleps_sip3d = new TH1D("otherleps_sip3d", "Sip3d for leptons not tagged as from a W"+g_sample_name, 6000,0,6);
  otherleps_sip3d->SetDirectory(rootdir);
  otherleps_sip3d->Sumw2();

  TH1D *otherleps_reliso04 = new TH1D("otherleps_reliso04", "RelIso04 for leptons not tagged as from a W"+g_sample_name, 6000,0,6);
  otherleps_reliso04->SetDirectory(rootdir);
  otherleps_reliso04->Sumw2();

  TH1D *otherleps_ptRatio = new TH1D("otherleps_ptRatio", "Pt Ratio for leptons not tagged as from a W"+g_sample_name, 6000,0,6);
  otherleps_ptRatio->SetDirectory(rootdir);
  otherleps_ptRatio->Sumw2();


  cout<<"Histograms initialized"<<endl;
  //cout<<__LINE__<<endl;
//===========================================
// Setup Stuff Pulled From External Files
//===========================================
  int eventsInFile;
  //Set up manual vertex reweighting.  
  if( conf->get("reweight") == "true" ){
    readyReweightHists();
  }
  if( conf->get("vpt_reweight") == "true" ){
    readyVPTReweight(savePath);
  }

  if(conf->get("pileup_reweight") == "true"){
    cout<<"Pileup reweighting with puWeight_Moriond2017.root"<<endl;
    g_pileup_hist_file = TFile::Open("auxFiles/puWeight_Moriond2017.root", "READ");
    //cout<<__LINE__<<endl;
    g_pileup_hist = (TH1D*)g_pileup_hist_file->Get("pileupWeight")->Clone("h_pileup_weight");
    //cout<<__LINE__<<endl;
    g_pileup_hist->SetDirectory(rootdir);
    //cout<<__LINE__<<endl;
    g_pileup_hist_file->Close();
  }

  //--------------------------------------------------------
  // 2D SUSY Scan ISR and BTag SF normalization Histograms
  //--------------------------------------------------------
  if(conf->get("susy_mc") == "true"){
    cout<<"Setting up normalization weights for ISR and Btag Scale Factors."<<endl;
    //cout<<__LINE__<<endl;
    if (conf->get("data_set") == "T5ZZ"){
      //cout<<__LINE__<<endl;
      g_SUSYsf_norm_file = TFile::Open("auxFiles/nsig_weights_t5zz.root", "READ");
    }
    else if (conf->get("data_set") == "TChiWZ"){
      //cout<<__LINE__<<endl;
      g_SUSYsf_norm_file = TFile::Open("auxFiles/nsig_weights_tchiwz.root", "READ");
    }
    else if (conf->get("data_set") == "TChiHZ"){
      //cout<<__LINE__<<endl;
      g_SUSYsf_norm_file = TFile::Open("auxFiles/nsig_weights_tchihz.root", "READ");
    }
    else if (conf->get("data_set") == "TChiHZ_TChiZZ"){
      //cout<<__LINE__<<endl;
      g_SUSYsf_norm_file = TFile::Open("auxFiles/nsig_weights_tchihz.root", "READ");
    }
    else if (conf->get("data_set") == "TChiZZ"){
      //cout<<__LINE__<<endl;
      g_SUSYsf_norm_file = TFile::Open("auxFiles/nsig_weights_tchizz.root", "READ");
    }
    else {
      //cout<<__LINE__<<endl;
      std::stringstream message;
      message<<"Can not pull normalization weights file for "<<conf->get("data_set")<<", no file configured for that dataset.";
      throw std::invalid_argument(message.str());
    }
    //cout<<__LINE__<<endl;

    g_isr_norm = (TH2D*)g_SUSYsf_norm_file->Get("h_avg_weight_isr")->Clone("h_isr_norm");
    //cout<<__LINE__<<endl;
    g_isr_norm_up = (TH2D*)g_SUSYsf_norm_file->Get("h_avg_weight_isr_UP")->Clone("h_isr_norm_up");
    //cout<<__LINE__<<endl;
    g_btagsf_norm = (TH2D*)g_SUSYsf_norm_file->Get("h_avg_weight_btagsf")->Clone("g_btagsf_norm");
    //cout<<__LINE__<<endl;
    g_btagsf_light_norm_up = (TH2D*)g_SUSYsf_norm_file->Get("h_avg_weight_btagsf_light_UP")->Clone("g_btagsf_light_norm_up");
    //cout<<__LINE__<<endl;
    g_btagsf_heavy_norm_up = (TH2D*)g_SUSYsf_norm_file->Get("h_avg_weight_btagsf_heavy_UP")->Clone("g_btagsf_heavy_norm_up");
    //cout<<__LINE__<<endl;
    

    g_isr_norm->SetDirectory(rootdir);
    //cout<<__LINE__<<endl;
    g_isr_norm_up->SetDirectory(rootdir);
    //cout<<__LINE__<<endl;
    g_btagsf_norm->SetDirectory(rootdir);
    //cout<<__LINE__<<endl;
    g_btagsf_light_norm_up->SetDirectory(rootdir);
    //cout<<__LINE__<<endl;
    g_btagsf_heavy_norm_up->SetDirectory(rootdir);
    //cout<<__LINE__<<endl;


    g_SUSYsf_norm_file->Close();
    //cout<<__LINE__<<endl;
  }

   
  //cout<<__LINE__<<endl;
  /*if( phys.isData() && conf->get("event_type")=="photon" ){
    cout<<"Pileup reweighting with "<<savePath+"L1PrescaleWeight_"+conf->get("signal_region")+".root"<<endl;
    g_l1prescale_file = TFile::Open(savePath+"L1PrescaleWeight_"+conf->get("signal_region")+".root", "READ");
    
    g_l1prescale_hist36 = (TH1D*)g_l1prescale_file->Get("rwt_nVert_HLT_Photon36_R9Id90_HE10_IsoM")->Clone("l1prescaleWeight36");
    g_l1prescale_hist36->SetDirectory(rootdir);

    g_l1prescale_hist30 = (TH1D*)g_l1prescale_file->Get("rwt_nVert_HLT_Photon30_R9Id90_HE10_IsoM")->Clone("l1prescaleWeight30");
    g_l1prescale_hist30->SetDirectory(rootdir);

    g_l1prescale_hist22 = (TH1D*)g_l1prescale_file->Get("rwt_nVert_HLT_Photon22_R9Id90_HE10_IsoM")->Clone("l1prescaleWeight22");
    g_l1prescale_hist22->SetDirectory(rootdir);

    g_l1prescale_file->Close();
  }*/
  //cout<<__LINE__<<endl;

  if( conf->get("rwt_photon_eff") == "true" ){
    cout<<"Reweighting for Effeciency with trigeff_Photon165_zmet2016.root"<<endl;
    TFile weight_eff_file("auxFiles/trigeff_Photon165_zmet2016.root", "READ");
    
    //barrel
    g_pt_eff_barrel = (TEfficiency*)weight_eff_file.Get("h_pt_eb_eff_jetht")->Clone("g_vpt_eff_barrel");
    g_pt_eff_barrel->SetDirectory(rootdir);

    //endcap
    g_pt_eff_endcap = (TEfficiency*)weight_eff_file.Get("h_pt_ee_eff_jetht")->Clone("g_vpt_eff_barrel");
    g_pt_eff_endcap->SetDirectory(rootdir);
    
    weight_eff_file.Close();
  }

  if( conf->get("rwt_muon_eff") == "true" ){
    cout<<"Reweighting for single muon trigger effeciency"<<endl;
    TFile weight_eff_file("/home/users/cwelke/analysis/CMSSW_8_0_22/V08-22-05/ZMET2015/makePlots/trigeff_Muon_pt_2016_withH.root", "READ");
    
    //barrel
    g_pt_eff_barrel = (TEfficiency*)weight_eff_file.Get("h_pt_denom_eb_ele27WPLoose_clone")->Clone("g_pt_eff_barrel");
    g_pt_eff_barrel->SetDirectory(rootdir);

    //endcap
    g_pt_eff_endcap = (TEfficiency*)weight_eff_file.Get("h_pt_denom_ee_ele27WPLoose_clone")->Clone("g_pt_eff_barrel");
    g_pt_eff_endcap->SetDirectory(rootdir);
    
    weight_eff_file.Close();
  }  

  //cout<<__LINE__<<endl;
  //set goodrun list
  if (conf->get("JSON") == "ICHEP"){
    const char* json_file = "auxFiles/golden_json_200716_12p9fb_snt.txt"; // ICHEP
    cout<<"Setting good run list: "<<json_file<<endl;
    set_goodrun_file(json_file);   
  }
  else if ((conf->get("JSON") == "18fb")){
    const char* json_file = "auxFiles/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_unblind18_sorted_snt.txt"; // 18.1 fb
    cout<<"Setting good run list: "<<json_file<<endl;
    set_goodrun_file(json_file);
  }
  else{
    const char* json_file = "auxFiles/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_snt.txt"; // 36.5 fb
    cout<<"Setting good run list: "<<json_file<<endl;
    set_goodrun_file(json_file);
  }

  //cout<<__LINE__<<endl;
  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  int eventCount=0;
  cout<<"DATASET: "<<conf->get("data_set")<<endl;
  if( nEvents >= 0 ) nEventsChain = nEvents;
  //cout<<__LINE__<<endl;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  //cout<<__LINE__<<endl;
  TIter fileIter(listOfFiles);
  //cout<<__LINE__<<endl;
//===========================================
// File Loop
//===========================================
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile file( currentFile->GetTitle() );
    //cout<<__LINE__<<endl;
    TTree *tree = (TTree*)file.Get("t");
    //cout<<__LINE__<<endl;
    if(fast) TTreeCache::SetLearnEntries(10); //What does this do?
    //cout<<__LINE__<<endl;
    if(fast) tree->SetCacheSize(128*1024*1024); //What does this do?
    //cout<<__LINE__<<endl;
    phys.Init(tree); //Loads in all the branches
    //cout<<__LINE__<<endl;
    eventsInFile = 0;
    //cout<<__LINE__<<endl;
    files_log<<"Running over new file: "<<currentFile->GetTitle()<<endl;
    cout<<"Running over new file: "<<currentFile->GetTitle()<<endl;
//===========================================
// Loop over Events in current file
//===========================================
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
      numEvents->Fill(0);
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      phys.GetEntry(event);
      ++nEventsTotal;
      //cout<<__LINE__<<endl;    
      // Progress
      WWW::progress( nEventsTotal, nEventsChain );
      setupGlobals();
      //eventsInFile++;
      //if (eventsInFile > 100) continue;
      //cout<<__LINE__<<endl;
//===========================================
// Debugging And Odd Corrections Before Cuts
//===========================================
      
      // ----------------
      // DEBUG MODE
      // ----------------
      printStats = false;
      printFail = false;

      //if (inspection_set.count(phys.evt()) != 0){
      /*if ( inspection_set_erl.count(make_tuple(phys.evt(), phys.run(), phys.lumi())) != 0){
        cout<<"evt: "<<phys.evt()<<" run: "<<phys.run()<<" lumi: "<<phys.lumi()<<" scale1fb: "<<phys.evt_scale1fb()<<endl;
        inspection_copy.erase(make_tuple(phys.evt(), phys.run(), phys.lumi()));
        printStats=true;
        printFail=true;
      }*/
      /*else{ //Use if you don't want care about events in your list that are not in the other's
        continue;
      }*/
//===========================================
// Cuts
//===========================================
      //cout<<__LINE__<<endl;      
      //Set up event weight
      /*if (event % 10000 == 0){
        cout<<"Weight: "<<weight<<endl;
      }*/

      if ( isDuplicate() ){
        //cout<<"Failed Duplicate"<<endl;
        numEvents->Fill(23);
        continue;
      } // check for duplicates
      //cout<<__LINE__<<endl;      

      if (! passFileSelections() ){
        //cout<<"Failed File Selections"<<endl;
        continue;
      }
      //cout<<__LINE__<<endl;

      if (! passBaseCut()){ 
        //cout<<"Failed Baseline"<<endl;
        continue; 
      }// Base Cut
      //cout<<__LINE__<<endl;      

      if (! hasGoodEvent()){
        //cout<<"Failed Good Event"<<endl;
        continue; // Event Type Specific Cuts
      }
      //cout<<__LINE__<<endl;      

      if (! passSignalRegionCuts()){ 
        //cout<<"Failed SR"<<endl;
        continue; // Signal Region Cuts
      }
      //cout<<__LINE__<<endl;
      
      if (conf->get("do_met_filters") == "true")
      {
        //cout<<"checking MET filters"<<endl;
        if (! passMETFilters()) continue; ///met filters
      }
      
      //double weight=1;
      double weight = getWeight();
      weight_log->Fill(log10(fabs(weight)));
      weight_log_flat->Fill(fabs(weight));

        // ----------------
        // DEBUG MODE
        // ----------------
        /*if (inspection_set_erl.count(make_tuple(phys.evt(), phys.run(), phys.lumi())) == 0){
          cout<<"NEW||evt: "<<phys.evt()<<" run: "<<phys.run()<<" lumi: "<<phys.lumi()<<" scale1fb: "<<phys.evt_scale1fb()<<" weight: "<<weight<<endl;
          //cout<<"Inspection Set Count "<<inspection_set_erl.count(make_tuple(phys.evt(), phys.run(), phys.lumi()))<<endl;
        }*/
        //When Debug mode is off, you can turn this on:
        //cout<<"evt: "<<phys.evt()<<" run: "<<phys.run()<<" lumi: "<<phys.lumi()<<" scale1fb: "<<phys.evt_scale1fb()<<" weight: "<<weight<<" extra_weight: "<< weight/phys.evt_scale1fb() <<endl;
//===========================================
// Analysis Code
//===========================================
      //cout<<__LINE__<<endl;
      //cout<<"Event Weight "<<weight<<endl;      
      //Fill in Histos
      double sumMETFilters = phys.Flag_HBHENoiseFilter()+phys.Flag_HBHEIsoNoiseFilter()+phys.Flag_CSCTightHaloFilter()+phys.Flag_EcalDeadCellTriggerPrimitiveFilter()+phys.Flag_goodVertices()+phys.Flag_eeBadScFilter();
      //cout<<__LINE__<<endl;      
      numMETFilters->Fill(sumMETFilters);

      /*if (weight < 0){
         cout<<"Negative Weight2: "<<weight<<" "<<phys.evt()<<endl;
      }
      */
      if (g_met != 0) {
        t1met->Fill(g_met, weight);
        t1met_widebin->Fill(g_met, weight);
      }
      if (g_ht != 0) {
        ht->Fill(g_ht, weight);
        ht_wide->Fill(g_ht, weight);
      }
      if (phys.gen_ht() != 0) gen_ht->Fill(phys.gen_ht(), weight);
      if (bosonPt() != 0){ 
        vpt->Fill(bosonPt(), weight); 
        vpt_flat->Fill(bosonPt(), weight); 
      }
      njets->Fill(g_njets, weight);
      nbtags_m->Fill(g_nBJetMedium, weight);
      nbtags_l->Fill(phys.nBJetLoose(), weight);
      nbtags_t->Fill(phys.nBJetTight(), weight);
      nVert->Fill(phys.nVert(), weight);
      nlep->Fill(phys.nlep(), weight);
      //cout<<"Filling nisotrack"<<endl;
      nisotrack->Fill(phys.nisoTrack_mt2(), weight);
      //cout<<__LINE__<<endl;
      if (g_mt2 != 0 ) mt2->Fill(g_mt2, weight);
      //cout<<__LINE__<<endl;
      if (g_mt2b != 0 ) mt2b->Fill(g_mt2b, weight);
      //cout<<__LINE__<<endl;
      if (g_njets > 0) dphi_jet1_met->Fill(acos(cos(g_met_phi - g_jets_p4.at(0).phi())), weight);
      //cout<<__LINE__<<endl;
      if (g_njets > 1) dphi_jet2_met->Fill(acos(cos(g_met_phi - g_jets_p4.at(1).phi())), weight);
      //cout<<__LINE__<<endl;
      
      if (conf->get("num_leptons") == "3"){
        pair<int, int> indicies = getMostZlikePair(phys.lep_p4());
        dilmass_zlike->Fill((phys.lep_p4().at(indicies.first) + phys.lep_p4().at(indicies.second)).M(), weight);
        
        indicies = getLeastZlikePair(phys.lep_p4());
        dilmass_zless->Fill((phys.lep_p4().at(indicies.first) + phys.lep_p4().at(indicies.second)).M(), weight);

        double m1=(phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M();
        double m2=(phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).M();
        double m3=(phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(2))).M();
        dilmass_sum->Fill(m1+m2+m3, weight);
        trilep_mass->Fill((phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).M(), weight);

        double trilep_phi = (phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1)) + phys.lep_p4().at(g_lep_inds.at(2))).phi();
        dPhi_3l_MET->Fill(acos(cos(g_met_phi - trilep_phi)), weight);

        lep1_pt->Fill(phys.lep_pt().at(g_lep_inds.at(0)), weight);
        lep2_pt->Fill(phys.lep_pt().at(g_lep_inds.at(1)), weight);
        lep3_pt->Fill(phys.lep_pt().at(g_lep_inds.at(2)), weight);

        lep1_eta->Fill(phys.lep_eta().at(g_lep_inds.at(0)), weight);
        lep2_eta->Fill(phys.lep_eta().at(g_lep_inds.at(1)), weight);
        lep3_eta->Fill(phys.lep_eta().at(g_lep_inds.at(2)), weight);

        int signs=(phys.lep_pdgId().at(g_lep_inds.at(0))/fabs(phys.lep_pdgId().at(g_lep_inds.at(0))));
        signs *= (phys.lep_pdgId().at(g_lep_inds.at(1))/fabs(phys.lep_pdgId().at(g_lep_inds.at(1))));
        signs *= (phys.lep_pdgId().at(g_lep_inds.at(2))/fabs(phys.lep_pdgId().at(g_lep_inds.at(2))));
        lep_signs->Fill(signs, weight);
      }
      else{
        lep1_pt->Fill(phys.lep_pt().at(g_lep_inds.at(0)), weight);
        lep2_pt->Fill(phys.lep_pt().at(g_lep_inds.at(1)), weight);

        lep1_eta->Fill(phys.lep_eta().at(g_lep_inds.at(0)), weight);
        lep2_eta->Fill(phys.lep_eta().at(g_lep_inds.at(1)), weight);

        dilmass_zlike->Fill((phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M(), weight);
        dilmass_zless->Fill((phys.lep_p4().at(g_lep_inds.at(0)) + phys.lep_p4().at(g_lep_inds.at(1))).M(), weight);
        dilmass_sum->Fill(phys.lep_p4().at(g_lep_inds.at(0)).M() + phys.lep_p4().at(g_lep_inds.at(1)).M(), weight);

        int signs=(phys.lep_pdgId().at(g_lep_inds.at(0))/fabs(phys.lep_pdgId().at(g_lep_inds.at(0))));
        signs *= (phys.lep_pdgId().at(g_lep_inds.at(1))/fabs(phys.lep_pdgId().at(g_lep_inds.at(1))));
        lep_signs->Fill(signs, weight);
      }

      lep_flavor->Fill(getLepFlavorString().Data(), weight);
      mjj_min_dphi->Fill(g_mjj_mindphi, weight);


      if (g_njets >= 2){
        pair<int, int> indicies = getMostZlikePair(g_jets_p4);
        mjj_zlike->Fill((g_jets_p4.at(indicies.first) + g_jets_p4.at(indicies.second)).M(), weight);
        dEta_jj_zlike->Fill(fabs(g_jets_p4.at(indicies.first).eta() - g_jets_p4.at(indicies.second).eta()), weight);

        indicies = getMostWlikePair(g_jets_p4);
        mjj_wlike->Fill((g_jets_p4.at(indicies.first) + g_jets_p4.at(indicies.second)).M(), weight);
        dEta_jj_wlike->Fill(fabs(g_jets_p4.at(indicies.first).eta() - g_jets_p4.at(indicies.second).eta()), weight);

        mt2_jl->Fill(getMT2_JL(), weight);
      }

      if (g_njets >= 1){
        jet1_pt->Fill(g_jets_p4.at(0).pt(), weight);
      }
      if (g_njets >= 2){
        jet2_pt->Fill(g_jets_p4.at(1).pt(), weight);
      }
      if (g_njets >= 3){
        jet3_pt->Fill(g_jets_p4.at(2).pt(), weight);
      }

      lep1_ip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(0))), weight);
      lep1_ip3derr->Fill(fabs(phys.lep_ip3derr().at(g_lep_inds.at(0))), weight);
      lep1_sip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(0)))/fabs(phys.lep_ip3derr().at(g_lep_inds.at(0))), weight);
      lep2_ip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(1))), weight);
      lep2_ip3derr->Fill(fabs(phys.lep_ip3derr().at(g_lep_inds.at(1))), weight);
      lep2_sip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(1)))/fabs(phys.lep_ip3derr().at(g_lep_inds.at(1))), weight);

      if (phys.nlep() > 2){
        lep3_ip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(2))), weight);
        lep3_ip3derr->Fill(fabs(phys.lep_ip3derr().at(g_lep_inds.at(2))), weight);
        lep3_sip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(2)))/fabs(phys.lep_ip3derr().at(g_lep_inds.at(2))), weight);
      }

      for (int c = 0; c < (int) g_lep_inds.size(); c++){
        if (phys.lep_isFromW().at(c)){
          Wleps_ip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(c))), weight);
          Wleps_sip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(c)))/fabs(phys.lep_ip3derr().at(g_lep_inds.at(c))), weight);
          Wleps_ip3derr->Fill(fabs(phys.lep_ip3derr().at(g_lep_inds.at(c))), weight);  
          Wleps_reliso04->Fill(fabs(phys.lep_relIso04EA().at(g_lep_inds.at(c))), weight);
          Wleps_ptRatio->Fill(fabs(phys.lep_ptRatio().at(g_lep_inds.at(c))), weight);
        }
        else if (phys.lep_isFromZ().at(g_lep_inds.at(c)) || phys.lep_isFromB().at(g_lep_inds.at(c)) || phys.lep_isFromC().at(g_lep_inds.at(c)) || phys.lep_isFromL().at(g_lep_inds.at(c)) || phys.lep_isFromLF().at(g_lep_inds.at(c)) ){
          otherleps_ip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(c))), weight);
          otherleps_sip3d->Fill(fabs(phys.lep_ip3d().at(g_lep_inds.at(c)))/fabs(phys.lep_ip3derr().at(g_lep_inds.at(c))), weight);
          otherleps_ip3derr->Fill(fabs(phys.lep_ip3derr().at(g_lep_inds.at(c))), weight);
          otherleps_reliso04->Fill(fabs(phys.lep_relIso04EA().at(g_lep_inds.at(c))), weight);
          otherleps_ptRatio->Fill(fabs(phys.lep_ptRatio().at(g_lep_inds.at(c))), weight);
        }
      }
//===========================================
// Signal Region Specific Histos
//===========================================
//===========================================
// Debugging And Odd Corrections After Cuts
//===========================================
      /*if (conf->get("rares") == "true"){
        //cout<<__LINE__<<endl;
        //cout<<"EVENT-LIST "<<eventCount<<" : "<<phys.evt()<<endl;
          //cout<<__LINE__<<endl;
        cout<<"EVENT-LIST "<<eventCount<<" : "<<phys.evt()<<" "<<g_met<<endl;
        eventCount++;
        if ( inVinceNotMine.count(phys.evt()) != 0){
          //printPass = true;
        }
        if ( inMineNotVince.count(phys.evt()) != 0){
          printFail = true;
        }
      }*/
      eventCount++;

      /*if ( inVinceNotMine.count(phys.evt()) != 0){
        printFail = true;
        cout<<"Event passed: "<<phys.evt()<<endl;
      }*/
      
      //if (printStats) {cout<<"Event: "<<phys.evt()<<endl;}
      //cout<<"Event_Run_Lumi: "<<phys.evt()<<" "<<phys.run()<<" "<<phys.lumi()<<endl;
      //if(g_met >= 300){
      //  cout<<"Event: "<<phys.evt()<<" MET: "<<g_met<<" g_njets: "<<g_njets<<" nbtags: "<<g_nBJetMedium<<" HT: "<<g_ht<<endl;
      //}
    }
    // Clean Up
    //cout<<__LINE__<<endl;
    delete tree;
    //cout<<__LINE__<<endl;
    file.Close();
  }

  // ----------------
  // DEBUG MODE
  // ----------------
  /*cout<<"Events that weren't in your babies:"<<endl;
  for (set<tuple<long,long,long>>::iterator it=inspection_copy.begin(); it!=inspection_copy.end(); ++it){
    cout<<"evt: "<<std::get<0>(*it)<<" run: "<<std::get<1>(*it)<<" lumi: "<<std::get<2>(*it)<<endl;
  }*/

  cout<<"Num events passed: "<<eventCount<<endl;
  files_log<<"Num events passed: "<<eventCount<<endl;
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }  
  //cout<<__LINE__<<endl;
//=======================================
// Write Out Histos
//=======================================
  output->cd();

  //Write out histograms to file
  numEvents->Write();
  //cout<<__LINE__<<endl;
  weight_log->Write();
  //cout<<__LINE__<<endl;
  weight_log_flat->Write();
  //cout<<__LINE__<<endl;
  numMETFilters->Write();
  //cout<<__LINE__<<endl;
  t1met->Write();
  //cout<<__LINE__<<endl;
  t1met_widebin->Write();
  //cout<<__LINE__<<endl;
  ht->Write();
  //cout<<__LINE__<<endl;
  ht_wide->Write();
  //cout<<__LINE__<<endl;
  gen_ht->Write();
  //cout<<__LINE__<<endl;
  vpt->Write();
  //cout<<__LINE__<<endl;
  vpt_flat->Write();
  //cout<<__LINE__<<endl;
  njets->Write();
  //cout<<__LINE__<<endl;
  nlep->Write();
  //cout<<__LINE__<<endl;
  nisotrack->Write();
  //cout<<__LINE__<<endl;
  nbtags_m->Write();
  //cout<<__LINE__<<endl;
  nbtags_l->Write();
  //cout<<__LINE__<<endl;
  nbtags_t->Write();
  //cout<<__LINE__<<endl;
  nVert->Write();
  //cout<<__LINE__<<endl;
  mt2->Write();
  //cout<<__LINE__<<endl;
  mt2b->Write();
  //cout<<__LINE__<<endl;
  dphi_jet1_met->Write();
  //cout<<__LINE__<<endl;
  dphi_jet2_met->Write();
  //cout<<__LINE__<<endl;
  mjj_min_dphi->Write();
  //cout<<__LINE__<<endl;
  dilmass_zlike->Write();
  //cout<<__LINE__<<endl;
  dilmass_zless->Write();
  //cout<<__LINE__<<endl;
  dilmass_sum->Write();
  //cout<<__LINE__<<endl;
  trilep_mass->Write();
  //cout<<__LINE__<<endl;
  dPhi_3l_MET->Write();
  //cout<<__LINE__<<endl;
  dEta_jj_zlike->Write();
  //cout<<__LINE__<<endl;
  mjj_zlike->Write();
  //cout<<__LINE__<<endl;
  dEta_jj_wlike->Write();
  //cout<<__LINE__<<endl;
  mjj_wlike->Write();
  //cout<<__LINE__<<endl;
  lep1_pt->Write();
  //cout<<__LINE__<<endl;
  lep2_pt->Write();
  //cout<<__LINE__<<endl;
  lep3_pt->Write();
  //cout<<__LINE__<<endl;
  lep1_eta->Write();
  //cout<<__LINE__<<endl;
  lep2_eta->Write();
  //cout<<__LINE__<<endl;
  lep3_eta->Write();
  //cout<<__LINE__<<endl;
  jet1_pt->Write();
  //cout<<__LINE__<<endl;
  jet2_pt->Write();
  //cout<<__LINE__<<endl;
  jet3_pt->Write();
  //cout<<__LINE__<<endl;
  mt2_jl->Write();
  //cout<<__LINE__<<endl;
  lep_signs->Write();
  //cout<<__LINE__<<endl;
  lep_flavor->Write();
  //cout<<__LINE__<<endl;

  lep1_ip3d->Write();
  //cout<<__LINE__<<endl;
  lep1_ip3derr->Write();
  //cout<<__LINE__<<endl;
  lep1_sip3d->Write();
  //cout<<__LINE__<<endl;
  lep2_ip3d->Write();
  //cout<<__LINE__<<endl;
  lep2_ip3derr->Write();
  //cout<<__LINE__<<endl;
  lep2_sip3d->Write();
  //cout<<__LINE__<<endl;
  lep3_ip3d->Write();
  //cout<<__LINE__<<endl;
  lep3_ip3derr->Write();
  //cout<<__LINE__<<endl;
  lep3_sip3d->Write();
  //cout<<__LINE__<<endl;
  Wleps_ip3d->Write();
  //cout<<__LINE__<<endl;
  Wleps_sip3d->Write();
  //cout<<__LINE__<<endl;
  Wleps_ip3derr->Write();
  //cout<<__LINE__<<endl;
  Wleps_reliso04->Write();
  //cout<<__LINE__<<endl;
  Wleps_ptRatio->Write();
  //cout<<__LINE__<<endl;
  otherleps_ip3d->Write();
  //cout<<__LINE__<<endl;
  otherleps_ip3derr->Write();
  //cout<<__LINE__<<endl;
  otherleps_sip3d->Write();
  //cout<<__LINE__<<endl;
  otherleps_reliso04->Write();
  //cout<<__LINE__<<endl;
  otherleps_ptRatio->Write();

  //close output file
  output->Write();
  output->Close();
  g_reweight_pairs.clear();
  files_log.close();
  //cout<<__LINE__<<endl;
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << num_events_veto_ttbar << " Events Vetod from TTBar" << endl;
  cout << num_events_veto_ttgamma << " Events Vetod from TTGamma" << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << nDuplicates << " Duplicates" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  
  //====================
  // Clean Everything Up
  //====================
  delete bmark;
  /*delete output;
  delete numEvents;
  delete weight_log;
  delete weight_log_flat;
  delete numMETFilters;
  delete t1met;
  delete t1met_widebin;
  delete rawmet;
  delete ht;
  delete ht_wide;
  delete gen_ht;
  delete vpt;
  delete vpt_flat;
  delete njets;
  delete nlep;
  delete nisotrack;
  delete nbtags_m;
  delete nbtags_l;
  delete nbtags_t;
  delete nVert;
  delete mt2;
  delete mt2b;
  delete dphi_jet1_met;
  delete dphi_jet2_met;
  delete dilmass;
  if (conf->get("signal_region") == "TChiHZ"){
    delete sum_mlb;
    delete m_bb_csv;
    delete m_bb_bpt;
    delete mt2j;
    delete mt2_fromb;
    delete mt2_hz;
    delete sum_pt_z_bb;
    //2D hists
    delete MT2_MT2B;
    delete MT2_MT2_fromb;
    delete MT2_MT2_HZ;
  }
  if (conf->get("GammaMuStudy") == "true"){
    delete MT_MuMET;
    delete dR_GammaMu;
    delete mu_pt;
  }*/

  return 0;
}
