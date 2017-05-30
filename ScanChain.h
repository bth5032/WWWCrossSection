// Usage:
// > root -b doAll.C

//
// 2016 MET study looper. Written by Bobak Hashemi May 13 2016
//

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

// Analysis Specific
#include "WWW.cc"

// CORE
//You can not include headers!!! This is not compiled code.
#include "External/dorky.cc"
#include "External/goodrun.cc"
#include "External/MT2Utility.cc"
#include "External/MT2.cc"

// Configuration parsing
#include "ConfigParser.C"
#include "ConfigHelper.C"

using namespace std;
//using namespace zmet;
using namespace duplicate_removal;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//Global Vars
ConfigParser *conf;
int nDuplicates=0;
int num_events_veto_ttbar=0;
int num_events_veto_ttgamma=0;
bool MCTriggerEmulation = true;

vector<pair <TH1D*, TString> > g_reweight_pairs;
TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
TH1D *g_pileup_hist, *g_l1prescale_hist22, *g_l1prescale_hist30, *g_l1prescale_hist36; 

//Btag and ISR Scale Factor overall normalization
TH2D *g_btagsf_norm, *g_btagsf_light_norm_up, *g_btagsf_heavy_norm_up;
TH2D *g_isr_norm, *g_isr_norm_up;
TFile *g_SUSYsf_norm_file;

TEfficiency *g_pt_eff_barrel, *g_pt_eff_endcap; 
TFile *g_weight_hist_file, *g_pileup_hist_file, *g_l1prescale_file;
TString g_sample_name;
TFile* currentFile = 0;
double g_scale_factor=1; //Holds scale factors for sample.

TH1I *numEvents; //Holds the number of events in the whole script and the number that pass various cuts 

// ----------------
// DEBUG MODE
// ----------------

set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(1600134,1,11269), make_tuple(944882,1,6655), make_tuple(1215665,1,8561), make_tuple(3793537,1,26715), make_tuple(3161564,1,22265), make_tuple(3065612,1,21589), make_tuple(3922112,1,27621), make_tuple(3867464,1,27236), make_tuple(393193,1,2769), make_tuple(2698570,1,19005), make_tuple(3895906,1,27436), make_tuple(12026619,1,20033), make_tuple(8316210,1,13853), make_tuple(20295173,1,33806), make_tuple(3129628,1,5213), make_tuple(8703007,1,14497), make_tuple(14958843,1,24917), make_tuple(518008,1,863), make_tuple(20250631,1,33733), make_tuple(14651822,1,24406), make_tuple(21849611,1,36396), make_tuple(24391258,1,40630), make_tuple(9230235,1,15376), make_tuple(24297457,1,40474), make_tuple(1453059,1,2420), make_tuple(2533975,1,4221), make_tuple(16039218,1,26717), make_tuple(5425457,1,9038), make_tuple(25306736,1,42155), make_tuple(26062108,1,43413), make_tuple(1062197,1,1769), make_tuple(16653930,1,27742), make_tuple(3545054,1,5905), make_tuple(13610367,1,22672), make_tuple(12965144,1,21596), make_tuple(22355385,1,37238), make_tuple(14025939,1,23363), make_tuple(10900371,1,18157), make_tuple(6866420,1,11438), make_tuple(7111057,1,11845), make_tuple(12705532,1,21164), make_tuple(13737873,1,22884), make_tuple(21828665,1,36361), make_tuple(20920805,1,34849), make_tuple(13543401,1,22560), make_tuple(22070301,1,36764), make_tuple(12355981,1,20582), make_tuple(14453977,1,24077), make_tuple(467990,1,779), make_tuple(21569057,1,35929), make_tuple(7968614,1,13274), make_tuple(21206179,1,35324), make_tuple(18871114,1,31435), make_tuple(9012627,1,15013), make_tuple(7003001,1,11665), make_tuple(26700794,1,44477), make_tuple(20125093,1,33523), make_tuple(5005685,1,8338), make_tuple(23145096,1,38554), make_tuple(22275218,1,37105), make_tuple(836965,1,1741), make_tuple(33130670,1,68879), make_tuple(24565216,1,51071), make_tuple(6424757,1,13357), make_tuple(27697266,1,57583), make_tuple(30227066,1,62842), make_tuple(1340659,1,2787), make_tuple(2098427,1,4363), make_tuple(8252029,1,17156), make_tuple(6466834,1,13445), make_tuple(4068934,1,8459), make_tuple(18928244,1,39351), make_tuple(15765673,1,32777), make_tuple(26246814,1,54567)};

//set<long> inspection_set = {99795815,998615983,998751102,999957050};

set<tuple<long,long,long>> inspection_copy = inspection_set_erl;

bool printStats = false;
bool printFail = false;

//=======================================================
// Global variables used for uncertainty fluctuations
//=======================================================
double g_dphi_metj1;
double g_dphi_metj2;
int g_njets;
double g_mbb;
double g_mjj_mindphi;
int g_nBJetMedium;
int g_nBJetLoose;
double g_met;
double g_met_phi;
double g_mt2;
double g_mt2b;
double g_ht;

vector<float> g_jets_csv;
vector<LorentzVector> g_jets_p4;
vector<LorentzVector> g_jets_medb_p4;


const int Z_MASS = 91;
const int W_MASS = 80;

const double BJET_CSV_TIGHT = 0.9535;
const double BJET_CSV_MED = 0.8484;
const double BJET_CSV_LOOSE = 0.5426;

double MAX_DR_JET_LEP_OVERLAP, JET_PT_MIN, JET_ETA_MAX, BJET_PT_MIN, BJET_ETA_MAX;
double W_JET_WINDOW_LOW, W_JET_WINDOW_HIGH, Z_VETO_WINDOW_LOW, Z_VETO_WINDOW_HIGH;

/* returns two most B-like jet indicies */
pair<int, int> getMostBlike();

/*Finds the pair of B tagged jets (csv medium) with dijet mass closest to the mass of the higgs*/
pair<int,int> getClosestBPairToHiggsMass();

/*Builds MT2b from two highest CSV jets*/
double getMT2B();

/*Builds Mbb from two highest CSV jets*/
double getMbb();

/*Builds MT2 (jet lepton sum) from two most W-like jets*/
double getMT2_JL();

/*This function gets the MT2 built out of the two Bjets in an event, no guarentee 
is made about selecting the highest csv jets*/
double getMT2ForBjets(bool select_highest_csv=false);

/*Builds MT2 for the two leading Bjets unless select_closest_higgs_mass is set, in which case it 
builds it out of the two bjets with dijet mass nearest the mass of the higgs.*/
double getMT2HiggsZ(bool select_highest_closest_higgs_mass=false);

/*Returns boson Pt, determines whether sample is gjets or zjets first*/
double bosonPt();

/* Builds the MT from the lepton at index id and the MET vector (assumes massless particles)*/
double getMTLepMET(short id=0);

/* Builds the delta R (sqrt(dPhi^2 + dEta^2)) between the lepton at index id and the leading photon*/
double getdRGammaLep(short id=0);

/*Searches the vector of lorentz vectors for the pair with mass nearest (furthest) from the target mass if 'close' is true (false)*/
pair<int, int> getPairWithMass(const vector<LorentzVector> &vecs, double target_mass, bool close);

/*These get the 4 vectors for objects with mass closest or furthest from the mass of the W or Z using getPairWithMass()*/
pair<int,int> getMostZlikePair(const vector<LorentzVector> &vecs);
pair<int,int> getLeastZlikePair(const vector<LorentzVector> &vecs);
pair<int,int> getMostWlikePair(const vector<LorentzVector> &vecs);

/*Loops through pairs of entries in the lep_pdgId vector and counts how many have opposite value*/
int getNumOSSFPairs();

/*Returns a string which gives the exact flavor composition of the leptons in the event*/
TString getLepFlavorString();

/*Returns the DeltaR between objects p1 and p2.*/
double DeltaR(const LorentzVector p1, const LorentzVector p2);

/* checks the gen record for a lepton with the index specified */
bool isCleanLepFromW(int index);
//=============================
// Triggers
//=============================
/*Checks that the event passes an "emulated photon trigger"*/
bool passPhotonEmulatedTrigger();

/*Ensures the event is within the efficiency plateu of the highest pt trigger it passed*/
bool passPhotonTriggers();

/*MC passes immediately, ensures data events were gathered from di-muon triggers*/
bool passMuonTriggers();

/*MC passes immediately, ensures data events were gathered from di-electron triggers*/
bool passElectronTriggers();

/*MC passes immediately, ensures data events were gathered from EMu triggers*/
bool passEMuTriggers();

/*MC passes immediately, ensures data events were gathered from SingleMu trigger*/
bool passSingleMuTriggers();

/*Helper method which chooses which above method to call. Calls EMu if the dil_flavor is emu, otherwise uses
the hyp_type to determine which to call. Events fail if they are hyp_type 2 and not tagged for emu*/
bool passLeptonHLTs();

//=============================
// Has Good Event Functions
//=============================
/*Lepton quality and Z mass cuts*/
bool hasGoodZ();

/*Photon quality cuts*/
bool hasGoodPhoton();

/*Method for testing whether the event has a good gamma mu pair trigger requirements are on the photon.
  It just checks muon quality stuff and then calls hasGoodPhoton()*/
bool hasGoodGammaMu();

/*Just a helper method that chooses which hasGood method to call based on the config event_type*/
bool hasGoodEvent();

//=============================
// Event Weight Assignment
//=============================

/*Goes through the chain of weight_from config options down to a config which does not have weight_from and
then adds a pair (config_name, hist_file) to the vector g_reweight_pairs.

For now this is depricated: NEEDS TO BE UPDATED WITH NEW CODE FIXES*/
void readyReweightHists();

/* Adds the vpt reweighting histogram to the g_reweight_pairs vector */
void readyVPTReweight(TString save_path);

/* Returns the trigger efficiency from g_pt_eff */
double getEff(const double &pt, const double &eta);

/*Loads the reweight hists from g_reweight_pairs and multiplies returns the weight associated with the proper
bin in the histogram*/
double getReweight();

/*This method stores fixes to the evt_scale1fb in the event of file corruptions. 
It's basically just a lookup table*/
double scale1fbFix();

/*Main function for determining the weights for each event*/
double getWeight();

/*Returns the weight associated with the photon prescales*/
double getPrescaleWeight();

//=============================
// Cuts
//=============================

/*Holds the cuts for all the signal regions, basically all the cuts that are turned on with a config option*/
bool passSignalRegionCuts();

/*Checks for a gen Neutrino (Real MET) and a gen Z (Real Z), only should be run when running
over samples tagged as "rares". This is only neccesary for the full prediction.*/
bool passRareCuts();

/*Front end method to "Dorky" duplicate removal*/
bool isDuplicate();

/*Checks for MET filters*/
bool passMETFilters();

/*Holds baseline cuts*/
bool passBaseCut();

/*Method which holds all the file specific selections, for instance cutting out the
  events with genht > 100 in the DY inclusive samples*/
bool passFileSelections();

//=============================
// Setup
//=============================
/*Takes in a p4 for a jet and determines whether that jet is less than MAX_DR_JET_LEP_OVERLAP for all leptons*/
bool isOverlapJet(const LorentzVector &jet_p4);

/*This function writes only elements of the given vector to g_jets_p4 variable, it can be passed the up and down variations as well.*/
void writeCleanedJets(const vector<LorentzVector> &vecs);

/*This function writes only elements of the given vector to g_jets_medb_p4 variable, it can be passed the up and down variations as well. Distingushed from the regular jet function because we need to keep track of the CSV values as well.*/
void writeCleanedBJets(const vector<LorentzVector> &vecs, const vector<float> &csvs);

/*Sets up global variables for the event which are the quantities that might be fluctuated in the process of computing uncertainty limits*/
void setupGlobals();

/*Obvi the event looper*/
int ScanChain( TChain* chain, ConfigParser *configuration, bool fast = true, int nEvents = -1);