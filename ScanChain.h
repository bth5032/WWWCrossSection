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
#include "TLorentzVector.h"

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

// Utility stuff
#include "Utils.C"

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

//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(40,1005,1), make_tuple(19,112,1), make_tuple(20,1133,1), make_tuple(66,788,2), make_tuple(97,791,2), make_tuple(18,797,1), make_tuple(42,812,1), make_tuple(100,812,2), make_tuple(60,815,2), make_tuple(22,825,1), make_tuple(35,842,1), make_tuple(6,847,1), make_tuple(59,84,2), make_tuple(5,852,1), make_tuple(9,860,1), make_tuple(59,891,2), make_tuple(15,899,1), make_tuple(3,906,1), make_tuple(28,113,1), make_tuple(96,1147,2), make_tuple(28,1148,1), make_tuple(63,1148,2), make_tuple(36,925,1), make_tuple(33,936,1), make_tuple(12,948,1), make_tuple(48,956,1), make_tuple(68,968,2), make_tuple(55,970,2), make_tuple(9,977,1), make_tuple(15,977,1), make_tuple(67,991,2), make_tuple(72,995,2), make_tuple(64,1154,2), make_tuple(81,1156,2), make_tuple(6,1160,1), make_tuple(78,1179,2), make_tuple(39,1186,1), make_tuple(82,1188,2), make_tuple(78,1216,2), make_tuple(12,1223,1), make_tuple(16,1237,1), make_tuple(59,1261,2), make_tuple(8,1262,1), make_tuple(70,1016,2), make_tuple(5,1020,1), make_tuple(69,1021,2), make_tuple(4,1289,1), make_tuple(95,1289,2), make_tuple(8,1294,1), make_tuple(81,1304,2), make_tuple(39,1306,1), make_tuple(87,1309,2), make_tuple(34,130,1), make_tuple(49,1310,1), make_tuple(15,1315,1), make_tuple(81,1327,2), make_tuple(38,1330,1), make_tuple(59,1330,2), make_tuple(8,1331,1), make_tuple(96,1332,2), make_tuple(83,133,2), make_tuple(12,1342,1), make_tuple(34,1342,1), make_tuple(64,1344,2), make_tuple(64,1348,2), make_tuple(24,1359,1), make_tuple(92,1363,2), make_tuple(14,137,1), make_tuple(66,137,2), make_tuple(7,1398,1), make_tuple(20,1398,1), make_tuple(71,1402,2), make_tuple(94,1438,2), make_tuple(5,1448,1), make_tuple(40,1451,1), make_tuple(81,1454,2), make_tuple(36,1467,1), make_tuple(100,1029,2), make_tuple(19,1032,1), make_tuple(27,1469,1), make_tuple(39,1476,1), make_tuple(90,1477,2), make_tuple(45,1487,1), make_tuple(20,1527,1), make_tuple(7,1529,1), make_tuple(26,1558,1), make_tuple(11,1561,1), make_tuple(53,1567,2), make_tuple(29,1589,1), make_tuple(88,1606,2), make_tuple(72,1625,2), make_tuple(5,1635,1), make_tuple(42,1646,1), make_tuple(10,1656,1), make_tuple(61,1669,2), make_tuple(10,1673,1), make_tuple(81,1682,2), make_tuple(26,1688,1), make_tuple(5,1696,1), make_tuple(39,169,1), make_tuple(57,171,2), make_tuple(56,1742,2), make_tuple(23,1747,1), make_tuple(45,1778,1), make_tuple(97,1782,2), make_tuple(45,1784,1), make_tuple(4,1789,1), make_tuple(35,17,1), make_tuple(85,1812,2), make_tuple(51,181,2), make_tuple(59,1839,2), make_tuple(75,1861,2), make_tuple(54,1862,2), make_tuple(91,1862,2), make_tuple(19,1867,1), make_tuple(77,1867,2), make_tuple(9,1870,1), make_tuple(80,1871,2), make_tuple(19,18,1), make_tuple(69,1903,2), make_tuple(87,1903,2), make_tuple(97,191,2), make_tuple(51,1925,2), make_tuple(5,1931,1), make_tuple(99,1935,2), make_tuple(4,1962,1), make_tuple(81,1962,2), make_tuple(88,1968,2), make_tuple(97,1969,2), make_tuple(52,198,2), make_tuple(24,2000,1), make_tuple(31,207,1), make_tuple(72,106,2), make_tuple(77,1074,2), make_tuple(46,1078,1), make_tuple(12,214,1), make_tuple(55,215,2), make_tuple(79,216,2), make_tuple(51,233,2), make_tuple(49,234,1), make_tuple(37,247,1), make_tuple(37,251,1), make_tuple(24,270,1), make_tuple(95,271,2), make_tuple(89,278,2), make_tuple(76,279,2), make_tuple(84,287,2), make_tuple(73,299,2), make_tuple(62,304,2), make_tuple(92,304,2), make_tuple(88,305,2), make_tuple(63,306,2), make_tuple(21,325,1), make_tuple(96,327,2), make_tuple(49,334,1), make_tuple(63,334,2), make_tuple(67,336,2), make_tuple(96,346,2), make_tuple(33,1082,1), make_tuple(32,348,1), make_tuple(76,349,2), make_tuple(40,375,1), make_tuple(57,375,2), make_tuple(19,383,1), make_tuple(47,3,1), make_tuple(13,409,1), make_tuple(51,427,2), make_tuple(64,427,2), make_tuple(73,428,2), make_tuple(35,442,1), make_tuple(23,475,1), make_tuple(30,47,1), make_tuple(24,1101,1), make_tuple(40,1104,1), make_tuple(11,4,1), make_tuple(23,508,1), make_tuple(71,511,2), make_tuple(73,517,2), make_tuple(51,525,2), make_tuple(66,530,2), make_tuple(14,538,1), make_tuple(53,549,2), make_tuple(94,558,2), make_tuple(94,563,2), make_tuple(64,567,2), make_tuple(79,572,2), make_tuple(22,57,1), make_tuple(88,583,2), make_tuple(13,586,1), make_tuple(43,590,1), make_tuple(14,592,1), make_tuple(32,603,1), make_tuple(74,607,2), make_tuple(9,624,1), make_tuple(39,629,1), make_tuple(87,630,2), make_tuple(6,635,1), make_tuple(61,63,2), make_tuple(97,1116,2), make_tuple(6,646,1), make_tuple(11,648,1), make_tuple(66,65,2), make_tuple(90,665,2), make_tuple(93,671,2), make_tuple(17,678,1), make_tuple(27,695,1), make_tuple(53,70,2), make_tuple(98,710,2), make_tuple(93,733,2), make_tuple(63,76,2), make_tuple(11,770,1), make_tuple(24,772,1), make_tuple(64,772,2), make_tuple(99,776,2), make_tuple(91,106,2), make_tuple(64,108,2), make_tuple(25,10,1), make_tuple(60,221,2), make_tuple(31,229,1), make_tuple(64,232,2), make_tuple(76,253,2), make_tuple(85,253,2), make_tuple(19,261,1), make_tuple(57,263,2), make_tuple(37,285,1), make_tuple(44,287,1), make_tuple(15,288,1), make_tuple(37,297,1), make_tuple(86,300,2), make_tuple(8,302,1), make_tuple(83,302,2), make_tuple(84,306,2), make_tuple(89,328,2), make_tuple(8,329,1), make_tuple(72,341,2), make_tuple(75,346,2), make_tuple(14,349,1), make_tuple(3,350,1), make_tuple(12,115,1), make_tuple(57,118,2), make_tuple(9,11,1), make_tuple(75,121,2), make_tuple(35,122,1), make_tuple(47,122,1), make_tuple(12,126,1), make_tuple(50,35,1), make_tuple(91,365,2), make_tuple(19,373,1), make_tuple(93,380,2), make_tuple(42,396,1), make_tuple(10,399,1), make_tuple(48,410,1), make_tuple(52,414,2), make_tuple(50,430,1), make_tuple(14,432,1), make_tuple(43,440,1), make_tuple(62,449,2), make_tuple(15,454,1), make_tuple(3,463,1), make_tuple(79,471,2), make_tuple(90,477,2), make_tuple(59,133,2), make_tuple(10,135,1), make_tuple(76,509,2), make_tuple(74,519,2), make_tuple(28,51,1), make_tuple(39,533,1), make_tuple(80,541,2), make_tuple(98,557,2), make_tuple(66,567,2), make_tuple(13,570,1), make_tuple(44,570,1), make_tuple(22,572,1), make_tuple(63,579,2), make_tuple(3,586,1), make_tuple(56,591,2), make_tuple(17,592,1), make_tuple(43,601,1), make_tuple(94,603,2), make_tuple(44,605,1), make_tuple(6,618,1), make_tuple(65,624,2), make_tuple(28,140,1), make_tuple(64,141,2), make_tuple(57,145,2), make_tuple(64,628,2), make_tuple(47,630,1), make_tuple(27,638,1), make_tuple(85,652,2), make_tuple(95,662,2), make_tuple(73,663,2), make_tuple(3,671,1), make_tuple(16,681,1), make_tuple(26,691,1), make_tuple(32,692,1), make_tuple(100,692,2), make_tuple(98,69,2), make_tuple(53,701,2), make_tuple(33,70,1), make_tuple(94,731,2), make_tuple(48,74,1), make_tuple(83,15,2), make_tuple(75,161,2), make_tuple(28,769,1), make_tuple(33,772,1), make_tuple(61,774,2), make_tuple(68,794,2), make_tuple(11,810,1), make_tuple(2,81,1), make_tuple(5,823,1), make_tuple(38,823,1), make_tuple(83,827,2), make_tuple(83,828,2), make_tuple(74,851,2), make_tuple(23,87,1), make_tuple(71,888,2), make_tuple(71,889,2), make_tuple(17,174,1), make_tuple(66,176,2), make_tuple(35,91,1), make_tuple(37,921,1), make_tuple(84,921,2), make_tuple(28,925,1), make_tuple(64,92,2), make_tuple(49,938,1), make_tuple(65,945,2), make_tuple(75,947,2), make_tuple(77,947,2), make_tuple(100,94,2), make_tuple(18,974,1), make_tuple(99,97,2), make_tuple(10,986,1), make_tuple(18,1,1), make_tuple(73,209,2), make_tuple(39,212,1)};

//set<long> inspection_set = {99795815,998615983,998751102,999957050};

//set<tuple<long,long,long>> inspection_copy = inspection_set_erl;

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

vector<short> g_lep_inds; //Holds poition of "analysis leptons"
short g_nlep;


const double Z_MASS = 91.1876;
const double W_MASS = 80.385;

const double BJET_CSV_TIGHT = 0.9535;
const double BJET_CSV_MED = 0.8484;
const double BJET_CSV_LOOSE = 0.5426;

double MAX_DR_JET_LEP_OVERLAP, JET_PT_MIN, JET_ETA_MAX, BJET_PT_MIN, BJET_ETA_MAX;
double W_JET_WINDOW_LOW, W_JET_WINDOW_HIGH, Z_VETO_WINDOW_LOW, Z_VETO_WINDOW_HIGH;

/* OLD Standards:
evt_type == flavor type. 0 == EE, 1 == MuMu, 2 == EMu
hyp_type == charge type. 0 == SS, 1 == OS, 2 == Photon
*/

enum flavor_type {EE = 0, EMu, MuMu};
const TString flavor_type_str[] = {"EE", "EMu", "MuMu"};
enum charge_type {SS = 0, OS, SFOS0, SFOS1, SFOS2, NA};
const TString charge_type_str[] = {"SS", "OS", "0SFOS", "1SFOS", "2SFOS", "NA"};

//Categories for Fake Rate Study
// The naming is the number of tight objects, followed by whether they are all real or whether there was some number of fakes.
// If none of the leps were tight, then there are the L_real and L_fake categories for when everything was loose. These are not
// broken down into how many are real and fake at the moment, just whether there were any reals or fakes.
enum FR_cat {T_real = 0, T_fake = 1, TT_real = 2, TT_fake = 3, TTT_real = 4, TTT_fake = 5, L_real = 6, L_fake = 7, CATERR = 8}; 
const TString FR_cats_str[] = {"T_real", "T_fake", "TT_real", "TT_fake", "TTT_real", "TTT_fake", "L_real", "L_fake", "CATERR"};

flavor_type FT;
charge_type CT;

//=============================
// Variable Computation
//=============================

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

/* Loops through lepton collection and finds the lep that has the max MT when combined with the MET */
double getMaxMTLepMET(bool verbose=false);

/* Builds the delta R (sqrt(dPhi^2 + dEta^2)) between the lepton at index id and the leading photon*/
double getdRGammaLep(short id=0);

/*Searches the vector of lorentz vectors for the pair with mass nearest (furthest) from the target mass if 'close' is true (false)*/
pair<int, int> getPairWithMass(const vector<LorentzVector> &vecs, double target_mass, bool close);

/*These get the 4 vectors for objects with mass closest or furthest from the mass of the W or Z using getPairWithMass()*/
pair<int,int> getMostZlikePair(const vector<LorentzVector> &vecs);
pair<int,int> getLeastZlikePair(const vector<LorentzVector> &vecs);
pair<int,int> getMostWlikePair(const vector<LorentzVector> &vecs);

/* Goes pairwise through jet collection and finds the pair of jets with the minium delta eta. */
pair<int,int> getClosestJetsInEta();

/* Goes pairwise through jet collection and finds the pair of jets with the minium delta R. */
pair<int,int> getClosestJetsInDR();

/*Loops through pairs of entries in the lep_pdgId vector and counts how many have opposite value*/
int getNumOSSFPairs();

/*Returns a string which gives the exact flavor composition of the leptons in the event*/
TString getLepFlavorString();

/*Returns the DeltaR between objects p1 and p2.*/
double DeltaR(const LorentzVector p1, const LorentzVector p2);

/* checks the gen record for a lepton with the index specified */
bool isCleanLepFromW(int index);

/* Checks whether the two leading analysis leptons are EE, EMu, or MuMu */
flavor_type getFlavorType();

/* Returns the charge type, SS or OS for 2 lep events, and SFOS0/1/2 for more than 2leps */
charge_type getChargeType();

/* Returns the category for the event in the fake rate study. Loops through leps, counts reals and how many pass the tight ID, returns appropriate category.*/
FR_cat getFRCategory();
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
/* Loops through lepton objects and adds indexed to g_lep_inds if they pass the tight selection (or fakable object when we are doing fake rate study)*/
void setLepIndexes();

/*Takes in a p4 for a jet and determines whether that jet is less than MAX_DR_JET_LEP_OVERLAP for all leptons*/
bool isOverlapJet(const LorentzVector &jet_p4);

/*This function writes only elements of the given vector to g_jets_p4 variable, it can be passed the up and down variations as well.*/
void writeCleanedJets(const vector<LorentzVector> &vecs);

/*This function writes only elements of the given vector to g_jets_medb_p4 variable, it can be passed the up and down variations as well. Distingushed from the regular jet function because we need to keep track of the CSV values as well.*/
void writeCleanedBJets(const vector<LorentzVector> &vecs, const vector<float> &csvs);

/*Sets up global variables for the event which are the quantities that might be fluctuated in the process of computing uncertainty limits*/
void setupGlobals();

/* Sets up constant valued variables */
void setupConstants();

/*Obvi the event looper*/
int ScanChain( TChain* chain, ConfigParser *configuration, bool fast = true, int nEvents = -1);