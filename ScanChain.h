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

// Needed for global vars for branches
//#include <functional>


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
int num_events_veto_WH=0;
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

bool FRS, LooseIso, FRS_use_veto;
//std::function<const std::vector<bool>&()> g_looseIDs; //holds loose IDs when doing the Fake Rate study.

// ----------------
// DEBUG MODE
// ----------------

//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(1067496,1,5339), make_tuple(1783520,1,8920), make_tuple(1591339,1,7959), make_tuple(362878,1,1815), make_tuple(805669,1,4030), make_tuple(1581054,1,7907), make_tuple(810619,1,4054), make_tuple(1968205,1,9844), make_tuple(1989888,1,9952), make_tuple(611008,1,3056), make_tuple(825488,1,4129), make_tuple(490909,1,2456), make_tuple(818639,1,4095), make_tuple(239158,1,1197), make_tuple(180411,1,903), make_tuple(629362,1,3148), make_tuple(808257,1,4043), make_tuple(108840,1,545), make_tuple(167012,1,836), make_tuple(262537,1,1313), make_tuple(262577,1,1314), make_tuple(1502941,1,7517), make_tuple(146118,1,731), make_tuple(964878,1,4826), make_tuple(589212,1,2947), make_tuple(607680,1,3039), make_tuple(776626,1,3884), make_tuple(28345,1,142), make_tuple(1887676,1,9441), make_tuple(1241420,1,6209), make_tuple(15478,1,78), make_tuple(113678,1,569), make_tuple(1728847,1,8647), make_tuple(19214,1,97), make_tuple(122725,1,614), make_tuple(969067,1,4847), make_tuple(388614,1,1944), make_tuple(1152260,1,5763), make_tuple(1453273,1,7268), make_tuple(1310096,1,6552), make_tuple(199591,1,999), make_tuple(1796630,1,8986), make_tuple(1760240,1,8804), make_tuple(492339,1,2463), make_tuple(1136073,1,5682), make_tuple(1253408,1,6269), make_tuple(1260013,1,6302), make_tuple(1864340,1,9324), make_tuple(1882534,1,9415), make_tuple(1718703,1,8596), make_tuple(1346852,1,6736), make_tuple(1045410,1,5229), make_tuple(1543892,1,7722), make_tuple(1703944,1,8522), make_tuple(161722,1,809), make_tuple(246153,1,1232), make_tuple(1108952,1,5546), make_tuple(1109445,1,5549), make_tuple(1116916,1,5586), make_tuple(1548746,1,7746), make_tuple(1375939,1,6882), make_tuple(344370,1,1723), make_tuple(1353835,1,6771), make_tuple(572944,1,2866), make_tuple(993485,1,4969), make_tuple(1475107,1,7378), make_tuple(143,1,1), make_tuple(489071,1,2446), make_tuple(149200,1,747), make_tuple(407201,1,2037), make_tuple(738443,1,3693), make_tuple(1528437,1,7644), make_tuple(62642,1,314), make_tuple(507786,1,2540), make_tuple(267049,1,1336), make_tuple(1017142,1,5087), make_tuple(1084531,1,5424), make_tuple(148197,1,742), make_tuple(578874,1,2895), make_tuple(835177,1,4177), make_tuple(493514,1,2469), make_tuple(1723647,1,8621), make_tuple(52319,1,262), make_tuple(415674,1,2079), make_tuple(594900,1,2976), make_tuple(624818,1,3125), make_tuple(515269,1,2577), make_tuple(1921546,1,9610), make_tuple(1993042,1,9968), make_tuple(1133080,1,5667), make_tuple(1486713,1,7436), make_tuple(1387628,1,6940), make_tuple(1118959,1,5596), make_tuple(383749,1,1920), make_tuple(667370,1,3338), make_tuple(182808,1,915), make_tuple(1417143,1,7088), make_tuple(1924373,1,9624)};

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
vector<bool> g_looseIDs;
vector<bool> g_tightIDs;


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

//Categories for different background samples
enum BG_cat {trueSS = 0, chargeFlip = 1, LLSS = 2, photonFake = 3, fake = 4, photonDoubleFake = 5, fake_photonFake = 6,  doubleFake = 7, otherPhotonFake = 8, trueWWW = 9, LL3l = 10, true3l = 11, photonTripleFake = 12, other = 13}; 
const TString BG_cats_str[] = {"trueSS", "chargeFlip", "LLSS", "photonFake", "fake", "photonDoubleFake", "fake_photonFake",  "doubleFake", "otherPhotonFake", "trueWWW", "LL3l", "true3l", "photonTripleFake", "other"};

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

/* Builds the MT from the lepton at index id and the MET vector (assumes massless particles) */
double getMTLepMET(short id=0);

/* Loops through lepton collection and finds the lep that has the max MT when combined with the MET */
double getMaxMTLepMET(bool verbose=false);

/* Builds the delta R (sqrt(dPhi^2 + dEta^2)) between the lepton at index id and the leading photon */
double getdRGammaLep(short id=0);

/* Returns the Cone Corrected pt for the FR Closure study */
double getFRConeCorrPt(int lep_index);

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

/* Returns background category for the event */
BG_cat getBGCategory();

/* Checks whether the higgs decayed to WW in VH MC */
bool passGenLevelWHWWW();
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