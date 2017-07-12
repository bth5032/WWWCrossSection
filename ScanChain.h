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

//HJ Debug
//---------
// EE
//---------
//tt1l
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(172521043,1,209625), make_tuple(34895595,1,42401), make_tuple(76708617,1,93207), make_tuple(163800554,1,195895), make_tuple(52822818,1,63173), make_tuple(1765116,1,2111), make_tuple(148359468,1,177429)};
//WZ
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(41059,1,206), make_tuple(179208,1,897), make_tuple(694076,1,3472), make_tuple(628873,1,3146), make_tuple(1984834,1,9927), make_tuple(1877546,1,9390), make_tuple(583273,1,2917), make_tuple(1715644,1,8581), make_tuple(249115,1,1246), make_tuple(1083013,1,5417)};
//Wjets
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(231344423,1,120807), make_tuple(38348441,1,21876), make_tuple(48896986,1,21854)};
//WWW
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(13,260,1), make_tuple(40,320,1), make_tuple(87,320,2), make_tuple(24,324,1), make_tuple(11,32,1), make_tuple(52,349,2), make_tuple(94,34,2), make_tuple(51,350,2), make_tuple(43,351,1), make_tuple(69,352,2), make_tuple(88,356,2), make_tuple(8,395,1), make_tuple(63,408,2), make_tuple(45,48,1), make_tuple(78,490,2), make_tuple(30,136,1), make_tuple(70,499,2), make_tuple(98,508,2), make_tuple(26,604,1), make_tuple(32,61,1), make_tuple(90,694,2), make_tuple(45,710,1), make_tuple(90,765,2), make_tuple(86,783,2), make_tuple(32,797,1), make_tuple(91,850,2), make_tuple(15,912,1), make_tuple(25,985,1), make_tuple(44,9,1), make_tuple(20,191,1), make_tuple(37,204,1), make_tuple(96,1003,2), make_tuple(15,796,1), make_tuple(3,805,1), make_tuple(44,84,1), make_tuple(22,972,1), make_tuple(88,1181,2), make_tuple(25,1213,1), make_tuple(97,1231,2), make_tuple(39,1250,1), make_tuple(36,1276,1), make_tuple(67,1388,2), make_tuple(19,1459,1), make_tuple(81,1036,2), make_tuple(73,1037,2), make_tuple(82,1525,2), make_tuple(86,1565,2), make_tuple(87,1678,2), make_tuple(100,1706,2), make_tuple(21,1710,1), make_tuple(83,1771,2), make_tuple(4,187,1), make_tuple(98,1950,2), make_tuple(17,235,1), make_tuple(58,256,2), make_tuple(36,295,1), make_tuple(80,36,2), make_tuple(75,382,2), make_tuple(80,397,2), make_tuple(17,417,1), make_tuple(95,445,2), make_tuple(57,507,2), make_tuple(20,529,1), make_tuple(24,552,1), make_tuple(73,580,2), make_tuple(94,590,2), make_tuple(54,591,2), make_tuple(52,620,2), make_tuple(62,656,2), make_tuple(20,667,1), make_tuple(28,742,1), make_tuple(65,759,2), make_tuple(36,775,1)};


//---------
// MM
//---------
//tt1l
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(115004552,1,139739), make_tuple(201956612,1,245391), make_tuple(200525870,1,243653), make_tuple(137699359,1,167314), make_tuple(85452148,1,103831), make_tuple(67497145,1,82014), make_tuple(194120538,1,232156), make_tuple(212303504,1,253901), make_tuple(182285055,1,218001), make_tuple(296084,1,355), make_tuple(108086802,1,129265), make_tuple(194380305,1,232466)};
//WZ
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(512746,1,2565), make_tuple(1499623,1,7500), make_tuple(258044,1,1291), make_tuple(285851,1,1430), make_tuple(374101,1,1871), make_tuple(743322,1,3718), make_tuple(608167,1,3042), make_tuple(596703,1,2985), make_tuple(988634,1,4945), make_tuple(917017,1,4587), make_tuple(774211,1,3872), make_tuple(884299,1,4423), make_tuple(177997,1,891), make_tuple(1758015,1,8792), make_tuple(660673,1,3305), make_tuple(4986,1,25), make_tuple(275595,1,1379), make_tuple(1415130,1,7078), make_tuple(1799119,1,8998), make_tuple(1973277,1,9869), make_tuple(1685877,1,8432), make_tuple(1633886,1,8172), make_tuple(1715381,1,8579), make_tuple(601792,1,3010), make_tuple(839402,1,4198), make_tuple(673268,1,3368), make_tuple(1307199,1,6538), make_tuple(1856925,1,9287), make_tuple(537646,1,2689), make_tuple(303177,1,1517), make_tuple(1651741,1,8261), make_tuple(1655157,1,8278), make_tuple(1878270,1,9394), make_tuple(1878497,1,9395), make_tuple(42556,1,213), make_tuple(1924691,1,9626), make_tuple(1541824,1,7711)};
//Wjets
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(76073390,1,34000)};
//WWW
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(91,106,2), make_tuple(64,108,2), make_tuple(25,10,1), make_tuple(60,221,2), make_tuple(31,229,1), make_tuple(64,232,2), make_tuple(76,253,2), make_tuple(85,253,2), make_tuple(19,261,1), make_tuple(57,263,2), make_tuple(37,285,1), make_tuple(44,287,1), make_tuple(15,288,1), make_tuple(37,297,1), make_tuple(86,300,2), make_tuple(8,302,1), make_tuple(83,302,2), make_tuple(84,306,2), make_tuple(89,328,2), make_tuple(8,329,1), make_tuple(72,341,2), make_tuple(75,346,2), make_tuple(14,349,1), make_tuple(3,350,1), make_tuple(12,115,1), make_tuple(57,118,2), make_tuple(9,11,1), make_tuple(75,121,2), make_tuple(35,122,1), make_tuple(47,122,1), make_tuple(12,126,1), make_tuple(50,35,1), make_tuple(91,365,2), make_tuple(19,373,1), make_tuple(93,380,2), make_tuple(42,396,1), make_tuple(10,399,1), make_tuple(48,410,1), make_tuple(52,414,2), make_tuple(50,430,1), make_tuple(14,432,1), make_tuple(43,440,1), make_tuple(62,449,2), make_tuple(15,454,1), make_tuple(3,463,1), make_tuple(79,471,2), make_tuple(90,477,2), make_tuple(59,133,2), make_tuple(10,135,1), make_tuple(76,509,2), make_tuple(74,519,2), make_tuple(28,51,1), make_tuple(39,533,1), make_tuple(80,541,2), make_tuple(98,557,2), make_tuple(66,567,2), make_tuple(13,570,1), make_tuple(44,570,1), make_tuple(22,572,1), make_tuple(63,579,2), make_tuple(3,586,1), make_tuple(56,591,2), make_tuple(17,592,1), make_tuple(43,601,1), make_tuple(94,603,2), make_tuple(44,605,1), make_tuple(6,618,1), make_tuple(65,624,2), make_tuple(28,140,1), make_tuple(64,141,2), make_tuple(57,145,2), make_tuple(64,628,2), make_tuple(47,630,1), make_tuple(27,638,1), make_tuple(85,652,2), make_tuple(95,662,2), make_tuple(73,663,2), make_tuple(3,671,1), make_tuple(16,681,1), make_tuple(26,691,1), make_tuple(32,692,1), make_tuple(100,692,2), make_tuple(98,69,2), make_tuple(53,701,2), make_tuple(33,70,1), make_tuple(94,731,2), make_tuple(48,74,1), make_tuple(83,15,2), make_tuple(75,161,2), make_tuple(28,769,1), make_tuple(33,772,1), make_tuple(61,774,2), make_tuple(68,794,2), make_tuple(11,810,1), make_tuple(2,81,1), make_tuple(5,823,1), make_tuple(38,823,1), make_tuple(83,827,2), make_tuple(83,828,2), make_tuple(74,851,2), make_tuple(23,87,1), make_tuple(71,888,2), make_tuple(71,889,2), make_tuple(17,174,1), make_tuple(66,176,2), make_tuple(35,91,1), make_tuple(37,921,1), make_tuple(84,921,2), make_tuple(28,925,1), make_tuple(64,92,2), make_tuple(49,938,1), make_tuple(65,945,2), make_tuple(75,947,2), make_tuple(77,947,2), make_tuple(100,94,2), make_tuple(18,974,1), make_tuple(99,97,2), make_tuple(10,986,1), make_tuple(18,1,1), make_tuple(73,209,2), make_tuple(39,212,1), make_tuple(40,1005,1), make_tuple(19,112,1), make_tuple(20,1133,1), make_tuple(66,788,2), make_tuple(97,791,2), make_tuple(18,797,1), make_tuple(42,812,1), make_tuple(100,812,2), make_tuple(60,815,2), make_tuple(22,825,1), make_tuple(35,842,1), make_tuple(6,847,1), make_tuple(59,84,2), make_tuple(5,852,1), make_tuple(9,860,1), make_tuple(59,891,2), make_tuple(15,899,1), make_tuple(3,906,1), make_tuple(28,113,1), make_tuple(96,1147,2), make_tuple(28,1148,1), make_tuple(63,1148,2), make_tuple(36,925,1), make_tuple(33,936,1), make_tuple(12,948,1), make_tuple(48,956,1), make_tuple(68,968,2), make_tuple(55,970,2), make_tuple(9,977,1), make_tuple(15,977,1), make_tuple(67,991,2), make_tuple(72,995,2), make_tuple(64,1154,2), make_tuple(81,1156,2), make_tuple(6,1160,1), make_tuple(78,1179,2), make_tuple(39,1186,1), make_tuple(82,1188,2), make_tuple(78,1216,2), make_tuple(12,1223,1), make_tuple(16,1237,1), make_tuple(59,1261,2), make_tuple(8,1262,1), make_tuple(70,1016,2), make_tuple(5,1020,1), make_tuple(69,1021,2), make_tuple(4,1289,1), make_tuple(95,1289,2), make_tuple(8,1294,1), make_tuple(81,1304,2), make_tuple(39,1306,1), make_tuple(87,1309,2), make_tuple(34,130,1), make_tuple(49,1310,1), make_tuple(15,1315,1), make_tuple(81,1327,2), make_tuple(38,1330,1), make_tuple(59,1330,2), make_tuple(8,1331,1), make_tuple(96,1332,2), make_tuple(83,133,2), make_tuple(12,1342,1), make_tuple(34,1342,1), make_tuple(64,1344,2), make_tuple(64,1348,2), make_tuple(24,1359,1), make_tuple(92,1363,2), make_tuple(14,137,1), make_tuple(66,137,2), make_tuple(7,1398,1), make_tuple(20,1398,1), make_tuple(71,1402,2), make_tuple(94,1438,2), make_tuple(5,1448,1), make_tuple(40,1451,1), make_tuple(81,1454,2), make_tuple(36,1467,1), make_tuple(100,1029,2), make_tuple(19,1032,1), make_tuple(27,1469,1), make_tuple(39,1476,1), make_tuple(90,1477,2), make_tuple(45,1487,1), make_tuple(20,1527,1), make_tuple(7,1529,1), make_tuple(26,1558,1), make_tuple(11,1561,1), make_tuple(53,1567,2), make_tuple(29,1589,1), make_tuple(88,1606,2), make_tuple(72,1625,2), make_tuple(5,1635,1), make_tuple(42,1646,1), make_tuple(10,1656,1), make_tuple(61,1669,2), make_tuple(10,1673,1), make_tuple(81,1682,2), make_tuple(26,1688,1), make_tuple(5,1696,1), make_tuple(39,169,1), make_tuple(57,171,2), make_tuple(56,1742,2), make_tuple(23,1747,1), make_tuple(45,1778,1), make_tuple(97,1782,2), make_tuple(45,1784,1), make_tuple(4,1789,1), make_tuple(35,17,1), make_tuple(85,1812,2), make_tuple(51,181,2), make_tuple(59,1839,2), make_tuple(75,1861,2), make_tuple(54,1862,2), make_tuple(91,1862,2), make_tuple(19,1867,1), make_tuple(77,1867,2), make_tuple(9,1870,1), make_tuple(80,1871,2), make_tuple(19,18,1), make_tuple(69,1903,2), make_tuple(87,1903,2), make_tuple(97,191,2), make_tuple(51,1925,2), make_tuple(5,1931,1), make_tuple(99,1935,2), make_tuple(4,1962,1), make_tuple(81,1962,2), make_tuple(88,1968,2), make_tuple(97,1969,2), make_tuple(52,198,2), make_tuple(24,2000,1), make_tuple(31,207,1), make_tuple(72,106,2), make_tuple(77,1074,2), make_tuple(46,1078,1), make_tuple(12,214,1), make_tuple(55,215,2), make_tuple(79,216,2), make_tuple(51,233,2), make_tuple(49,234,1), make_tuple(37,247,1), make_tuple(37,251,1), make_tuple(24,270,1), make_tuple(95,271,2), make_tuple(89,278,2), make_tuple(76,279,2), make_tuple(84,287,2), make_tuple(73,299,2), make_tuple(62,304,2), make_tuple(92,304,2), make_tuple(88,305,2), make_tuple(63,306,2), make_tuple(21,325,1), make_tuple(96,327,2), make_tuple(49,334,1), make_tuple(63,334,2), make_tuple(67,336,2), make_tuple(96,346,2), make_tuple(33,1082,1), make_tuple(32,348,1), make_tuple(76,349,2), make_tuple(40,375,1), make_tuple(57,375,2), make_tuple(19,383,1), make_tuple(47,3,1), make_tuple(13,409,1), make_tuple(51,427,2), make_tuple(64,427,2), make_tuple(73,428,2), make_tuple(35,442,1), make_tuple(23,475,1), make_tuple(30,47,1), make_tuple(24,1101,1), make_tuple(40,1104,1), make_tuple(11,4,1), make_tuple(23,508,1), make_tuple(71,511,2), make_tuple(73,517,2), make_tuple(51,525,2), make_tuple(66,530,2), make_tuple(14,538,1), make_tuple(53,549,2), make_tuple(94,558,2), make_tuple(94,563,2), make_tuple(64,567,2), make_tuple(79,572,2), make_tuple(22,57,1), make_tuple(88,583,2), make_tuple(13,586,1), make_tuple(43,590,1), make_tuple(14,592,1), make_tuple(32,603,1), make_tuple(74,607,2), make_tuple(9,624,1), make_tuple(39,629,1), make_tuple(87,630,2), make_tuple(6,635,1), make_tuple(61,63,2), make_tuple(97,1116,2), make_tuple(6,646,1), make_tuple(11,648,1), make_tuple(66,65,2), make_tuple(90,665,2), make_tuple(93,671,2), make_tuple(17,678,1), make_tuple(27,695,1), make_tuple(53,70,2), make_tuple(98,710,2), make_tuple(93,733,2), make_tuple(63,76,2), make_tuple(11,770,1), make_tuple(24,772,1), make_tuple(64,772,2), make_tuple(99,776,2)};
//---------
// EM
//---------
//tt1l
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(88453374,1,107477), make_tuple(88716550,1,107797), make_tuple(160631391,1,195178), make_tuple(44908560,1,53708), make_tuple(170922018,1,204412), make_tuple(77128802,1,92241), make_tuple(151493647,1,181177), make_tuple(37771506,1,45173), make_tuple(17018342,1,20353), make_tuple(65100836,1,77857), make_tuple(135279842,1,161786), make_tuple(133470329,1,159622)};
//WZ
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(230491,1,1153), make_tuple(60717,1,304), make_tuple(1199407,1,5999), make_tuple(898282,1,4493), make_tuple(268223,1,1342), make_tuple(1827384,1,9139), make_tuple(133699,1,669), make_tuple(689933,1,3451), make_tuple(246981,1,1236), make_tuple(92704,1,464), make_tuple(1107528,1,5539), make_tuple(1359736,1,6801), make_tuple(1882657,1,9416), make_tuple(997493,1,4989), make_tuple(1673228,1,8368), make_tuple(1544227,1,7723), make_tuple(447705,1,2239), make_tuple(1603355,1,8019), make_tuple(557903,1,2791), make_tuple(795099,1,3977), make_tuple(909674,1,4550), make_tuple(585137,1,2927), make_tuple(736755,1,3685), make_tuple(1086930,1,5436), make_tuple(1087555,1,5439), make_tuple(1411240,1,7058), make_tuple(576006,1,2881), make_tuple(1624872,1,8127), make_tuple(1315667,1,6580)};
//Wjets
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(58063318,1,33123), make_tuple(9408033,1,4205), make_tuple(152914932,1,68342)};
//WWW
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(95,110,2), make_tuple(4,222,1), make_tuple(90,223,2), make_tuple(3,231,1), make_tuple(98,24,2), make_tuple(95,254,2), make_tuple(55,258,2), make_tuple(55,264,2), make_tuple(72,268,2), make_tuple(18,313,1), make_tuple(67,31,2), make_tuple(75,337,2), make_tuple(57,341,2), make_tuple(32,353,1), make_tuple(37,116,1), make_tuple(95,119,2), make_tuple(48,365,1), make_tuple(20,367,1), make_tuple(22,36,1), make_tuple(38,384,1), make_tuple(63,397,2), make_tuple(91,402,2), make_tuple(70,404,2), make_tuple(30,418,1), make_tuple(99,427,2), make_tuple(58,434,2), make_tuple(11,455,1), make_tuple(39,469,1), make_tuple(18,483,1), make_tuple(84,134,2), make_tuple(90,13,2), make_tuple(76,505,2), make_tuple(21,528,1), make_tuple(83,531,2), make_tuple(82,538,2), make_tuple(59,545,2), make_tuple(79,545,2), make_tuple(44,553,1), make_tuple(70,554,2), make_tuple(49,572,1), make_tuple(28,576,1), make_tuple(74,582,2), make_tuple(65,584,2), make_tuple(25,586,1), make_tuple(42,590,1), make_tuple(68,619,2), make_tuple(45,622,1), make_tuple(7,150,1), make_tuple(66,153,2), make_tuple(97,629,2), make_tuple(24,630,1), make_tuple(19,645,1), make_tuple(73,646,2), make_tuple(12,658,1), make_tuple(15,659,1), make_tuple(73,682,2), make_tuple(47,68,1), make_tuple(73,695,2), make_tuple(63,707,2), make_tuple(8,723,1), make_tuple(74,72,2), make_tuple(71,733,2), make_tuple(77,739,2), make_tuple(71,752,2), make_tuple(62,760,2), make_tuple(2,15,1), make_tuple(81,163,2), make_tuple(52,763,2), make_tuple(3,767,1), make_tuple(95,775,2), make_tuple(46,776,1), make_tuple(69,77,2), make_tuple(35,782,1), make_tuple(32,790,1), make_tuple(27,803,1), make_tuple(44,805,1), make_tuple(30,827,1), make_tuple(70,82,2), make_tuple(79,847,2), make_tuple(14,869,1), make_tuple(69,86,2), make_tuple(64,88,2), make_tuple(89,168,2), make_tuple(25,16,1), make_tuple(56,910,2), make_tuple(49,922,1), make_tuple(87,931,2), make_tuple(36,936,1), make_tuple(46,946,1), make_tuple(52,994,2), make_tuple(75,995,2), make_tuple(32,21,1), make_tuple(86,220,2), make_tuple(51,1000,2), make_tuple(15,1127,1), make_tuple(76,1127,2), make_tuple(78,1127,2), make_tuple(21,1133,1), make_tuple(48,1133,1), make_tuple(53,1134,2), make_tuple(27,781,1), make_tuple(32,785,1), make_tuple(51,788,2), make_tuple(89,802,2), make_tuple(46,803,1), make_tuple(36,805,1), make_tuple(93,825,2), make_tuple(70,835,2), make_tuple(100,839,2), make_tuple(67,841,2), make_tuple(56,869,2), make_tuple(68,873,2), make_tuple(65,891,2), make_tuple(53,896,2), make_tuple(98,906,2), make_tuple(40,910,1), make_tuple(72,92,2), make_tuple(96,92,2), make_tuple(97,931,2), make_tuple(71,941,2), make_tuple(71,991,2), make_tuple(54,992,2), make_tuple(64,1158,2), make_tuple(41,1167,1), make_tuple(58,1180,2), make_tuple(2,1204,1), make_tuple(98,1218,2), make_tuple(56,1254,2), make_tuple(4,1270,1), make_tuple(12,1279,1), make_tuple(47,1022,1), make_tuple(67,1284,2), make_tuple(96,1306,2), make_tuple(18,1309,1), make_tuple(8,1362,1), make_tuple(29,137,1), make_tuple(13,1384,1), make_tuple(47,1388,1), make_tuple(32,1390,1), make_tuple(38,1415,1), make_tuple(59,1439,2), make_tuple(20,1454,1), make_tuple(79,1034,2), make_tuple(48,1035,1), make_tuple(15,1500,1), make_tuple(82,1505,2), make_tuple(89,1514,2), make_tuple(64,1536,2), make_tuple(46,1541,1), make_tuple(42,154,1), make_tuple(25,155,1), make_tuple(16,1569,1), make_tuple(87,1571,2), make_tuple(96,1578,2), make_tuple(17,157,1), make_tuple(98,1623,2), make_tuple(96,1626,2), make_tuple(8,1635,1), make_tuple(82,1672,2), make_tuple(99,1682,2), make_tuple(36,168,1), make_tuple(2,1691,1), make_tuple(38,1694,1), make_tuple(13,1702,1), make_tuple(71,1705,2), make_tuple(37,1718,1), make_tuple(29,1744,1), make_tuple(10,176,1), make_tuple(22,1772,1), make_tuple(64,1794,2), make_tuple(97,1805,2), make_tuple(27,1829,1), make_tuple(73,1832,2), make_tuple(53,1835,2), make_tuple(54,1839,2), make_tuple(31,183,1), make_tuple(93,1842,2), make_tuple(83,184,2), make_tuple(15,1057,1), make_tuple(8,1063,1), make_tuple(3,1065,1), make_tuple(7,1067,1), make_tuple(12,1869,1), make_tuple(67,1875,2), make_tuple(89,187,2), make_tuple(3,1918,1), make_tuple(80,1919,2), make_tuple(4,1931,1), make_tuple(61,1932,2), make_tuple(6,193,1), make_tuple(4,1970,1), make_tuple(22,20,1), make_tuple(58,257,2), make_tuple(71,268,2), make_tuple(19,286,1), make_tuple(38,300,1), make_tuple(97,302,2), make_tuple(59,318,2), make_tuple(28,31,1), make_tuple(25,321,1), make_tuple(92,328,2), make_tuple(68,332,2), make_tuple(19,335,1), make_tuple(38,340,1), make_tuple(29,349,1), make_tuple(29,355,1), make_tuple(29,359,1), make_tuple(2,364,1), make_tuple(69,366,2), make_tuple(3,394,1), make_tuple(25,395,1), make_tuple(84,39,2), make_tuple(55,403,2), make_tuple(20,404,1), make_tuple(59,407,2), make_tuple(38,409,1), make_tuple(57,439,2), make_tuple(13,443,1), make_tuple(42,450,1), make_tuple(46,45,1), make_tuple(27,461,1), make_tuple(44,467,1), make_tuple(19,473,1), make_tuple(57,47,2), make_tuple(64,47,2), make_tuple(72,1101,2), make_tuple(26,493,1), make_tuple(79,506,2), make_tuple(88,531,2), make_tuple(43,534,1), make_tuple(15,542,1), make_tuple(30,566,1), make_tuple(56,566,2), make_tuple(89,574,2), make_tuple(84,580,2), make_tuple(82,581,2), make_tuple(68,604,2), make_tuple(5,609,1), make_tuple(79,619,2), make_tuple(94,623,2), make_tuple(95,628,2), make_tuple(35,62,1), make_tuple(58,641,2), make_tuple(33,1109,1), make_tuple(91,110,2), make_tuple(53,1117,2), make_tuple(15,653,1), make_tuple(85,654,2), make_tuple(72,677,2), make_tuple(32,683,1), make_tuple(60,686,2), make_tuple(80,688,2), make_tuple(88,69,2), make_tuple(81,700,2), make_tuple(5,710,1), make_tuple(54,718,2), make_tuple(26,71,1), make_tuple(100,721,2), make_tuple(95,722,2), make_tuple(36,732,1), make_tuple(74,743,2), make_tuple(88,743,2), make_tuple(48,744,1), make_tuple(9,753,1), make_tuple(38,769,1), make_tuple(66,772,2), make_tuple(33,775,1), make_tuple(63,775,2)};
//---------
// 0SFOS
//---------
//tt1l
//set<tuple<long,long,long>>  inspection_set_erl = {};
//WZ
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(872570,1,4364), make_tuple(935087,1,4677)};
//Wjets
//set<tuple<long,long,long>>  inspection_set_erl = {};
//WWW
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(21,107,1), make_tuple(18,257,1), make_tuple(37,259,1), make_tuple(89,339,2), make_tuple(34,389,1), make_tuple(42,405,1), make_tuple(58,451,2), make_tuple(58,130,2), make_tuple(18,499,1), make_tuple(55,669,2), make_tuple(39,715,1), make_tuple(86,732,2), make_tuple(89,74,2), make_tuple(33,771,1), make_tuple(61,786,2), make_tuple(42,836,1), make_tuple(30,842,1), make_tuple(89,846,2), make_tuple(74,885,2), make_tuple(62,906,2), make_tuple(75,969,2), make_tuple(42,975,1), make_tuple(39,994,1), make_tuple(92,191,2), make_tuple(28,192,1), make_tuple(23,1000,1), make_tuple(66,1122,2), make_tuple(4,816,1), make_tuple(25,864,1), make_tuple(18,909,1), make_tuple(46,1147,1), make_tuple(93,94,2), make_tuple(11,953,1), make_tuple(41,993,1), make_tuple(2,1203,1), make_tuple(50,1219,1), make_tuple(98,1222,2), make_tuple(58,1278,2), make_tuple(36,1342,1), make_tuple(20,1390,1), make_tuple(65,1398,2), make_tuple(48,1401,1), make_tuple(59,1411,2), make_tuple(54,1442,2), make_tuple(24,1032,1), make_tuple(40,1469,1), make_tuple(20,1515,1), make_tuple(10,1623,1), make_tuple(56,164,2), make_tuple(62,1674,2), make_tuple(55,1048,2), make_tuple(7,1051,1), make_tuple(66,1717,2), make_tuple(63,1759,2), make_tuple(18,1781,1), make_tuple(83,1790,2), make_tuple(5,1826,1), make_tuple(24,1872,1), make_tuple(37,1875,1), make_tuple(57,187,2), make_tuple(42,1911,1), make_tuple(76,191,2), make_tuple(96,1951,2), make_tuple(61,1987,2), make_tuple(64,1078,2), make_tuple(92,239,2), make_tuple(8,268,1), make_tuple(4,276,1), make_tuple(8,285,1), make_tuple(7,313,1), make_tuple(79,408,2), make_tuple(16,432,1), make_tuple(80,503,2), make_tuple(80,509,2), make_tuple(57,519,2), make_tuple(68,529,2), make_tuple(4,541,1), make_tuple(96,562,2), make_tuple(73,609,2), make_tuple(16,614,1), make_tuple(78,621,2), make_tuple(100,635,2), make_tuple(91,1114,2), make_tuple(43,716,1)};
//---------
// 1SFOS
//---------
//tt1l
//set<tuple<long,long,long>>  inspection_set_erl = {};
//WZ
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(2370770,1,6783), make_tuple(592822,1,2965), make_tuple(1443784,1,7221), make_tuple(1572480,1,7865), make_tuple(1970501,1,9855), make_tuple(364022,1,1821), make_tuple(703929,1,3521), make_tuple(205436,1,1028), make_tuple(1405506,1,7030), make_tuple(1779341,1,8899), make_tuple(1301949,1,6512), make_tuple(1007702,1,5040), make_tuple(1484283,1,7423), make_tuple(104990,1,526), make_tuple(588446,1,2943), make_tuple(1514407,1,7574), make_tuple(1571931,1,7862), make_tuple(330671,1,1654), make_tuple(885844,1,4431), make_tuple(380432,1,1903), make_tuple(792194,1,3962), make_tuple(1056221,1,5283), make_tuple(212328,1,1062), make_tuple(830942,1,4156), make_tuple(286436,1,1433), make_tuple(34194,1,171), make_tuple(189632,1,949), make_tuple(395638,1,1979), make_tuple(9530,1,48), make_tuple(276258,1,1382), make_tuple(89606,1,449), make_tuple(108617,1,544), make_tuple(108829,1,545), make_tuple(263564,1,1319), make_tuple(1110504,1,5554), make_tuple(1276620,1,6385), make_tuple(304936,1,1525), make_tuple(553738,1,2770), make_tuple(456940,1,2286), make_tuple(707557,1,3539), make_tuple(1130442,1,5654), make_tuple(453566,1,2269), make_tuple(642676,1,3215), make_tuple(1529809,1,7651), make_tuple(319565,1,1599), make_tuple(822503,1,4114), make_tuple(822622,1,4114), make_tuple(1578408,1,7894), make_tuple(1689555,1,8450), make_tuple(1248871,1,6246), make_tuple(664977,1,3326), make_tuple(138548,1,693), make_tuple(334819,1,1675), make_tuple(1412373,1,7064), make_tuple(1615967,1,8082), make_tuple(1151176,1,5758), make_tuple(1089805,1,5451), make_tuple(1107735,1,5540), make_tuple(542959,1,2716), make_tuple(1293568,1,6470), make_tuple(307869,1,1540), make_tuple(1298320,1,6493), make_tuple(376309,1,1882), make_tuple(531609,1,2659), make_tuple(1079861,1,5401), make_tuple(1985626,1,9931), make_tuple(1633391,1,8169), make_tuple(165280,1,827), make_tuple(1540502,1,7705), make_tuple(1884172,1,9423), make_tuple(1889900,1,9452), make_tuple(269992,1,1351), make_tuple(959573,1,4799), make_tuple(1747155,1,8738), make_tuple(1721745,1,8611), make_tuple(1245035,1,6227), make_tuple(1461139,1,7308), make_tuple(1544368,1,7724), make_tuple(550658,1,2754), make_tuple(955641,1,4780), make_tuple(1910050,1,9553), make_tuple(223562,1,1118), make_tuple(833764,1,4170), make_tuple(683537,1,3419), make_tuple(673266,1,3368), make_tuple(1167363,1,5839), make_tuple(1610550,1,8055), make_tuple(660413,1,3303), make_tuple(1190638,1,5955), make_tuple(364839,1,1825), make_tuple(555761,1,2780), make_tuple(690324,1,3453), make_tuple(151655,1,759), make_tuple(245401,1,1228), make_tuple(444485,1,2223), make_tuple(444505,1,2223), make_tuple(501108,1,2507), make_tuple(501782,1,2510), make_tuple(794450,1,3974), make_tuple(904613,1,4525), make_tuple(279138,1,1396), make_tuple(124025,1,621), make_tuple(1052305,1,5263), make_tuple(1797484,1,8990), make_tuple(1648309,1,8244), make_tuple(1716336,1,8584), make_tuple(1486415,1,7434), make_tuple(1120047,1,5602), make_tuple(1198546,1,5994), make_tuple(1495523,1,7480), make_tuple(1421348,1,7109), make_tuple(1487223,1,7438), make_tuple(1494869,1,7476), make_tuple(532366,1,2663), make_tuple(710567,1,3554), make_tuple(1118795,1,5596), make_tuple(565482,1,2829), make_tuple(1013563,1,5069), make_tuple(936777,1,4685), make_tuple(1558019,1,7792), make_tuple(1032448,1,5164), make_tuple(1261740,1,6311), make_tuple(1600052,1,8002), make_tuple(135870,1,680), make_tuple(337368,1,1688), make_tuple(519030,1,2596), make_tuple(502633,1,2514), make_tuple(1290414,1,6454), make_tuple(1752057,1,8763), make_tuple(144449,1,723), make_tuple(625878,1,3131), make_tuple(1308265,1,6543), make_tuple(576052,1,2881), make_tuple(1079510,1,5399), make_tuple(1446672,1,7235), make_tuple(1754018,1,8772), make_tuple(1873028,1,9368), make_tuple(1954423,1,9775)};
//Wjets
//set<tuple<long,long,long>>  inspection_set_erl = {};
//WWW
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(42,221,1), make_tuple(39,257,1), make_tuple(2,275,1), make_tuple(14,289,1), make_tuple(23,331,1), make_tuple(56,332,2), make_tuple(43,333,1), make_tuple(59,119,2), make_tuple(96,361,2), make_tuple(47,36,1), make_tuple(13,402,1), make_tuple(64,406,2), make_tuple(3,410,1), make_tuple(6,420,1), make_tuple(20,427,1), make_tuple(86,428,2), make_tuple(33,448,1), make_tuple(3,484,1), make_tuple(58,501,2), make_tuple(81,504,2), make_tuple(9,522,1), make_tuple(42,530,1), make_tuple(70,588,2), make_tuple(55,592,2), make_tuple(4,607,1), make_tuple(82,141,2), make_tuple(72,150,2), make_tuple(83,638,2), make_tuple(95,64,2), make_tuple(41,667,1), make_tuple(55,679,2), make_tuple(70,689,2), make_tuple(98,701,2), make_tuple(10,703,1), make_tuple(79,706,2), make_tuple(82,746,2), make_tuple(5,750,1), make_tuple(12,751,1), make_tuple(33,760,1), make_tuple(69,765,2), make_tuple(72,770,2), make_tuple(23,793,1), make_tuple(25,845,1), make_tuple(34,85,1), make_tuple(62,887,2), make_tuple(16,170,1), make_tuple(84,918,2), make_tuple(42,947,1), make_tuple(76,966,2), make_tuple(72,983,2), make_tuple(51,19,2), make_tuple(70,19,2), make_tuple(70,820,2), make_tuple(20,833,1), make_tuple(20,842,1), make_tuple(53,849,2), make_tuple(62,866,2), make_tuple(26,880,1), make_tuple(96,1142,2), make_tuple(64,935,2), make_tuple(51,948,2), make_tuple(48,968,1), make_tuple(60,983,2), make_tuple(51,990,2), make_tuple(88,993,2), make_tuple(10,116,1), make_tuple(45,1183,1), make_tuple(71,1201,2), make_tuple(65,1248,2), make_tuple(81,126,2), make_tuple(17,1020,1), make_tuple(58,1024,2), make_tuple(84,1287,2), make_tuple(41,1326,1), make_tuple(88,1328,2), make_tuple(54,1342,2), make_tuple(94,134,2), make_tuple(29,1360,1), make_tuple(85,1371,2), make_tuple(88,138,2), make_tuple(86,1413,2), make_tuple(58,1415,2), make_tuple(21,1424,1), make_tuple(59,1440,2), make_tuple(60,1440,2), make_tuple(73,1444,2), make_tuple(85,1040,2), make_tuple(51,150,2), make_tuple(4,1550,1), make_tuple(97,1670,2), make_tuple(94,1674,2), make_tuple(14,1054,1), make_tuple(28,1704,1), make_tuple(28,173,1), make_tuple(61,1741,2), make_tuple(59,1779,2), make_tuple(3,1781,1), make_tuple(12,1810,1), make_tuple(56,1839,2), make_tuple(9,183,1), make_tuple(36,1057,1), make_tuple(62,1064,2), make_tuple(20,1895,1), make_tuple(48,190,1), make_tuple(69,1924,2), make_tuple(49,1951,1), make_tuple(88,1977,2), make_tuple(96,1981,2), make_tuple(15,234,1), make_tuple(33,249,1), make_tuple(22,270,1), make_tuple(7,307,1), make_tuple(92,332,2), make_tuple(62,335,2), make_tuple(77,1095,2), make_tuple(77,380,2), make_tuple(100,38,2), make_tuple(79,419,2), make_tuple(24,442,1), make_tuple(33,486,1), make_tuple(66,1096,2), make_tuple(21,511,1), make_tuple(45,552,1), make_tuple(8,566,1), make_tuple(85,575,2), make_tuple(94,592,2), make_tuple(35,59,1), make_tuple(91,601,2), make_tuple(92,606,2), make_tuple(31,613,1), make_tuple(12,618,1), make_tuple(58,627,2), make_tuple(56,640,2), make_tuple(17,1117,1), make_tuple(50,667,1), make_tuple(18,686,1), make_tuple(17,687,1), make_tuple(24,689,1), make_tuple(85,691,2), make_tuple(29,725,1), make_tuple(64,739,2), make_tuple(63,745,2), make_tuple(10,752,1), make_tuple(43,757,1)};
//---------
// 2SFOS
//---------
//tt1l
//set<tuple<long,long,long>>  inspection_set_erl = {};
//WZ
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(1067496,1,5339), make_tuple(1783520,1,8920), make_tuple(1591339,1,7959), make_tuple(362878,1,1815), make_tuple(805669,1,4030), make_tuple(1581054,1,7907), make_tuple(810619,1,4054), make_tuple(1968205,1,9844), make_tuple(1989888,1,9952), make_tuple(611008,1,3056), make_tuple(825488,1,4129), make_tuple(490909,1,2456), make_tuple(818639,1,4095), make_tuple(239158,1,1197), make_tuple(180411,1,903), make_tuple(629362,1,3148), make_tuple(808257,1,4043), make_tuple(108840,1,545), make_tuple(167012,1,836), make_tuple(262537,1,1313), make_tuple(262577,1,1314), make_tuple(1502941,1,7517), make_tuple(146118,1,731), make_tuple(964878,1,4826), make_tuple(589212,1,2947), make_tuple(607680,1,3039), make_tuple(776626,1,3884), make_tuple(28345,1,142), make_tuple(1887676,1,9441), make_tuple(1241420,1,6209), make_tuple(15478,1,78), make_tuple(113678,1,569), make_tuple(1728847,1,8647), make_tuple(19214,1,97), make_tuple(122725,1,614), make_tuple(969067,1,4847), make_tuple(388614,1,1944), make_tuple(1152260,1,5763), make_tuple(1453273,1,7268), make_tuple(1310096,1,6552), make_tuple(199591,1,999), make_tuple(1796630,1,8986), make_tuple(1760240,1,8804), make_tuple(492339,1,2463), make_tuple(1136073,1,5682), make_tuple(1253408,1,6269), make_tuple(1260013,1,6302), make_tuple(1864340,1,9324), make_tuple(1882534,1,9415), make_tuple(1718703,1,8596), make_tuple(1346852,1,6736), make_tuple(1045410,1,5229), make_tuple(1543892,1,7722), make_tuple(1703944,1,8522), make_tuple(161722,1,809), make_tuple(246153,1,1232), make_tuple(1108952,1,5546), make_tuple(1109445,1,5549), make_tuple(1116916,1,5586), make_tuple(1548746,1,7746), make_tuple(1375939,1,6882), make_tuple(344370,1,1723), make_tuple(1353835,1,6771), make_tuple(572944,1,2866), make_tuple(993485,1,4969), make_tuple(1475107,1,7378), make_tuple(143,1,1), make_tuple(489071,1,2446), make_tuple(149200,1,747), make_tuple(407201,1,2037), make_tuple(738443,1,3693), make_tuple(1528437,1,7644), make_tuple(62642,1,314), make_tuple(507786,1,2540), make_tuple(267049,1,1336), make_tuple(1017142,1,5087), make_tuple(1084531,1,5424), make_tuple(148197,1,742), make_tuple(578874,1,2895), make_tuple(835177,1,4177), make_tuple(493514,1,2469), make_tuple(1723647,1,8621), make_tuple(52319,1,262), make_tuple(415674,1,2079), make_tuple(594900,1,2976), make_tuple(624818,1,3125), make_tuple(515269,1,2577), make_tuple(1921546,1,9610), make_tuple(1993042,1,9968), make_tuple(1133080,1,5667), make_tuple(1486713,1,7436), make_tuple(1387628,1,6940), make_tuple(1118959,1,5596), make_tuple(383749,1,1920), make_tuple(667370,1,3338), make_tuple(182808,1,915), make_tuple(1417143,1,7088), make_tuple(1924373,1,9624)};
//Wjets
//set<tuple<long,long,long>>  inspection_set_erl = {};
//WWW
//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(93,108,2), make_tuple(9,258,1), make_tuple(41,262,1), make_tuple(66,401,2), make_tuple(43,138,1), make_tuple(61,545,2), make_tuple(92,566,2), make_tuple(49,591,1), make_tuple(52,622,2), make_tuple(96,654,2), make_tuple(21,664,1), make_tuple(69,731,2), make_tuple(63,156,2), make_tuple(99,779,2), make_tuple(70,784,2), make_tuple(66,805,2), make_tuple(86,80,2), make_tuple(4,817,1), make_tuple(16,824,1), make_tuple(92,864,2), make_tuple(62,882,2), make_tuple(26,889,1), make_tuple(35,173,1), make_tuple(54,928,2), make_tuple(31,977,1), make_tuple(55,999,2), make_tuple(33,782,1), make_tuple(21,795,1), make_tuple(16,868,1), make_tuple(75,901,2), make_tuple(19,1147,1), make_tuple(26,91,1), make_tuple(98,940,2), make_tuple(80,980,2), make_tuple(47,988,1), make_tuple(16,990,1), make_tuple(63,115,2), make_tuple(54,1195,2), make_tuple(69,129,2), make_tuple(36,132,1), make_tuple(9,1414,1), make_tuple(35,1459,1), make_tuple(82,1477,2), make_tuple(83,1502,2), make_tuple(73,1549,2), make_tuple(5,1566,1), make_tuple(58,1571,2), make_tuple(40,1579,1), make_tuple(42,1681,1), make_tuple(10,1707,1), make_tuple(61,1864,2), make_tuple(24,187,1), make_tuple(9,1956,1), make_tuple(45,225,1), make_tuple(22,256,1), make_tuple(74,262,2), make_tuple(76,294,2), make_tuple(80,108,2), make_tuple(95,35,2), make_tuple(75,367,2), make_tuple(85,396,2), make_tuple(31,40,1), make_tuple(3,413,1), make_tuple(59,48,2), make_tuple(2,525,1), make_tuple(18,595,1), make_tuple(80,600,2), make_tuple(60,647,2), make_tuple(93,64,2), make_tuple(97,718,2)};


//set<tuple<long,long,long>>  inspection_set_erl = {make_tuple(12048681,1,37025), make_tuple(4356035,1,22535), make_tuple(4897488,1,25336), make_tuple(1881256,1,9732), make_tuple(2327643,1,12041), make_tuple(4611795,1,23858), make_tuple(3974725,1,20563), make_tuple(4629795,1,23951), make_tuple(3833457,1,19832), make_tuple(2962380,1,15325), make_tuple(5188918,1,26844), make_tuple(4861517,1,25150), make_tuple(3690950,1,19094), make_tuple(4355161,1,22530), make_tuple(1589698,1,8224), make_tuple(5701134,1,29493), make_tuple(2489220,1,12877), make_tuple(5484982,1,28375), make_tuple(5780797,1,29905), make_tuple(1004906,1,5199), make_tuple(556292,1,2878), make_tuple(4579806,1,23693), make_tuple(3703246,1,19157), make_tuple(5177168,1,26783), make_tuple(6444163,1,33337), make_tuple(970728,1,5022), make_tuple(5741549,1,29702), make_tuple(4842606,1,25052)};

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

vector<short> g_lep_inds; //Holds poition of "analysis leptons"
short g_nlep;
vector<bool> g_looseIDs;


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