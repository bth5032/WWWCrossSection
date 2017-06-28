import ROOT as r
import argparse, os, sys, root_utils

FR_enum = {"T_real": 0, "T_fake": 1, "TT_real": 2, "TT_fake": 3, "TTT_real": 4, "TTT_fake": 5, "L_real": 6, "L_fake": 7, "CATERR": 8};

real_tight_bins = {
  "2lepSSEE"  : [FR_enum["TT_real"]],
  "2lepSSEMu" : [FR_enum["TT_real"]],
  "2lepSSMuMu": [FR_enum["TT_real"]],
  "3lep_0SFOS": [FR_enum["TTT_real"]],
  "3lep_1SFOS": [FR_enum["TTT_real"]],
  "3lep_2SFOS": [FR_enum["TTT_real"]]
}

real_loose_bins = {
  "2lepSSEE"  : [FR_enum["T_real"], FR_enum["L_real"]],
  "2lepSSEMu" : [FR_enum["T_real"], FR_enum["L_real"]],
  "2lepSSMuMu": [FR_enum["T_real"], FR_enum["L_real"]],
  "3lep_0SFOS": [FR_enum["TT_real"],FR_enum["T_real"],FR_enum["L_real"]],
  "3lep_1SFOS": [FR_enum["TT_real"],FR_enum["T_real"],FR_enum["L_real"]],
  "3lep_2SFOS": [FR_enum["TT_real"],FR_enum["T_real"],FR_enum["L_real"]]
}

fake_tight_bins = {
  "2lepSSEE"  : [FR_enum["TT_fake"]],
  "2lepSSEMu" : [FR_enum["TT_fake"]],
  "2lepSSMuMu": [FR_enum["TT_fake"]],
  "3lep_0SFOS": [FR_enum["TTT_fake"]],
  "3lep_1SFOS": [FR_enum["TTT_fake"]],
  "3lep_2SFOS": [FR_enum["TTT_fake"]]
}

fake_loose_bins = {
  "2lepSSEE"  : [FR_enum["T_fake"], FR_enum["L_fake"]],
  "2lepSSEMu" : [FR_enum["T_fake"], FR_enum["L_fake"]],
  "2lepSSMuMu": [FR_enum["T_fake"], FR_enum["L_fake"]],
  "3lep_0SFOS": [FR_enum["TT_fake"],FR_enum["T_fake"],FR_enum["L_fake"]],
  "3lep_1SFOS": [FR_enum["TT_fake"],FR_enum["T_fake"],FR_enum["L_fake"]],
  "3lep_2SFOS": [FR_enum["TT_fake"],FR_enum["T_fake"],FR_enum["L_fake"]]
}

def PrintTable(yields, study_dir):
  print("Using Config: %s" % study_dir)
  print("\\begin{table}[ht!]")
  print("\\begin{center}")
  print("\\begin{tabular}{l|c|c|c|c}")
  print("Signal Region & Real/Tight & Real/Loose & Fake/Tight & Fake/Loose \\\\ \\hline")
  
  print("SSee          & %0.2f      & %0.2f      & %0.2f      & %0.2f     \\\\" % (yields["2lepSSEE"]["rt"],yields["2lepSSEE"]["rl"],yields["2lepSSEE"]["ft"],yields["2lepSSEE"]["fl"]) )
  print("SSe$\mu$      & %0.2f      & %0.2f      & %0.2f      & %0.2f     \\\\" % (yields["2lepSSEMu"]["rt"],yields["2lepSSEMu"]["rl"],yields["2lepSSEMu"]["ft"],yields["2lepSSEMu"]["fl"]) )
  print("SS$\mu\mu$    & %0.2f      & %0.2f      & %0.2f      & %0.2f     \\\\" % (yields["2lepSSMuMu"]["rt"],yields["2lepSSMuMu"]["rl"],yields["2lepSSMuMu"]["ft"],yields["2lepSSMuMu"]["fl"]) )
  print("3lep 0SFOS    & %0.2f      & %0.2f      & %0.2f      & %0.2f     \\\\" % (yields["3lep_0SFOS"]["rt"],yields["3lep_0SFOS"]["rl"],yields["3lep_0SFOS"]["ft"],yields["3lep_0SFOS"]["fl"]) )
  print("3lep 1SFOS    & %0.2f      & %0.2f      & %0.2f      & %0.2f     \\\\" % (yields["3lep_1SFOS"]["rt"],yields["3lep_1SFOS"]["rl"],yields["3lep_1SFOS"]["ft"],yields["3lep_1SFOS"]["fl"]) )
  print("3lep 2SFOS    & %0.2f      & %0.2f      & %0.2f      & %0.2f     \\\\" % (yields["3lep_2SFOS"]["rt"],yields["3lep_2SFOS"]["rl"],yields["3lep_2SFOS"]["ft"],yields["3lep_2SFOS"]["fl"]) )  

  print("\\end{tabular}")
  print("\\end{center}")
  print("\\end{table}")

def getYieldsFromSample(hist_loc, SR):
  """Takes in Signal region and the location to a histogram, looks up the proper bins for what is called real tight, real loose, fake tight, fake loose, and returns the counts for the sample in a tuple (rt, rl, ft, fl)"""
  f = r.TFile(hist_loc)
  h = f.Get("fr_counts")
  rt = rl = ft = fl = 0

  for b in real_tight_bins[SR]:
    rt+=h.GetBinContent(h.FindBin(b))
  for b in real_loose_bins[SR]:
    rl+=h.GetBinContent(h.FindBin(b))
  for b in fake_tight_bins[SR]:
    ft+=h.GetBinContent(h.FindBin(b))
  for b in fake_loose_bins[SR]:
    fl+=h.GetBinContent(h.FindBin(b))    

  return (rt, rl, ft, fl)

def getAllYields(samples, study_dir):
  """Constructs and returns the yields dictionary used in PrintTable. Goes through each SR and adds the yields for each sample in that SR and organizes them in the dict."""
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRStudy/%s/" % study_dir
  SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]

  yields = {}

  for sr in SRs:
    rt = rl = ft = fl = 0
    for s in samples:
      rt_, rl_, ft_, fl_ = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr)
      rt += rt_
      rl += rl_
      ft += ft_
      fl += fl_

    yields[sr] = {"rt": rt, "rl": rl, "ft": ft, "fl": fl}

  return yields

def main():
  parser = argparse.ArgumentParser(add_help=False)

  parser.add_argument("-s", "--study_dir", help="The config directory name in FRStudy, e.g. LooseIso", type=str, default="Baseline")
  
  parser.add_argument("--all", help="Use all histograms", action="store_true")
  parser.add_argument("--frbg", help="Use standard FR BG samples", action="store_true")
  parser.add_argument("--signal", help="Use signal samples", action="store_true")
  parser.add_argument("--bg", help="Use all BG samples", action="store_true")
  
  parser.add_argument("--ttbar_dilep", help="Use TTBar to dilepton sample", action="store_true")
  parser.add_argument("--ttbar_1lep", help="Use TTBar to single lepton samples", action="store_true")
  parser.add_argument("--dy", help="Use DY sample", action="store_true")
  parser.add_argument("--wz", help="Use WZ sample", action="store_true")
  parser.add_argument("--ww", help="Use WW sample", action="store_true")
  parser.add_argument("--zz", help="Use ZZ sample", action="store_true")
  parser.add_argument("--vvv", help="Use VVV sample", action="store_true")
  parser.add_argument("--wjets", help="Use W+Jets sample", action="store_true")
  parser.add_argument("--ttv", help="Use TTV sample", action="store_true")
  parser.add_argument("--singletop", help="Use TTV sample", action="store_true")
  parser.add_argument("--other", help="Use Other (VH) sample", action="store_true")
  parser.add_argument("--www", help="Use WWW sample", action="store_true")
  
  parser.add_argument("--help", help="Print help message and quit", action="store_true")

  args=parser.parse_args()

  samples = []

  #Signal
  if (args.www or args.all or args.signal):
    samples.append("WWW")

  #BG
  if (args.ttbar_dilep or args.all or args.bg):
    samples.append("TTBar2l")
  if (args.wz or args.all or args.bg):
    samples.append("WZ")
  if (args.zz or args.all or args.bg):
    samples.append("ZZ")
  if (args.ww or args.all or args.bg):
    samples.append("WW")
  if (args.ttv or args.all or args.bg):
    samples.append("TTV")    
  if (args.vvv or args.all or args.bg):
    samples.append("VVV")

  #FR BG
  if (args.ttbar_1lep or args.all or args.bg or args.frbg):
    samples.append("TTBar1l")
  if (args.singletop or args.all or args.bg or args.frbg):
    samples.append("SingleTop")
  if (args.wjets or args.all or args.bg or args.frbg):
    samples.append("WJets")
  if (args.dy or args.all or args.bg or args.frbg):
    samples.append("DY")
  if (args.other or args.all or args.bg or args.frbg):
    samples.append("Other")

  if (args.help):
    print("going to print help!")
    parser.print_help()
    exit()

  print("Going to use %s to make table..." %samples)

  yields = getAllYields(samples, args.study_dir)
  PrintTable(yields, args.study_dir)

if __name__ == "__main__":
  main()