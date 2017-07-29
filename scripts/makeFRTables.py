import argparse, os, sys, root_utils, math

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-s", "--study_dir", help="The config directory name in FRStudy, e.g. LooseIso", type=str, default="Baseline")
parser.add_argument("--sample_table", help="Make table in format broken up by samples", action="store_true")
parser.add_argument("--txt", help="Print table in txt format", action="store_true")

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
parser.add_argument("--wh", help="Use WH sample", action="store_true")

parser.add_argument("-h", "--help", help="Print help message and quit", action="store_true")

args=parser.parse_args()

if (args.help):
  parser.print_help()
  exit(0)

import ROOT as r

signals = ["WWW", "Wh"]
bgs = ["TTBar2l", "WZ", "ZZ", "WW", "TTV", "VVV", "TTBar1l", "SingleTop", "WJets", "DY", "Other"]

FR_enum = {"T_real": 0, "T_fake": 1, "TT_real": 2, "TT_fake": 3, "TTT_real": 4, "TTT_fake": 5, "L_real": 6, "L_fake": 7, "CATERR": 8};

real_tight_bins = {
  "2lepSSEE"  : [FR_enum["TT_real"]],
  "2lepSSEMu" : [FR_enum["TT_real"]],
  "2lepSSMuMu": [FR_enum["TT_real"]],
  "3lep_0SFOS": [FR_enum["TTT_real"]],
  "3lep_1SFOS": [FR_enum["TTT_real"]],
  "3lep_2SFOS": [FR_enum["TTT_real"]]
}

"""real_loose_bins = {
  "2lepSSEE"  : [FR_enum["T_real"], FR_enum["L_real"]],
  "2lepSSEMu" : [FR_enum["T_real"], FR_enum["L_real"]],
  "2lepSSMuMu": [FR_enum["T_real"], FR_enum["L_real"]],
  "3lep_0SFOS": [FR_enum["TT_real"],FR_enum["T_real"],FR_enum["L_real"]],
  "3lep_1SFOS": [FR_enum["TT_real"],FR_enum["T_real"],FR_enum["L_real"]],
  "3lep_2SFOS": [FR_enum["TT_real"],FR_enum["T_real"],FR_enum["L_real"]]
}"""

real_loose_bins = {
  "2lepSSEE"  : [FR_enum["T_real"]],
  "2lepSSEMu" : [FR_enum["T_real"]],
  "2lepSSMuMu": [FR_enum["T_real"]],
  "3lep_0SFOS": [FR_enum["TT_real"]],
  "3lep_1SFOS": [FR_enum["TT_real"]],
  "3lep_2SFOS": [FR_enum["TT_real"]]
}

fake_tight_bins = {
  "2lepSSEE"  : [FR_enum["TT_fake"]],
  "2lepSSEMu" : [FR_enum["TT_fake"]],
  "2lepSSMuMu": [FR_enum["TT_fake"]],
  "3lep_0SFOS": [FR_enum["TTT_fake"]],
  "3lep_1SFOS": [FR_enum["TTT_fake"]],
  "3lep_2SFOS": [FR_enum["TTT_fake"]]
}

"""fake_loose_bins = {
  "2lepSSEE"  : [FR_enum["T_fake"], FR_enum["L_fake"]],
  "2lepSSEMu" : [FR_enum["T_fake"], FR_enum["L_fake"]],
  "2lepSSMuMu": [FR_enum["T_fake"], FR_enum["L_fake"]],
  "3lep_0SFOS": [FR_enum["TT_fake"],FR_enum["T_fake"],FR_enum["L_fake"]],
  "3lep_1SFOS": [FR_enum["TT_fake"],FR_enum["T_fake"],FR_enum["L_fake"]],
  "3lep_2SFOS": [FR_enum["TT_fake"],FR_enum["T_fake"],FR_enum["L_fake"]]
}"""

fake_loose_bins = {
  "2lepSSEE"  : [FR_enum["T_fake"]],
  "2lepSSEMu" : [FR_enum["T_fake"]],
  "2lepSSMuMu": [FR_enum["T_fake"]],
  "3lep_0SFOS": [FR_enum["TT_fake"]],
  "3lep_1SFOS": [FR_enum["TT_fake"]],
  "3lep_2SFOS": [FR_enum["TT_fake"]]
}

def PrintSRTable(yields, study_dir, latex):
  print("Using Config: %s" % study_dir)
  print("\\begin{table}[ht!]")
  print("\\begin{center}")
  print("\\begin{tabular}{l|c|c|c|c}")
  print("Signal Region & Real/Tight & Real/Loose & Fake/Tight & Fake/Loose \\\\ \\hline")
  
  print("SSee          & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (yields["2lepSSEE"]["rt"],yields["2lepSSEE"]["rt_unc"],yields["2lepSSEE"]["rl"],yields["2lepSSEE"]["rl_unc"],yields["2lepSSEE"]["ft"],yields["2lepSSEE"]["ft_unc"],yields["2lepSSEE"]["fl"],yields["2lepSSEE"]["fl_unc"]) )
  print("SSe$\mu$      & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (yields["2lepSSEMu"]["rt"],yields["2lepSSEMu"]["rt_unc"],yields["2lepSSEMu"]["rl"],yields["2lepSSEMu"]["rl_unc"],yields["2lepSSEMu"]["ft"],yields["2lepSSEMu"]["ft_unc"],yields["2lepSSEMu"]["fl"],yields["2lepSSEMu"]["fl_unc"]) )
  print("SS$\mu\mu$    & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (yields["2lepSSMuMu"]["rt"],yields["2lepSSMuMu"]["rt_unc"],yields["2lepSSMuMu"]["rl"],yields["2lepSSMuMu"]["rl_unc"],yields["2lepSSMuMu"]["ft"],yields["2lepSSMuMu"]["ft_unc"],yields["2lepSSMuMu"]["fl"],yields["2lepSSMuMu"]["fl_unc"]) )
  print("3lep 0SFOS    & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (yields["3lep_0SFOS"]["rt"],yields["3lep_0SFOS"]["rt_unc"],yields["3lep_0SFOS"]["rl"],yields["3lep_0SFOS"]["rl_unc"],yields["3lep_0SFOS"]["ft"],yields["3lep_0SFOS"]["ft_unc"],yields["3lep_0SFOS"]["fl"],yields["3lep_0SFOS"]["fl_unc"]) )
  print("3lep 1SFOS    & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (yields["3lep_1SFOS"]["rt"],yields["3lep_1SFOS"]["rt_unc"],yields["3lep_1SFOS"]["rl"],yields["3lep_1SFOS"]["rl_unc"],yields["3lep_1SFOS"]["ft"],yields["3lep_1SFOS"]["ft_unc"],yields["3lep_1SFOS"]["fl"],yields["3lep_1SFOS"]["fl_unc"]) )
  print("3lep 2SFOS    & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (yields["3lep_2SFOS"]["rt"],yields["3lep_2SFOS"]["rt_unc"],yields["3lep_2SFOS"]["rl"],yields["3lep_2SFOS"]["rl_unc"],yields["3lep_2SFOS"]["ft"],yields["3lep_2SFOS"]["ft_unc"],yields["3lep_2SFOS"]["fl"],yields["3lep_2SFOS"]["fl_unc"]) ) 

  print("\\end{tabular}")
  print("\\end{center}")
  print("\\end{table}")

def PrintSampleTable(yields, study_dir, latex):
  print("Using Config: %s" % study_dir)
  SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]

  for sr in SRs:
    if (latex):
      print("")
      print("\\begin{table}[ht!]")
      print("\\begin{center}")
      print("")
      print("\\textbf{%s}" % sr.replace('_', ' '))
      print("")
      print("\\begin{tabular}{l|c|c|c|c}")
      print("Sample Name & Real/Tight & Real/Loose & Fake/Tight & Fake/Loose \\\\ \\hline")
   
      for s in bgs: #for each BG sample
        if s in yields[sr]:
          print("%s          & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (s, yields[sr][s]["rt"],yields[sr][s]["rt_unc"],yields[sr][s]["rl"],yields[sr][s]["rl_unc"],yields[sr][s]["ft"],yields[sr][s]["ft_unc"],yields[sr][s]["fl"],yields[sr][s]["fl_unc"]) )
      print("\\hline")
      for s in signals: #for each signal sample
        if s in yields[sr]:
          print("%s          & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (s, yields[sr][s]["rt"],yields[sr][s]["rt_unc"],yields[sr][s]["rl"],yields[sr][s]["rl_unc"],yields[sr][s]["ft"],yields[sr][s]["ft_unc"],yields[sr][s]["fl"],yields[sr][s]["fl_unc"]) )
      print("\\hline")
      s="signal"
      print("%s          & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (s, yields[sr][s]["rt"],yields[sr][s]["rt_unc"],yields[sr][s]["rl"],yields[sr][s]["rl_unc"],yields[sr][s]["ft"],yields[sr][s]["ft_unc"],yields[sr][s]["fl"],yields[sr][s]["fl_unc"]) )
      s="bg"
      print("%s          & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (s, yields[sr][s]["rt"],yields[sr][s]["rt_unc"],yields[sr][s]["rl"],yields[sr][s]["rl_unc"],yields[sr][s]["ft"],yields[sr][s]["ft_unc"],yields[sr][s]["fl"],yields[sr][s]["fl_unc"]) )
   
      print("\\end{tabular}")
      print("\\end{center}")
      print("\\end{table}")
      print("")
    else:
      print("+---------------------------------------------+")
      print("SR: %s" % sr.replace('_', ' '))
      print("+---------------------------------------------+")
      print("Sample  Real/Tight  Real/Loose  Fake/Tight  Fake/Loose ")
      for s in yields[sr]: #for each sample
        print("%s %0.2f+/-%0.2f %0.2f+/-%0.2f %0.2f+/-%0.2f %0.2f+/-%0.2f" % (s, yields[sr][s]["rt"],yields[sr][s]["rt_unc"],yields[sr][s]["rl"],yields[sr][s]["rl_unc"],yields[sr][s]["ft"],yields[sr][s]["ft_unc"],yields[sr][s]["fl"],yields[sr][s]["fl_unc"]) )
      print("")

def getYieldsFromSample(hist_loc, SR):
  """Takes in Signal region and the location to a histogram, looks up the proper bins for what is called real tight, real loose, fake tight, fake loose, and returns the counts for the sample in a tuple (rt, rl, ft, fl)"""
  f = r.TFile(hist_loc)
  h = f.Get("fr_counts")
  rt = rl = ft = fl = rt_unc = rl_unc = ft_unc = fl_unc = 0 

  for b in real_tight_bins[SR]:
    rt+=h.GetBinContent(h.FindBin(b))
    rt_unc+=h.GetBinError(h.FindBin(b))*h.GetBinError(h.FindBin(b))
  for b in real_loose_bins[SR]:
    rl+=h.GetBinContent(h.FindBin(b))
    rl_unc+=h.GetBinError(h.FindBin(b))*h.GetBinError(h.FindBin(b))
  for b in fake_tight_bins[SR]:
    ft+=h.GetBinContent(h.FindBin(b))
    ft_unc+=h.GetBinError(h.FindBin(b))*h.GetBinError(h.FindBin(b))
  for b in fake_loose_bins[SR]:
    fl+=h.GetBinContent(h.FindBin(b))    
    fl_unc+=h.GetBinError(h.FindBin(b))*h.GetBinError(h.FindBin(b))

  rt_unc = math.sqrt(rt_unc)
  rl_unc = math.sqrt(rl_unc)
  ft_unc = math.sqrt(ft_unc)
  fl_unc = math.sqrt(fl_unc)

  return (rt, rl, ft, fl, rt_unc, rl_unc, ft_unc, fl_unc)

def getCombinedYields(samples, study_dir):
  """Constructs and returns the yields dictionary used in PrintTable. Goes through each SR and adds the yields for each sample in that SR and organizes them in the dict."""
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRStudy/%s/" % study_dir
  SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]

  yields = {}

  for sr in SRs:
    rt = rl = ft = fl = rt_unc = rl_unc = ft_unc = fl_unc = 0 
    for s in samples:
      rt_, rl_, ft_, fl_, rt_unc_, rl_unc_, ft_unc_, fl_unc_ = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr)
      rt += rt_
      rl += rl_
      ft += ft_
      fl += fl_
      rt_unc+=rt_unc_*rt_unc_
      rl_unc+=rl_unc_*rl_unc_
      ft_unc+=ft_unc_*ft_unc_
      fl_unc+=fl_unc_*fl_unc_

    rt_unc = math.sqrt(rt_unc)
    rl_unc = math.sqrt(rl_unc)
    ft_unc = math.sqrt(ft_unc)
    fl_unc = math.sqrt(fl_unc)

    yields[sr] = {"rt": rt, "rl": rl, "ft": ft, "fl": fl, "rt_unc": rt_unc, "rl_unc": rl_unc, "ft_unc": ft_unc, "fl_unc": fl_unc}

  return yields

def getSampleYields(samples, study_dir):
  """Constructs and returns the yields dictionary used in PrintSampleTable. Goes through each SR and adds the yields to the dict for each sample (they are kept seperate as opposed to getCombinedYields) in that SR and organizes them in the dict."""
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRStudy/%s/" % study_dir
  SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]

  yields = {}

  for sr in SRs:
    yields[sr] = {}
    yields[sr]["signal"] = {"rt": 0, "rl": 0, "ft": 0, "fl": 0, "rt_unc": 0, "rl_unc": 0, "ft_unc": 0, "fl_unc": 0}
    yields[sr]["bg"] = {"rt": 0, "rl": 0, "ft": 0, "fl": 0, "rt_unc": 0, "rl_unc": 0, "ft_unc": 0, "fl_unc": 0}
    for s in samples:
      rt, rl, ft, fl, rt_unc, rl_unc, ft_unc, fl_unc = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr)
      yields[sr][s] = {"rt": rt, "rl": rl, "ft": ft, "fl": fl, "rt_unc": rt_unc, "rl_unc": rl_unc, "ft_unc": ft_unc, "fl_unc": fl_unc}
      if s in signals:
        yields[sr]["signal"]["rt"] += rt
        yields[sr]["signal"]["rl"] += rl
        yields[sr]["signal"]["ft"] += ft
        yields[sr]["signal"]["fl"] += fl
        yields[sr]["signal"]["rt_unc"] += rt_unc*rt_unc
        yields[sr]["signal"]["rl_unc"] += rl_unc*rl_unc
        yields[sr]["signal"]["ft_unc"] += ft_unc*ft_unc
        yields[sr]["signal"]["fl_unc"] += fl_unc*fl_unc
      if s in bgs:
        yields[sr]["bg"]["rt"] += rt
        yields[sr]["bg"]["rl"] += rl
        yields[sr]["bg"]["ft"] += ft
        yields[sr]["bg"]["fl"] += fl
        yields[sr]["bg"]["rt_unc"] += rt_unc*rt_unc
        yields[sr]["bg"]["rl_unc"] += rl_unc*rl_unc
        yields[sr]["bg"]["ft_unc"] += ft_unc*ft_unc
        yields[sr]["bg"]["fl_unc"] += fl_unc*fl_unc

    yields[sr]["signal"]["rt_unc"] = math.sqrt(yields[sr]["signal"]["rt_unc"])
    yields[sr]["signal"]["rl_unc"] = math.sqrt(yields[sr]["signal"]["rl_unc"])
    yields[sr]["signal"]["ft_unc"] = math.sqrt(yields[sr]["signal"]["ft_unc"])
    yields[sr]["signal"]["fl_unc"] = math.sqrt(yields[sr]["signal"]["fl_unc"])
    yields[sr]["bg"]["rt_unc"] = math.sqrt(yields[sr]["signal"]["rt_unc"])
    yields[sr]["bg"]["rl_unc"] = math.sqrt(yields[sr]["signal"]["rl_unc"])
    yields[sr]["bg"]["ft_unc"] = math.sqrt(yields[sr]["signal"]["ft_unc"])
    yields[sr]["bg"]["fl_unc"] = math.sqrt(yields[sr]["signal"]["fl_unc"])




  return yields

def main():
  samples = []
  latex=True

  #Signal
  if (args.www or args.all or args.signal):
    samples.append("WWW")
  if (args.wh or args.all or args.signal):
    samples.append("Wh")

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

  if (args.txt):
    latex=False


  print("Going to use %s to make table..." %samples)
  
  if (not args.sample_table):
    yields = getCombinedYields(samples, args.study_dir)
    PrintSRTable(yields, args.study_dir, latex)
  else:
    yields = getSampleYields(samples, args.study_dir)
    PrintSampleTable(yields, args.study_dir, latex)

if __name__ == "__main__":
  main()