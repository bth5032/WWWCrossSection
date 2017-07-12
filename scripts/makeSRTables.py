import argparse, os, sys, root_utils, math

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-x", "--config", help="Config directory in which to look for histograms", type=str, default="Sync")
parser.add_argument("-r", "--raw", help="Get raw counts instead of yields from table", action="store_true")
parser.add_argument("-h", "--help", help="Print help message and quit", action="store_true")

args=parser.parse_args()

if (args.help):
  parser.print_help()
  exit(0)

import ROOT as r

SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]

def PrintSRTable(yields, bg_samples, signal_samples, pretty_names):
  print("\\begin{table}[ht!]")
  print("\\begin{center}")
  print("\\begin{tabular}{l|c|c|c|c|c|c}")
  print("Sample & ee & e$\mu$ & $\mu\mu$ & 0SFOS & 1SFOS & 2SFOS \\\\ \\hline")
  
  #Setup for sum of BG and Signal counting.
  bg_count = {}
  sig_count = {}

  for sr in SRs:
    bg_count[sr] = {}
    sig_count[sr] = {}
    bg_count[sr]["cv"] = 0
    bg_count[sr]["unc"] = 0
    sig_count[sr]["cv"] = 0
    sig_count[sr]["unc"] = 0

  for sample in bg_samples:
    print("%s & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (pretty_names[sample], yields[sample]["2lepSSEE"]["cv"], yields[sample]["2lepSSEE"]["unc"], yields[sample]["2lepSSEMu"]["cv"], yields[sample]["2lepSSEMu"]["unc"], yields[sample]["2lepSSMuMu"]["cv"], yields[sample]["2lepSSMuMu"]["unc"], yields[sample]["3lep_0SFOS"]["cv"], yields[sample]["3lep_0SFOS"]["unc"], yields[sample]["3lep_1SFOS"]["cv"], yields[sample]["3lep_1SFOS"]["unc"], yields[sample]["3lep_2SFOS"]["cv"], yields[sample]["3lep_2SFOS"]["unc"] ) )
    
    #Every time I get a new sample, add its count and unc (in quadrature) to the bg_count dict.
    for sr in SRs:
      bg_count[sr]["cv"] += yields[sample][sr]["cv"]
      bg_count[sr]["unc"] += yields[sample][sr]["unc"]*yields[sample][sr]["unc"]
  

  print("\\hline")
  
  for sample in signal_samples:
    print("%s & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (pretty_names[sample], yields[sample]["2lepSSEE"]["cv"], yields[sample]["2lepSSEE"]["unc"], yields[sample]["2lepSSEMu"]["cv"], yields[sample]["2lepSSEMu"]["unc"], yields[sample]["2lepSSMuMu"]["cv"], yields[sample]["2lepSSMuMu"]["unc"], yields[sample]["3lep_0SFOS"]["cv"], yields[sample]["3lep_0SFOS"]["unc"], yields[sample]["3lep_1SFOS"]["cv"], yields[sample]["3lep_1SFOS"]["unc"], yields[sample]["3lep_2SFOS"]["cv"], yields[sample]["3lep_2SFOS"]["unc"] ) )  
    #Every time I get a new sample, add its count and unc (in quadrature) to the sig_count dict.
    for sr in SRs:
      sig_count[sr]["cv"] += yields[sample][sr]["cv"]
      sig_count[sr]["unc"] += yields[sample][sr]["unc"]*yields[sample][sr]["unc"]
  
  #After all the samples are counted, go through and take the sqrt of the uncs so that they are added in quadrature properly
  for sr in SRs:
    bg_count[sr]["unc"] = math.sqrt(bg_count[sr]["unc"])
    sig_count[sr]["unc"] = math.sqrt(sig_count[sr]["unc"])
  
  #Print the BG and Signal sums
  print("\\hline")
  print("%s & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % ("Background Sum", bg_count["2lepSSEE"]["cv"], bg_count["2lepSSEE"]["unc"], bg_count["2lepSSEMu"]["cv"], bg_count["2lepSSEMu"]["unc"], bg_count["2lepSSMuMu"]["cv"], bg_count["2lepSSMuMu"]["unc"], bg_count["3lep_0SFOS"]["cv"], bg_count["3lep_0SFOS"]["unc"], bg_count["3lep_1SFOS"]["cv"], bg_count["3lep_1SFOS"]["unc"], bg_count["3lep_2SFOS"]["cv"], bg_count["3lep_2SFOS"]["unc"] ) )
  print("%s & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % ("Signal Sum", sig_count["2lepSSEE"]["cv"], sig_count["2lepSSEE"]["unc"], sig_count["2lepSSEMu"]["cv"], sig_count["2lepSSEMu"]["unc"], sig_count["2lepSSMuMu"]["cv"], sig_count["2lepSSMuMu"]["unc"], sig_count["3lep_0SFOS"]["cv"], sig_count["3lep_0SFOS"]["unc"], sig_count["3lep_1SFOS"]["cv"], sig_count["3lep_1SFOS"]["unc"], sig_count["3lep_2SFOS"]["cv"], sig_count["3lep_2SFOS"]["unc"] ) )

  print("\\end{tabular}")
  print("\\end{center}")
  print("\\end{table}")

def getYieldsFromSample(hist_loc):
  """Takes in the location to a histogram and gets the integral and error of the MET histogram through its entire range."""
  f = r.TFile(hist_loc)
  h = f.Get("type1MET").Clone("met")
  unc=r.Double()
  cv=h.IntegralAndError(1,-1, unc)
  return (cv, float(unc))

def getRawYieldsFromSample(hist_loc):
  """Takes in the location to a histogram and gets the integral and error of the MET histogram through its entire range."""
  f = r.TFile(hist_loc)
  h = f.Get("numEvents").Clone("n_events")
  cv = h.GetBinContent(1)
  cv -= h.Integral(2,-1)
  unc = math.sqrt(cv)
  return (cv, float(unc))

def getSampleYields(samples, base_hists_path, raw):
  """Constructs and returns the yields dictionary used in PrintSampleTable. Goes through each Sample and stores the yields in a dict"""
  yields = {}

  for sample in samples:
    yields[sample] = {}
    for sr in SRs:
      if raw:
        cv, unc = getRawYieldsFromSample(base_hists_path+sr+"/"+sample+".root")
      else:
        cv, unc = getYieldsFromSample(base_hists_path+sr+"/"+sample+".root")
      yields[sample][sr] = {"cv": cv, "unc": unc}

  return yields

def main():
  bg_samples=["VVV", "WZ", "WW", "ZZ", "WJets", "DY", "TTV", "TTBar2l", "TTBar1l", "SingleTop", "Other"]
  signal_samples=["WWW", "Wh"]

  pretty_names={
  "WZ": "WZ", 
  "WW": "WW", 
  "WJets": "Wjets", 
  "ZZ": "ZZ", 
  "TTBar1l": "TTBar $\\to$ 1 lep", 
  "TTBar2l": "TTBar $\\to$ 2 lep", 
  "DY": "Zjets", 
  "VVV": "VVV", 
  "TTV": "TTV", 
  "SingleTop": "Single Top", 
  "Other": "ZH",
  "WWW": "WWW (signal)",
  "Wh": "WH $\\to$ WWW$^*$ (signal)"}

  samples=bg_samples+signal_samples


  raw=False
  if (args.raw):
    raw=True

  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/%s/" % args.config
  
  yields = getSampleYields(samples, base_hists_path, raw)
  PrintSRTable(yields, bg_samples, signal_samples, pretty_names)

if __name__ == "__main__":
  main()