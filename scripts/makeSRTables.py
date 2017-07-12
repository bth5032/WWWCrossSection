import argparse, os, sys, root_utils, math

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-x", "--config", help="Config directory in which to look for histograms", type=str, default="Sync")
parser.add_argument("-h", "--help", help="Print help message and quit", action="store_true")

args=parser.parse_args()

if (args.help):
  parser.print_help()
  exit(0)

import ROOT as r

def PrintSRTable(yields):
  print("\\begin{table}[ht!]")
  print("\\begin{center}")
  print("\\begin{tabular}{l|c|c|c|c|c|c}")
  print("Sample & EE & $\mu\mu$ & E$\mu$ & 0SFOS & 1SFOS & 2SFOS \\\\ \\hline")
  
  for sample in yields.keys():
      print("%s & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f & %0.2f $\pm$ %0.2f \\\\" % (sample, yields[sample]["ee"]["cv"], yields[sample]["ee"]["unc"], yields[sample]["mm"]["cv"], yields[sample]["mm"]["unc"], yields[sample]["em"]["cv"], yields[sample]["em"]["unc"], yields[sample]["0sfos"]["cv"], yields[sample]["0sfos"]["unc"], yields[sample]["1sfos"]["cv"], yields[sample]["1sfos"]["unc"], yields[sample]["2sfos"]["cv"], yields[sample]["2sfos"]["unc"] ) )

  print("\\end{tabular}")
  print("\\end{center}")
  print("\\end{table}")

def getYieldsFromSample(hist_loc):
  """Takes in the location to a histogram and gets the integral and error of the MET histogram through its entire range."""
  f = r.TFile(hist_loc)
  h = f.Get("t1met")
  unc=ROOT.Double()
  cv=h.IntegralAndError(1,-1, unc)
  return (cv, float(unc))

def getSampleYields(samples, base_hists_path):
  """Constructs and returns the yields dictionary used in PrintSampleTable. Goes through each Sample and stores the yields in a dict"""
  SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]

  yields = {}

  for sample in samples:
    yields[sample] = {}
    for sr in SRs:
      cv, unc = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr)
      yields[sample][sr] = {"cv": cv, "unc": unc}

  return yields

def main():
  samples=["WZ", "WW", "WJets", "ZZ", "TTBar1l", "TTBar2l", "DY", "VVV", "TTV", "SingleTop", "Other", "WWW"]
  latex=True
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/%s/" % args.config

  if (args.txt):
    latex=False
  
  yields = getSampleYields(samples, base_hists_path)
  PrintSRTable(yields)

if __name__ == "__main__":
  main()