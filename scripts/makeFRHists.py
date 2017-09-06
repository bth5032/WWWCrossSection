#!/usr/bin/env python

"""Creates Histogram for the Fake Rate Contribution to the Signal Regions using FR Histogram and Output Yields from given region."""

import argparse, os, sys, math, array, ctypes, itertools
import ROOT as r
r.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-s", "--study_dir", help="The config directory name in FRClosure, e.g. Btag", type=str, default="Btag")
parser.add_argument("-u", "--usage", help="Print help message and quit", action="store_true")

args=parser.parse_args()

if (args.usage):
  parser.print_help()
  exit(0)

def err_ratio(num, den, Dnum, Dden):
  """Returns the 1 sigma on a ratio"""
  return abs((float(den)*Dnum - float(num)*Dden)/(den*den))

def makeHistos(pt_bins, eta_bins_m, eta_bins_e, sample_files, FR_Hist_path):
  fr_file = r.TFile(FR_Hist_path)

  for s_file in sample_files:
    print(s_file)

def main():
  samples = ["TTBar1l", "WJets"]
  SRs = ["2lepSS", "2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/"

  FR_Hist_path = "auxFiles/fakerate_pt_v_eta.root" 

  pt_bins = [10, 15, 20, 25, 30, 35, 50, 120]
  eta_bins_m = [0, 1.2, 2.1, 2.4]
  eta_bins_e = [0, 0.8, 1.479, 2.5]
  
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/%s" % args.study_dir

  for signal_region in SRs:
    hist_paths= [ "%s/%s/%s.root" % (base_hists_path, signal_region, s) for s in samples]
    makeHistos(pt_bins, eta_bins_e, eta_bins_m, hist_paths, FR_Hist_path)

if __name__ == "__main__":
  main()