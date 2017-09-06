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

def buildIntervals(pt_bins, eta_bins):
  """Returns all intervals in pt and eta bins"""
  pt_intervals = [(pt_bins[i], pt_bins[i+1] - 0.001) for i in xrange(0, len(pt_bins) - 1) ]
  eta_intervals = [(eta_bins[i], eta_bins[i+1] - 0.001) for i in xrange(0, len(eta_bins) - 1) ]
  return [(i[0][0], i[0][1], i[1][0], i[1][1]) for i in itertools.product(pt_intervals, eta_intervals)]

def getYield(hist, pt_bins, eta_bins, h_fr):
  """Computes the number of objects in the pt/eta bins and applies the FR given with the error."""

  err = r.Double()
  for pt_low, pt_high, eta_low, eta_high in buildIntervals(pt_bins, eta_bins):
    x1=hist.GetXaxis().FindBin(pt_low)
    x2=hist.GetXaxis().FindBin(pt_high)
    y1=hist.GetYaxis().FindBin(eta_low)
    y2=hist.GetYaxis().FindBin(eta_high)
    print(x1,x2,y1,y2)
    print(pt_low, pt_high, eta_low, eta_high)
    y = hist.IntegralAndError(x1, x2, y1, y2, err)
    print("Count in (pt, eta) in (%0.2f, %0.2f, %0.2f, %0.2f) = %0.2f +/- %0.2f" % (pt_low, pt_high, eta_low, eta_high, y, err) )

  return (0,0)

def makeHistos(pt_bins, eta_bins_m, eta_bins_e, sample_files, FR_Hist_path):
  """Writes out the completed histogram for a list of sample files for the Fake Rate in that region."""
  fr_file = r.TFile(FR_Hist_path)
  h_m_FRs = fr_file.Get("muon_fakerate_conecorrpt_v_eta").Clone("h_m_FRs")
  h_e_FRs = fr_file.Get("elec_fakerate_conecorrpt_v_eta").Clone("h_e_FRs")

  output_dir = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/Prediction/%s" % args.study_dir

  for s_file in sample_files:
    tfile = r.TFile(s_file)
    s = s_file[s_file.rfind('/')+1:s_file.index(".root")] #extract sample name from filename
    
    print(s_file)

    h_loose_e = tfile.Get("loose_loose_pteta_e").Clone("h_loose_e_%s" % s)
    h_loose_m = tfile.Get("loose_loose_pteta_m").Clone("h_loose_m_%s" % s)

    num_e, err_e = getYield(h_loose_e, pt_bins, eta_bins_m, h_e_FRs)
    num_m, err_m = getYield(h_loose_m, pt_bins, eta_bins_m, h_m_FRs)

def main():
  samples = ["TTBar1l", "WJets"]
  #SRs = ["2lepSS", "2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
  SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu"]
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/"

  FR_Hist_path = "auxFiles/fakerate_pt_v_eta.root" 

  pt_bins = [10, 15, 20, 25, 30, 35, 50, 120, 6001]
  eta_bins_m = [0, 1.2, 2.1, 2.4, 3]
  eta_bins_e = [0, 0.8, 1.479, 2.5, 3]
  
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/%s" % args.study_dir

  for signal_region in SRs:
    hist_paths= [ "%s/%s/%s.root" % (base_hists_path, signal_region, s) for s in samples]
    makeHistos(pt_bins, eta_bins_e, eta_bins_m, hist_paths, FR_Hist_path)

if __name__ == "__main__":
  main()