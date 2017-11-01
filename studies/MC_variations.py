#!/usr/bin/env python

"""Creates Histogram for the Fake Rate Clsoure Test.

Before this script is run, the configs in FRClosure/<Conf>/ need to be made...
This script will output new histograms in the FRClosure/<Conf>/Combined/. 
One with the Fake rate prediction called 'fakes.root' and one for each sample."""

import argparse, os, sys, math, array, ctypes, itertools
import ROOT as r
r.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-u", "--usage", help="Print help message and quit", action="store_true")

args=parser.parse_args()

if (args.usage):
  parser.print_help()
  exit(0)

variations = {}

def err_ratio(num, den, Dnum, Dden):
  """Returns the 1 sigma on a ratio"""
  return abs((float(den)*Dnum - float(num)*Dden)/(den*den))
  
def fillYields(filename, sr):
  """Computes the number of objects in the pt/eta bins and applies the FR given with the error."""

  f = r.TFile(filename, "read")
  hist = f.Get("MC_variations").Clone("h_%s" % sr)

  variations[sr] = {
    "tot" : hist->GetBinContnent(1),
    "high_FR" : hist->GetBinContnent(2),
    "low_FR" : hist->GetBinContnent(2),
    "f0p5r0p5" : hist->GetBinContnent(9),
    "f2r2" : hist->GetBinContnent(5),
    "alpha_up" : hist->GetBinContnent(10),
    "alpha_dn" : hist->GetBinContnent(11),
    "pdf_up" : hist->GetBinContnent(12),
    "pdf_dn" : hist->GetBinContnent(13)
  }

  for i in xrange(3,11):
    count = hist->GetBinContnent(i)
    if count > variations[sr]["high_FR"]:
      print("high for %s came from %d" % (sr, i) )
      variations[sr]["high_FR"] = count
    else if count < variations[sr]["low_FR"]:
      variations[sr]["low_FR"] = count
      print("low for %s came from %d" % (sr, i) )

def printLatexTable(SRs):
  outfile = open("outputs/MC_Variation.tex", "w+")

  pretty_SR_names = ["2lepSSEE": "$ee$","2lepSSEMu": "$e \mu$","2lepSSMuMu":  "$\mu \mu$","3lep_0SFOS": "0SFOS","3lep_1SFOS": "1SFOS","3lep_2SFOS": "2SFOS"]

  outfile.write("\\begin{table}[ht!]\n")
  outfile.write("\\begin{center}\n")
  outfile.write("\\begin{tabular}{l"+"|c"*len(line1.keys())+"}\n")
  
  var_names = ["tot", "high_FR", "low_FR", "f0p5r0p5", "f2r2", "alpha_up", "alpha_dn", "pdf_up", "pdf_dn"]

  header="Variation "
  for v in var_names:
    header+="& %s " % v
  header+="\\\\ \\hline\n"

  outfile.write(header)

  ##==============
  ## Variations
  ##==============
  for sr in SRs:
    line = "%s" % pretty_SR_names[sr]
    for v in var_names:
      line+=" & $%0.2f (%0.2f)$" % (variations[sr][v], variations[sr][v]/variations[sr][0])
    line+= "\\\\ \n"
    outfile.write(line)
    
  outfile.write("\\end{tabular} \n")
  outfile.write("\\end{center} \n")
  outfile.write("\\end{table} \n")

def main():
  SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/Sync/SRs/"
  

  for signal_region in SRs:
    fillYields("%s/%s/WWW.root" % (base_hists_path, signal_region), signal_region)

  printLatexTable(SRs)

if __name__ == "__main__":
  main()