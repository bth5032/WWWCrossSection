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
parser.add_argument("-s", "--sample", help="Sample to run over (WWW, Wh, Signal)", type=str, default="Signal")

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
    "tot" : hist.GetBinContent(1),
    "tot_unc" : hist.GetBinError(1),
    "high_FR" : hist.GetBinContent(2),
    "low_FR" : hist.GetBinContent(2),
    "f0p5r0p5" : hist.GetBinContent(9),
    "f2r2" : hist.GetBinContent(5),
    "alpha_up" : hist.GetBinContent(10),
    "alpha_dn" : hist.GetBinContent(11),
    "pdf_up" : hist.GetBinContent(12),
    "pdf_dn" : hist.GetBinContent(13)
  }

  for i in xrange(3,11):
    count = hist.GetBinContent(i)
    if count > variations[sr]["high_FR"]:
      print("high for %s came from %d" % (sr, i) )
      variations[sr]["high_FR"] = count
    elif count < variations[sr]["low_FR"]:
      variations[sr]["low_FR"] = count
      print("low for %s came from %d" % (sr, i) )

def printLatexTable(SRs, sample):
  outfile = open("outputs/MC_Variation_%s.tex" % sample, "w+")

  pretty_SR_names = {"2lepSSEE": "$ee$","2lepSSEMu": "$e \mu$","2lepSSMuMu":  "$\mu \mu$","3lep_0SFOS": "0SFOS","3lep_1SFOS": "1SFOS","3lep_2SFOS": "2SFOS"}
  var_names = ["tot", "high_FR", "low_FR", "alpha_up", "alpha_dn", "pdf_up", "pdf_dn"]
  pretty_var_names = {"tot": "Yield", "high_FR": "F\\&R High", "low_FR": "F\\&R Low", "f2r2": "F2R2", "f0p5r0p5": "F$0.5$R$0.5$", "alpha_up": "$\\alpha_s$ Low", "alpha_dn": "$\\alpha_s$ High", "pdf_up": "PDF Low", "pdf_dn": "PDF High"}
  
  outfile.write("\\begin{table}[ht!]\n")
  outfile.write("\\begin{center}\n")
  outfile.write("\\begin{tabular}{l"+"|c"*len(var_names)+"}\n")

  header="Variation "
  for v in var_names:
    header+="& %s " % pretty_var_names[v]
  header+="\\\\ \\hline\n"

  outfile.write(header)

  ##==============
  ## Variations
  ##==============
  for sr in SRs:
    line = "%s" % pretty_SR_names[sr]
    for v in var_names:
      if v == "tot":
        line+=" & $%0.2f \pm %0.2f$" % (variations[sr][tot], variations[sr]["tot_unc"])
      else:  
        line+=" & $%0.2f (%0.2f)$" % (variations[sr][v], variations[sr][v]/variations[sr]["tot"])
    line+= "\\\\ \n"
    outfile.write(line)

  outfile.write("\\end{tabular} \n")
  outfile.write("\\end{center} \n")
  outfile.write("\\end{table} \n")

def main():
  SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/Sync/SRs/"
  sample = args.sample

  for signal_region in SRs:
    fillYields("%s/%s/%s.root" % (base_hists_path, signal_region, sample), signal_region)

  printLatexTable(SRs, sample)

if __name__ == "__main__":
  main()