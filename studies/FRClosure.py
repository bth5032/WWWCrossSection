#!/usr/bin/env python

"""Plots and outputs the rate of predicting BG samples from Loose/Tight into Tight/Tight vs. QCD"""

import argparse, os, sys, math, array, ctypes, itertools
import ROOT as r
r.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-s", "--study_dir", help="The config directory name in FRClosure, e.g. LooseIso", type=str, default="LooseIso")
parser.add_argument("-p", "--pt_bins", help="Choose Pt bins for the study, ex: [10, 15, 20, 25, 30, 35, 50, 120]", type=str, default="[10, 15, 20, 25, 30, 35, 50, 120]")
parser.add_argument("-e", "--eta_bins", help="Choose Eta bins for the study, ex: [0, 1.2, 2.1, 2.4]", type=str, default="[0, 0.8, 1.479, 2.5]")
parser.add_argument("--txt", help="Print table in txt format", action="store_true")

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
parser.add_argument("--zh", help="Use ZH sample", action="store_true")
parser.add_argument("--www", help="Use WWW sample", action="store_true")
parser.add_argument("--wh", help="Use WH sample", action="store_true")

parser.add_argument("-l", "--lepton", help="Choose which lepton is used to select the pt bin. 1 for leading, 2 for subleading, -1 for trailing.", type=int, default=-1)

parser.add_argument("-u", "--usage", help="Print help message and quit", action="store_true")

args=parser.parse_args()

#SRs = ["2lepSS", "2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
SRs = ["2lepSS", "2lepSSEE","2lepSSEMu","2lepSSMuMu"]
pretty_SR_names = {"2lepSS": "SS",
"2lepSSEE": "SSee",
"2lepSSEMu": "SSe$\mu$",
"2lepSSMuMu": "SS$\mu\mu$",
"3lep_0SFOS": "0SFOS",
"3lep_1SFOS": "1SFOS",
"3lep_2SFOS": "2SFOS"}
lepton=-1

base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/"
base_plot_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/"

if (args.usage):
  parser.print_help()
  exit(0)

def err_ratio(num, den, Dnum, Dden):
  """Returns the 1 sigma on a ratio"""
  return abs((float(den)*Dnum - float(num)*Dden)/(den*den))

def getCell(nTight, nLoose, DnTight, DnLoose):
  """Returns (fr, Dfr) for the numerator and denomenator"""
  try:
    return (nTight, DnTight, nLoose, DnLoose,float(nTight)/nLoose, err_ratio(nTight, nLoose, DnTight, DnLoose))
  except Exception as e:
    #print(e)
    return (nTight, DnTight, nLoose, DnLoose, -1, -1)

def parsePtBins(bin_string):
  """Makes a list of ints from a string of the form [x0,x1,x2,...xN]"""
  s_bins = bin_string[1:-1].split(',')
  bins = []
  for s in s_bins:
    bins.append(float(s))
  if (bins[-1] < 6000):
    bins.append(6001.0)
  return bins

def parseEtaBins(bin_string):
  """Makes a list of ints from a string of the form [x0,x1,x2,...xN]"""
  if bin_string == "muon":
    return [0, 1.2, 2.1, 2.4]
  elif bin_string == "electron":
    return [0, 0.8, 1.479, 2.5]

  s_bins = bin_string[1:-1].split(',')
  bins = []
  for s in s_bins:
    bins.append(float(s))
  if (bins[-1] < 2.4):
    bins.append(2.4)
  return bins

def makePlots(bg_t, bg_l, qcd_t, qcd_l, pt_bins, sr):
  """Plots the QCD and BG stack as a function of the pt bins"""
  canvas=r.TCanvas("c_%s" % sr,"", 2000, 2000)
  canvas.cd()

  bg=bg_t.Clone("stack_%s" % sr)
  bg.Divide(bg_l)

  qcd=qcd_t.Clone("qcd_%s" % sr)
  qcd.Divide(qcd_l)

  pt_b_array = array.array('d', [float(a) for a in pt_bins])

  qcd = qcd.Rebin(len(pt_bins)-1,"qcd_rebinned_%s" % sr,  pt_b_array)
  bg = bg.Rebin(len(pt_bins)-1,"bg_rebinned_%s" % sr,  pt_b_array)

  bg.SetFillColor(r.kBlue)
  bg.Draw("HISTE")
  qcd.Draw("same")

  canvas.SaveAs("%s/%s.png" % (base_plot_path, sr))
  canvas.SaveAs("%s/%s.pdf" % (base_plot_path, sr)) 

def PrintSRTable(yields, pt_bins, eta_bins, latex):
  """Prints the yeilds tables and make plots"""
  print(yields)
  print("Using Config: %s" % args.study_dir)
  print("\\begin{table}[ht!]")
  print("\\begin{center}")
  head="\\begin{tabular}{l|l"
  title = "Signal Region & Source "
  for i in xrange(len(pt_bins) - 1):
    title+=" & %0.0f - %0.0f " % (pt_bins[i], pt_bins[i+1])
    head+="|c"
  title+= " \\\\ \\hline"
  head+="}"
  print(head)
  print(title)
  
  for sr in SRs:
    for j in xrange(len(eta_bins) - 1):
      row = "%s $(\eta \\in [%f,%f))$ &  BGs " % (pretty_SR_names[sr], eta_bins[j], eta_bins[j+1])
      for i in xrange(len(pt_bins) - 1):
        row += " & $\\frac{%0.2f \pm %0.2f}{%0.2f \pm %0.2f} (%0.2f \pm %0.2f)$ " % getCell(yields[sr]["bg_t"]["%f_%f" % (pt_bins[i], eta_bins[j])], yields[sr]["bg_l"]["%f_%f" % (pt_bins[i], eta_bins[j])], yields[sr]["bg_t_unc"]["%f_%f" % (pt_bins[i], eta_bins[j])], yields[sr]["bg_l_unc"]["%f_%f" % (pt_bins[i], eta_bins[j])])
      row+= " \\\\"
      print(row)

    #makePlots(yields[sr]["bg_sum_t"], yields[sr]["bg_sum_l"], yields[sr]["h_qcd_t"], yields[sr]["h_qcd_l"], pt_bins, sr)

  print("\\end{tabular}")
  print("\\end{center}")
  print("\\end{table}")

def getYieldsFromSample(hist_loc, SR, pt_bins, eta_bins, sample):
  """returns a tuple of lists ([Tight_bin1, ..., Tight_binN], [Loose_bin1, ..., Loose_binN]) for the number of events in"""
  f = r.TFile(hist_loc)

  if lepton == 1:
    h_tight_pteta = f.Get("tight_l1pteta").Clone("tight_pteta_%s_%s" % (SR, sample))
    h_loose_pteta = f.Get("loose_l1pteta").Clone("tight_pteta_%s_%s" % (SR, sample))
  elif lepton == 2:
    h_tight_pteta = f.Get("tight_l2pteta").Clone("tight_pteta_%s_%s" % (SR, sample))
    h_loose_pteta = f.Get("loose_l2pteta").Clone("tight_pteta_%s_%s" % (SR, sample))
  else:
    if("2lep" in SR):
      h_tight_pteta = f.Get("tight_l2pteta").Clone("tight_pteta_%s_%s" % (SR, sample))
      h_loose_pteta = f.Get("loose_l2pteta").Clone("tight_pteta_%s_%s" % (SR, sample))
    else:
      h_tight_pteta = f.Get("tight_l3pteta").Clone("tight_pteta_%s_%s" % (SR, sample))
      h_loose_pteta = f.Get("loose_l3pteta").Clone("tight_pteta_%s_%s" % (SR, sample))

  h_tight_pteta.SetDirectory(0)
  h_loose_pteta.SetDirectory(0)

  t = {}
  t_unc = {}
  l = {}
  l_unc = {}

  #loop over all pt intervals
  for i in xrange(len(pt_bins) - 1):
    for j in xrange(len(eta_bins) - 1):
      pt_low = h_tight_pteta.GetXaxis().FindBin(pt_bins[i])
      pt_high = h_tight_pteta.GetXaxis().FindBin(pt_bins[i+1]-0.0001)
      eta_low = h_tight_pteta.GetYaxis().FindBin(eta_bins[j])
      eta_high = h_tight_pteta.GetYaxis().FindBin(eta_bins[j+1]-.00001)
      bin_string = "%f_%f" % (pt_bins[i], eta_bins[j])

      t_unc_=r.Double()
      t_=h_tight_pteta.IntegralAndError(pt_low, pt_high, eta_low, eta_high, t_unc_)

      l_unc_=r.Double()
      l_=h_loose_pteta.IntegralAndError(pt_low, pt_high, eta_low, eta_high, l_unc_)
      
      t[bin_string] = t_
      t_unc[bin_string] = t_unc_
      l[bin_string] = l_
      l_unc[bin_string] = l_unc_

  print("Returning for %s from %s: tight=%s+/-%s, loose=%s+/-%s" % (sample, hist_loc, t, t_unc, l, l_unc))

  return (t, l, t_unc, l_unc, h_tight_pteta, h_loose_pteta)

def getCombinedYields(samples, pt_bins, eta_bins):
  """Constructs and returns the yields dictionary used in PrintTable. Goes through each SR and adds the yields for each sample in that SR and organizes them in the dict."""

  yields = {}
  #bg_stack_t = r.THStack("bg_stack_tight", "Background Stack for Tight Events")
  #bg_stack_l = r.THStack("bg_stack_loose", "Background Stack for Loose Events")

  for sr in SRs:
    t = {}
    l = {}
    t_unc = {}
    l_unc = {}
    
    for i in xrange(len(pt_bins) - 1):
      for j in xrange(len(eta_bins) -1):
        t["%f_%f" % (pt_bins[i], eta_bins[j])] = 0
        l["%f_%f" % (pt_bins[i], eta_bins[j])] = 0
        t_unc["%f_%f" % (pt_bins[i], eta_bins[j])] = 0
        l_unc["%f_%f" % (pt_bins[i], eta_bins[j])] = 0
      
    for s in samples:
      t_, l_, t_unc_, l_unc_, h_tight_pt, h_loose_pt = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr, pt_bins, eta_bins, s)
      #print("integral of thing %f" % h_tight_pt.Integral())
      for i in xrange(len(pt_bins) - 1):
        for j in xrange(len(eta_bins) -1):
          bin_string = "%f_%f" % (pt_bins[i], eta_bins[j])
          t[bin_string] += t_[bin_string]
          l[bin_string] += l_[bin_string]
          t_unc[bin_string] += t_unc_[bin_string]*t_unc_[bin_string]
          l_unc[bin_string] += l_unc_[bin_string]*l_unc_[bin_string]
          #bg_stack_t.Add(h_tight_pt)
          #bg_stack_l.Add(h_loose_pt)
          #if "bg_sum_t" in locals():
          #  bg_sum_t.Add(h_tight_pt)
          #  bg_sum_l.Add(h_loose_pt)
          #else:
          #  bg_sum_t = h_tight_pt.Clone("bg_sum_t_%s" % sr)
          #  bg_sum_t.SetDirectory(0)
          #  bg_sum_l = h_loose_pt.Clone("bg_sum_l_%s" % sr)
          #  bg_sum_l.SetDirectory(0)

      for i in xrange(len(pt_bins) - 1):
        for j in xrange(len(eta_bins) -1):
          bin_string = "%f_%f" % (pt_bins[i], eta_bins[j])
          t_unc[bin_string] = math.sqrt(t_unc[bin_string])
          l_unc[bin_string] = math.sqrt(l_unc[bin_string])

    #Get numbers for QCD as well
    #qcd_t, qcd_l, qcd_t_unc, qcd_l_unc, h_qcd_t, h_qcd_l = getYieldsFromSample(base_hists_path+sr+"/QCD.root", sr, pt_bins, "QCD")

    #yields[sr] = {"bg_t": t, "bg_l": l, "bg_t_unc": t_unc, "bg_l_unc": l_unc, "qcd_t": qcd_t, "qcd_l": qcd_l, "qcd_t_unc": qcd_t_unc, "qcd_l_unc": qcd_l_unc, "bg_stack_t": bg_stack_t, "bg_stack_l": bg_stack_l, "bg_sum_t": bg_sum_t, "bg_sum_l": bg_sum_l, "h_qcd_t": h_qcd_t, "h_qcd_l": h_qcd_l}
    yields[sr] = {"bg_t": t, "bg_l": l, "bg_t_unc": t_unc, "bg_l_unc": l_unc}

  return yields

def main():
  samples = []
  latex=True
  global base_hists_path, base_plot_path

  if (args.ttbar_dilep or args.ttbar_1lep or args.dy or args.wz or args.ww or args.zz or args.vvv or args.wjets or args.ttv or args.singletop or args.zh or args.www or args.wh):
    #Asking for specific samples
    #Signal
    if (args.www):
      samples.append("WWW")
    if (args.wh):
      samples.append("WH")

    #BG
    if (args.ttbar_dilep):
      samples.append("TTBar2l")
    if (args.wz):
      samples.append("WZ")
    if (args.zz):
      samples.append("ZZ")
    if (args.ww):
      samples.append("WW")
    if (args.ttv):
      samples.append("TTV")    
    if (args.vvv):
      samples.append("VVV")

    #FR BG
    if (args.ttbar_1lep):
      samples.append("TTBar1l")
    if (args.singletop):
      samples.append("SingleTop")
    if (args.wjets):
      samples.append("WJets")
    if (args.dy):
      samples.append("DY")
    if (args.zh):
      samples.append("ZH")
  else:
    #Asking for default samples
    samples = ["WJets", "TTBar1l"]


  if (args.txt):
    latex=False


  print("Going to use %s to make table from %s..." % (samples, args.study_dir))
  
  pt_bins = parsePtBins(args.pt_bins)
  eta_bins = parseEtaBins(args.eta_bins)
  lepton=args.lepton
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/%s/" % args.study_dir
  base_plot_path = "/home/users/bhashemi/public_html/WWWCrossSection/FRClosure/%s/" % args.study_dir

  yields = getCombinedYields(samples, pt_bins, eta_bins)
  #print(yields)
  PrintSRTable(yields, pt_bins, eta_bins, latex)

if __name__ == "__main__":
  main()