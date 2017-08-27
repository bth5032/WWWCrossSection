#!/usr/bin/env python

"""Plots and outputs the rate of predicting BG samples from Loose/Tight into Tight/Tight vs. QCD"""

import argparse, os, sys, root_utils, math

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-s", "--study_dir", help="The config directory name in FRClosure, e.g. LooseIso", type=str, default="LooseIso")
parser.add_argument("-b", "--bins", help="Choose Pt bins for the study, ex: [0,10,20,30]", type=str, default="[0,25,50,75,100,150]")
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

SRs = ["2lepSS", "2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
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

import ROOT as r

def err_ratio(num, den, Dnum, Dden):
  """Returns the 1 sigma on a ratio"""
  return abs((den*Dnum - num*Dden)/(den*den))

def getCell(nTight, nLoose, DnTight, DnLoose):
  """Returns (fr, Dfr) for the numerator and denomenator"""
  try:
    return (nTight, DnTight, nLoose, DnLoose,float(num)/den, err_ratio(num, den, Dnum, Dden))
  except:
    return (-1, -1, -1, -1, -1, -1)

def parsePtBins(bin_string):
  """Makes a list of ints from a string of the form [x0,x1,x2,...xN]"""
  s_bins = bin_string[1:-1].split(',')
  bins = []
  for s in s_bins:
    bins.append(float(s))
  if (bins[-1] < 6000):
    bins.append(6001.0)
  return bins

def makePlots(stack_t, stack_l, qcd_t, qcd_l, pt_bins, sr):
  """Plots the QCD and BG stack as a function of the pt bins"""
  a=r.TCanvas("c","", 2000, 2000)
  a.cd()

  stack=stack_t.Clone("stack_%s" % sr)
  stack.Divide(stack_l)

  qcd=qcd_t.Clone("qcd_%s" % sr)
  qcd.Divide(qcd_l)

  stack.Draw()
  qcd.Draw("same")

  a.SaveAs("%s_%s.png" % (base_plot_path, sr))
  a.SaveAs("%s_%s.pdf" % (base_plot_path, sr)) 

def PrintSRTable(yields, study_dir, pt_bins, latex):
  """Prints the yeilds tables and make plots"""
  print("Using Config: %s" % study_dir)
  print("\\begin{table}[ht!]")
  print("\\begin{center}")
  head="\\begin{tabular}{l"
  title = "Source "
  for i in xrange(len(pt_bins) - 1):
    title+=" & %0.0f - %0.0f " % (pt_bins[i], pt_bins[i+1])
    head+="|c"
  title+= " \\\\ \\hline"
  head+="}"
  print(head)
  print(title)
  
  for sr in SRs:
    row = "%s (BGs)        " % pretty_SR_names[sr]
    for i in xrange(len(pt_bins) - 1):
      row += " & \\frac{%0.2f $\pm$ %0.2f}{%0.2f $\pm$ %0.2f} (%0.2f $\pm$ %0.2f) " % getCell(yields[sr][bg_t][i], yields[sr][bg_l][i], yields[sr][bg_t_unc][i], yields[sr][bg_l_unc][i])
    row+= " \\\\"
    print(row)
    row = "%s (QCD)        " % pretty_SR_names[sr]
    for i in xrange(len(pt_bins) - 1):
      row += " & \\frac{%0.2f $\pm$ %0.2f}{%0.2f $\pm$ %0.2f} (%0.2f $\pm$ %0.2f) " % getCell(yields[sr][qcd_t][i], yields[sr][qcd_l][i], yields[sr][qcd_t_unc][i], yields[sr][qcd_l_unc][i])
    row+= " \\\\"
    print(row)

    makePlots(yields[sr][])

  print("\\end{tabular}")
  print("\\end{center}")
  print("\\end{table}")


def getYieldsFromSample(hist_loc, SR, pt_bins, sample):
  """returns a tuple of lists ([Tight_bin1, ..., Tight_binN], [Loose_bin1, ..., Loose_binN]) for the number of events in"""
  f = r.TFile(hist_loc)

  if lepton == 1:
    h_tight_pt = f.Get("tight_lep1pt").Clone("tight_pt_%s_%s" % (SR, sample))
    h_loose_pt = f.Get("tight_lep1pt").Clone("tight_pt_%s_%s" % (SR, sample))
  elif lepton == 2:
    h_tight_pt = f.Get("tight_lep2pt").Clone("tight_pt_%s_%s" % (SR, sample))
    h_loose_pt = f.Get("tight_lep2pt").Clone("tight_pt_%s_%s" % (SR, sample))
  else:
    if("2lep" in SR):
      h_tight_pt = f.Get("tight_lep2pt").Clone("tight_pt_%s_%s" % (SR, sample))
      h_loose_pt = f.Get("tight_lep2pt").Clone("tight_pt_%s_%s" % (SR, sample))
    else:
      h_tight_pt = f.Get("tight_lep3pt").Clone("tight_pt_%s_%s" % (SR, sample))
      h_loose_pt = f.Get("tight_lep3pt").Clone("tight_pt_%s_%s" % (SR, sample))

  t = []
  t_unc = []
  l = []
  l_unc = []

  #loop over all pt intervals
  for i in xrange(len(pt_bins) - 1):
    low = h_tight_pt.FindBin(pt_bins[i])
    high = h_tight_pt.FindBin(pt_bins[i+1])

    t_unc_=r.Double()
    t_=h_tight_pt.IntegralAndError(low, high, t_unc_)

    l_unc_=r.Double()
    l_=h_loose_pt.IntegralAndError(low, high, l_unc_)
    
    t.append(rt_)
    t_unc.append(rt_unc_)
    l.append(rl_)
    l_unc.append(rl_unc_)

  return (t, l, t_unc, l_unc, h_tight_pt, h_loose_pt)

def getCombinedYields(samples, pt_bins):
  """Constructs and returns the yields dictionary used in PrintTable. Goes through each SR and adds the yields for each sample in that SR and organizes them in the dict."""

  yields = {}
  bg_stack_t = r.THStack("bg_stack_tight", "Background Stack for Tight Events")
  bg_stack_l = r.THStack("bg_stack_loose", "Background Stack for Loose Events")

  for sr in SRs:
    t = []
    l = []
    t_unc = []
    l_unc = []
    
    qcd_t = []
    qcd_l = []
    qcd_t_unc = []
    qcd_l_unc = []
    
    for i in xrange(len(pt_bins) - 1):
      t.append(0)
      l.append(0)
      t_unc.append(0)
      l_unc.append(0)
      qcd_t.append(0)
      qcd_l.append(0)
      qcd_t_unc.append(0)
      qcd_l_unc.append(0)
      
    for s in samples:
      t_, l_, t_unc_, l_unc_, h_tight_pt, h_loose_pt = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr, pt_bins, s)
      for i in xrange(len(pt_bins) - 1):
        t[i] += t_[i]
        l[i] += l_[i]
        t_unc[i] += t_unc_[i]*t_unc_[i]
        l_unc[i] += l_unc_[i]*l_unc_[i]
      bg_stack_t.Add(h_tight_pt)
      bg_stack_l.Add(h_loose_pt)

    #Get numbers for QCD as well
    t_, l_, t_unc_, l_unc_, h_qcd_t, h_qcd_l = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr, pt_bins, "QCD")
    for i in xrange(len(pt_bins) - 1):
      qcd_t[i] += qcd_t_[i]
      qcd_l[i] += qcd_l_[i]
      qcd_t_unc[i] += qcd_t_unc_[i]*qcd_t_unc_[i]
      qcd_l_unc[i] += qcd_l_unc_[i]*qcd_l_unc_[i]

    for i in xrange(len(pt_bins) - 1):
      t_unc[i] = math.sqrt(t_unc[i])
      l_unc[i] = math.sqrt(l_unc[i])
      qcd_t_unc[i] = math.sqrt(qcd_t_unc[i])
      qcd_l_unc[i] = math.sqrt(qcd_l_unc[i])

    yields[sr] = {"bg_t": t, "bg_l": l, "bg_t_unc": t_unc, "bg_l_unc": l_unc, "qcd_t": qcd_t, "qcd_l": qcd_l, "qcd_t_unc": qcd_t_unc, "qcd_l_unc": qcd_l_unc, "h_bg_t": bg_stack_t, "h_bg_l": bg_stack_l, "h_qcd_t": h_qcd_t, "h_qcd_l": h_qcd_l}

  return yields

def main():
  samples = []
  latex=True


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


  print("Going to use %s to make table..." %samples)
  
  pt_bins = parsePtBins(args.bins)
  lepton=args.lepton
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/%s/" % study_dir
  base_plot_path = "/home/users/bhashemi/public_html/WWWCrossSection/%s/" % study_dir

  yields = getCombinedYields(samples, args.study_dir, pt_bins)
  #print(yields)
  PrintSRTable(yields, args.study_dir, pt_bins, latex)

if __name__ == "__main__":
  main()