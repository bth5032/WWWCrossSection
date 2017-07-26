import argparse, os, sys, root_utils, math

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-s", "--study_dir", help="The config directory name in FRStudy, e.g. LooseIso", type=str, default="LooseIso")
parser.add_argument("-b", "--bins", help="Choose Pt bins for the study, ex: [0,10,20,30]", type=str, default="[0,25,50,75,100,150]")
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

parser.add_argument("-h", "--help", help="Print help message and quit", action="store_true")

args=parser.parse_args()

SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
pretty_SR_names = {"2lepSSEE": "SSee",
"2lepSSEMu": "SSe$\mu$",
"2lepSSMuMu": "SS$\mu\mu$",
"3lep_0SFOS": "0SFOS",
"3lep_1SFOS": "1SFOS",
"3lep_2SFOS": "2SFOS"}

if (args.help):
  parser.print_help()
  exit(0)

import ROOT as r

def err_ratio(num, den, Dnum, Dden):
  """Returns the 1 sigma on a ratio"""
  return abs((den*Dnum - num*Dden)/(den*den))

def getCell(num, den, Dnum, Dden):
  """Returns (fr, Dfr) for the numerator and denomenator"""
  try:
    return (float(num)/den, err_ratio(num, den, Dnum, Dden))
  except:
    return (-1, -1)

def parsePtBins(bin_string):
  """Makes a list of ints from a string of the form [x0,x1,x2,...xN]"""
  s_bins = bin_string[1:-1].split(',')
  bins = []
  for s in s_bins:
    bins.append(float(s))
  if (bins[-1] < 6000):
    bins.append(6001.0)
  return bins

def PrintSRTable(yields, study_dir, pt_bins, latex):
  """Prints the yeilds tables for tight_fake/loose_fake and tight/loose (inclusive)"""
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
    row = "%s (fakes)         " % pretty_SR_names[sr]
    for i in xrange(len(pt_bins) - 1):
      row += " & %0.2f $\pm$ %0.2f " % (getCell(yields[sr]["t"][i], yields[sr]["l"][i], yields[sr]["t_unc"][i], yields[sr]["l_unc"][i]))
    row+= " \\\\"
    print(row)

  print("\\end{tabular}")
  print("\\end{center}")
  print("\\end{table}")

def getYieldsFromSample(hist_loc, SR, pt_bins):
  """Takes in Signal region and the location to a histogram, looks up the proper bins for what is called real tight, real loose, fake tight, fake loose, and returns the counts for the sample in a tuple (rt, rl, ft, fl)"""
  f = r.TFile(hist_loc)
  if (SR == "3lep_0SFOS" or SR == "3lep_1SFOS" or SR == "3lep_2SFOS"):
    h_loose_pt = f.Get("loose_lep3pt").Clone("h_loose_pt_%s" % SR)
    h_tight_pt = f.Get("tight_lep3pt").Clone("h_tight_pt_%s" % SR)
  else:
    h_loose_pt = f.Get("loose_lep2pt").Clone("h_loose_pt_%s" % SR)
    h_tight_pt = f.Get("tight_lep2pt").Clone("h_tight_pt_%s" % SR)

  t = []
  t_unc = []
  l = []
  l_unc = []

  #loop over all pt intervals
  for i in xrange(len(pt_bins) - 1):
    low = h_tight_real_pt.FindBin(pt_bins[i])
    high = h_tight_real_pt.FindBin(pt_bins[i+1])

    t_unc_=r.Double()
    t_=h_tight_real_pt.IntegralAndError(low, high, rt_unc_)

    l_unc_=r.Double()
    l_=h_loose_real_pt.IntegralAndError(low, high, rl_unc_)

    t.append(t_)
    t_unc.append(t_unc_)
    l.append(l_)
    l_unc.append(l_unc_)

  return (t, l, t_unc, l_unc)

def getCombinedYields(samples, study_dir, pt_bins):
  """Constructs and returns the yields dictionary used in PrintTable. Goes through each SR and adds the yields for each sample in that SR and organizes them in the dict."""
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRStudy/%s/" % study_dir

  yields = {}

  for sr in SRs:
    t = []
    l = []
    t_unc = []
    l_unc = []
    for i in xrange(len(pt_bins) - 1):
      t.append(0)
      l.append(0)
      t_unc.append(0)
      l_unc.append(0)

    for s in samples:
      t_, l_, t_unc_, l_unc_ = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr, pt_bins)
      for i in xrange(len(pt_bins) - 1):
        t[i] += t_[i]
        l[i] += l_[i]
        t_unc[i] += t_unc_[i]*t_unc_[i]
        l_unc[i] += l_unc_[i]*l_unc_[i]

    for i in xrange(len(pt_bins) - 1):
      t_unc[i] = math.sqrt(t_unc[i])
      l_unc[i] = math.sqrt(l_unc[i])

    yields[sr] = {"t": t, "l": l, "t_unc": t_unc, "l_unc": l_unc}

  return yields

def main():
  samples = []
  latex=True

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

  if (args.txt):
    latex=False


  print("Going to use %s to make table..." %samples)
  
  pt_bins = parsePtBins(args.bins)

  yields = getCombinedYields(samples, args.study_dir, pt_bins)
  #print(yields)
  PrintSRTable(yields, args.study_dir, pt_bins, latex)

if __name__ == "__main__":
  main()