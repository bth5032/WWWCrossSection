import argparse, os, sys, root_utils, math

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-s", "--study_dir", help="The config directory name in FRStudy, e.g. LooseIso", type=str, default="LooseIso")
parser.add_argument("-b", "--bins", help="Choose Pt bins for the study, ex: [0,10,20,30]", type=str, default="[0,25,50,75,100,150]")
parser.add_argument("--txt", help="Print table in txt format", action="store_true")
parser.add_argument("-f", "--fake_rate", help="give fake rate either as a single number or list, ex: 0.1 or [0.1,0.2,0.3,0.2] ", type=str, default="0.1")

parser.add_argument("--all", help="Use all histograms", action="store_true")
parser.add_argument("--frbg", help="Use standard FR BG samples", action="store_true")
parser.add_argument("--signal", help="Use signal samples", action="store_true")
parser.add_argument("--bg", help="Use all BG samples", action="store_true")

parser.add_argument("--ttbar_dilep", help="Use TTBar to dilepton sample", action="store_true")
parser.add_argument("--ttbar_1lep", help="Use TTBar to single lepton samples", action="store_true")
parser.add_argument("--dy", help="Use DY sample", action="store_true")
parser.add_argument("--wz", help="Use WZ sample", action="store_true")
parser.add_argument("--wh", help="Use WH sample", action="store_true")
parser.add_argument("--ww", help="Use WW sample", action="store_true")
parser.add_argument("--zz", help="Use ZZ sample", action="store_true")
parser.add_argument("--vvv", help="Use VVV sample", action="store_true")
parser.add_argument("--wjets", help="Use W+Jets sample", action="store_true")
parser.add_argument("--ttv", help="Use TTV sample", action="store_true")
parser.add_argument("--singletop", help="Use TTV sample", action="store_true")
parser.add_argument("--other", help="Use Other (VH) sample", action="store_true")
parser.add_argument("--www", help="Use WWW sample", action="store_true")

parser.add_argument("-u", "--usage", help="Print help message and quit", action="store_true")

args=parser.parse_args()

SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
pretty_SR_names = {"2lepSSEE": "SSee",
"2lepSSEMu":  "SSe$\mu$  ",
"2lepSSMuMu": "SS$\mu\mu$",
"3lep_0SFOS": "0SFOS     ",
"3lep_1SFOS": "1SFOS     ",
"3lep_2SFOS": "2SFOS     "}

if (args.usage):
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

def parseFRBins(fr_string, pt_bins):
  """Takes in the pt_bins list and the string given on the command line for the fake rates. Checks whether FR is a single number or PT dependent and ensures the proper number of FRs are given if it's a list."""
  l = []
  
  try: #Check for single number
    single = float(fr_string)
    for a in xrange(len(pt_bins) - 1):
      l.append(single)
  except ValueError: #The string must be a list
    fr_bins = fr_string[1:-1].split(',')
    l = []
    for fr in fr_bins:
      l.append(float(fr))
  
  if (len(l) != (len(pt_bins) - 1) ): #Check the length of the FR list is equal to the number of pt_bin intervals.
    print("The length of the fake rate list is not equal to the number of fake rates needed. FR list: %s, pt_bins: %s" % (l, pt_bins))
    exit(1)

  return l

def PrintSRTable(yields, study_dir, pt_bins, fr_bins, latex):
  """Prints the yeilds tables for tight_fake/loose_fake and tight/loose (inclusive)"""
  print("Using Config: %s" % study_dir)
  print("\\begin{table}[ht!]")
  print("\\begin{center}")
  head="\\begin{tabular}{l"
  title = "$p_{T}$ bin: "
  for i in xrange(len(pt_bins) - 1):
    title+=" & %0.0f - %0.0f [GeV]" % (pt_bins[i], pt_bins[i+1])
    head+="|c"
  title+= " \\\\ \\hline"
  head+="}"
  print(head)
  print(title)
  
  for sr in SRs:
    row_loose = "%s (loose):     " % pretty_SR_names[sr]
    row_tight = "%s (tight):     " % pretty_SR_names[sr]
    row_fr    = "%s (fake rate): " % pretty_SR_names[sr]
    for i in xrange(len(pt_bins) - 1):
      fr, fr_unc = getCell(yields[sr]["tf"][i], yields[sr]["lf"][i], yields[sr]["tf_unc"][i], yields[sr]["lf_unc"][i])
      row_loose += " & %0.2f $\pm$ %0.2f (%0.2f $\pm$ %0.2f)" % (yields[sr]["l"][i], yields[sr]["l_unc"][i], yields[sr]["lf"][i], yields[sr]["lf_unc"][i])
      row_tight += " & %0.2f $\pm$ %0.2f (%0.2f $\pm$ %0.2f)" % (yields[sr]["t"][i], yields[sr]["t_unc"][i], yields[sr]["tf"][i], yields[sr]["tf_unc"][i])
      row_fr    += " & %0.2f $\pm$ %0.2f " % (fr, fr_unc)

    row_loose += "\\\\"
    row_tight += "\\\\"
    row_fr += "\\\\ \\hline"
    print(row_loose)
    print(row_tight)
    print(row_fr)

  print("\\end{tabular}")
  print("\\end{center}")
  print("\\end{table}")

def getYieldsFromSample(hist_loc, SR, pt_bins):
  """Takes in Signal region and the location to a histogram, looks up the proper bins for what is called real tight, real loose, fake tight, fake loose, and returns the counts for the sample in a tuple (rt, rl, ft, fl)"""
  f = r.TFile(hist_loc)
  if (SR == "3lep_0SFOS" or SR == "3lep_1SFOS" or SR == "3lep_2SFOS"):
    h_loose_pt = f.Get("loose_lep3pt").Clone("h_loose_pt_%s" % SR)
    h_tight_pt = f.Get("tight_lep3pt").Clone("h_tight_pt_%s" % SR)
    h_tight_fake_pt = f.Get("tight_fake_lep3pt").Clone("h_tight_fake_pt_%s" % SR)
    h_loose_fake_pt = f.Get("loose_fake_lep3pt").Clone("h_loose_fake_pt_%s" % SR)
  else:
    h_loose_pt = f.Get("loose_lep2pt").Clone("h_loose_pt_%s" % SR)
    h_tight_pt = f.Get("tight_lep2pt").Clone("h_tight_pt_%s" % SR)
    h_tight_fake_pt = f.Get("tight_fake_lep2pt").Clone("h_tight_fake_pt_%s" % SR)
    h_loose_fake_pt = f.Get("loose_fake_lep2pt").Clone("h_loose_fake_pt_%s" % SR)

  t = []
  t_unc = []
  l = []
  l_unc = []
  tf = []
  tf_unc = []
  lf = []
  lf_unc = []

  #loop over all pt intervals
  for i in xrange(len(pt_bins) - 1):
    low = h_tight_pt.FindBin(pt_bins[i])
    high = h_tight_pt.FindBin(pt_bins[i+1])

    t_unc_=r.Double()
    t_=h_tight_pt.IntegralAndError(low, high, t_unc_)

    l_unc_=r.Double()
    l_=h_loose_pt.IntegralAndError(low, high, l_unc_)

    tf_unc_=r.Double()
    tf_=h_tight_fake_pt.IntegralAndError(low, high, tf_unc_)

    lf_unc_=r.Double()
    lf_=h_loose_fake_pt.IntegralAndError(low, high, lf_unc_)

    t.append(t_)
    t_unc.append(t_unc_)
    l.append(l_)
    l_unc.append(l_unc_)
    tf.append(tf_)
    tf_unc.append(tf_unc_)
    lf.append(lf_)
    lf_unc.append(lf_unc_)

  return (t, l, tf, lf, t_unc, l_unc, tf_unc, lf_unc)

def getCombinedYields(samples, study_dir, pt_bins):
  """Constructs and returns the yields dictionary used in PrintTable. Goes through each SR and adds the yields for each sample in that SR and organizes them in the dict."""
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRStudy/%s/" % study_dir

  yields = {}

  for sr in SRs:
    t = []
    l = []
    tf = []
    lf = []
    t_unc = []
    l_unc = []
    tf_unc = []
    lf_unc = []
    for i in xrange(len(pt_bins) - 1):
      t.append(0)
      l.append(0)
      tf.append(0)
      lf.append(0)
      t_unc.append(0)
      l_unc.append(0)
      tf_unc.append(0)
      lf_unc.append(0)

    for s in samples:
      t_, l_, tf_, lf_, t_unc_, l_unc_, tf_unc_, lf_unc_ = getYieldsFromSample(base_hists_path+sr+"/"+s+".root", sr, pt_bins)
      for i in xrange(len(pt_bins) - 1):
        t[i] += t_[i]
        l[i] += l_[i]
        tf[i] += tf_[i]
        lf[i] += lf_[i]
        t_unc[i] += t_unc_[i]*t_unc_[i]
        l_unc[i] += l_unc_[i]*l_unc_[i]
        tf_unc[i] += tf_unc_[i]*tf_unc_[i]
        lf_unc[i] += lf_unc_[i]*lf_unc_[i]

    for i in xrange(len(pt_bins) - 1):
      t_unc[i] = math.sqrt(t_unc[i])
      l_unc[i] = math.sqrt(l_unc[i])
      tf_unc[i] = math.sqrt(tf_unc[i])
      lf_unc[i] = math.sqrt(lf_unc[i])

    yields[sr] = {"t": t, "l": l, "tf": tf, "lf": lf, "t_unc": t_unc, "l_unc": l_unc, "tf_unc": tf_unc, "lf_unc": lf_unc}

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
  
  pt_bins = parsePtBins(args.bins)
  fr_bins = parseFRBins(args.fake_rate, pt_bins)

  yields = getCombinedYields(samples, args.study_dir, pt_bins)
  #print(yields)
  PrintSRTable(yields, args.study_dir, pt_bins, fr_bins, latex)

if __name__ == "__main__":
  main()