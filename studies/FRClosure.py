#!/usr/bin/env python

"""Creates Histogram for the Fake Rate Clsoure Test.

Before this script is run, the configs in FRClosure/<Conf>/ need to be made...
This script will output new histograms in the FRClosure/<Conf>/Combined/. 
One with the Fake rate prediction called 'fakes.root' and one for each sample."""

import argparse, os, sys, math, array, ctypes, itertools
import ROOT as r
r.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument("-s", "--study_dir", help="The config directory name in FRClosure, e.g. Btag", type=str, default="Btag")
parser.add_argument("-c", "--closure_error", help="Set the closure error on the fakerate method", type=float, default=0.3)
parser.add_argument("-r", "--signal_regions", help="Choose SS, 3lep, or all signal regions", type=str, default="SS")
parser.add_argument("-u", "--usage", help="Print help message and quit", action="store_true")

args=parser.parse_args()

if (args.usage):
  parser.print_help()
  exit(0)

def err_ratio(num, den, Dnum, Dden):
  """Returns the 1 sigma on a ratio"""
  return abs((float(den)*Dnum - float(num)*Dden)/(den*den))

def frErr(epsilon, Depsilon):
  """Computes the fake rate and error from the number in Phil's Hist, which is NTight/(NLoose+NTight)"""
  fr = epsilon/float(1-epsilon)
  Dfr = Depsilon*(1/float(1-epsilon))*(1+fr)
  return (fr, Dfr)

def buildIntervals(pt_bins, eta_bins):
  """Returns all intervals in pt and eta bins"""
  pt_intervals = [(pt_bins[i], pt_bins[i+1] - 0.001) for i in xrange(0, len(pt_bins) - 1) ]
  eta_intervals = [(eta_bins[i], eta_bins[i+1] - 0.001) for i in xrange(0, len(eta_bins) - 1) ]
  return [(i[0][0], i[0][1], i[1][0], i[1][1]) for i in itertools.product(pt_intervals, eta_intervals)]
  
def getYield(hist, pt_bins, eta_bins, h_fr):
  """Computes the number of objects in the pt/eta bins and applies the FR given with the error."""

  cv = 0
  dev = 0

  err = r.Double()
  for pt_low, pt_high, eta_low, eta_high in buildIntervals(pt_bins, eta_bins):
    x1=hist.GetXaxis().FindBin(pt_low)
    x2=hist.GetXaxis().FindBin(pt_high)
    y1=hist.GetYaxis().FindBin(eta_low)
    y2=hist.GetYaxis().FindBin(eta_high)
    #print(x1,x2,y1,y2)
    #print(pt_low, pt_high, eta_low, eta_high)
    y = hist.IntegralAndError(x1, x2, y1, y2, err)
    #print("Count in (pt, eta) in (%0.2f, %0.2f, %0.2f, %0.2f) = %0.2f +/- %0.2f" % (pt_low, pt_high, eta_low, eta_high, y, err) )
    
    #Get FR from Histogram
    epsilon_x = h_fr.GetXaxis().FindBin(pt_low)
    epsilon_y = h_fr.GetYaxis().FindBin(eta_low)
    epsilon = h_fr.GetBinContent(epsilon_x, epsilon_y)
    epsilon_err = h_fr.GetBinError(epsilon_x, epsilon_y)

    fr, fr_err = frErr(epsilon, ep)

    #print("Got fake rate: %0.2f +/ %0.2f" % (fr, fr_err))

    #print("Bin count is: %0.2f +/- %0.2f" % (fr*y, math.sqrt( (fr*err)**2 + (y*fr_err)**2 ) ) )
    yield_bin = fr*y
    err_bin = (fr*err)**2 + (y*fr_err)**2
    cv+=yield_bin
    dev += err_bin

  return (cv, math.sqrt(dev))
  
def makeHistos(pt_bins, eta_bins_m, eta_bins_e, sample_files, FR_Hist_path, signal_region):
  """Writes out the completed histogram for a list of sample files for the Fake Rate in that region."""
  fr_file = r.TFile(FR_Hist_path)
  h_m_FRs = fr_file.Get("muon_fakerate_conecorrpt_v_eta").Clone("h_m_FRs")
  h_e_FRs = fr_file.Get("elec_fakerate_conecorrpt_v_eta").Clone("h_e_FRs")

  output_dir = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/%s" % args.study_dir

  total_count = 0
  total_error = 0

  #Loop over samples in a particlar signal region
  for s_file in sample_files:
    tfile = r.TFile(s_file)
    s = s_file[s_file.rfind('/')+1:s_file.index(".root")] #extract sample name from filename
    
    h_loose_e = tfile.Get("loose_loose_pteta_e").Clone("h_loose_e_%s" % s)
    h_loose_m = tfile.Get("loose_loose_pteta_m").Clone("h_loose_m_%s" % s)

    num_e, err_e = getYield(h_loose_e, pt_bins, eta_bins_m, h_e_FRs)
    num_m, err_m = getYield(h_loose_m, pt_bins, eta_bins_m, h_m_FRs)

    print("count for %s: %0.2f +/- %0.2f " % (s, num_e+num_m, math.sqrt(err_e**2 + err_m**2)) )
    total_count += num_e+num_m
    total_error += err_e**2 + err_m**2

  print("Total Fake Count (plus closure error): %0.2f +/- %0.2f (%0.2f)" % (total_count, math.sqrt(total_error), math.sqrt(total_error+(total_count*args.closure_error)**2) ) ) 
  total_error += (total_count*args.closure_error)**2
  total_error = math.sqrt(total_error)


  outfile = r.TFile("%s/%s/Fakes.root" % (output_dir, signal_region), "RECREATE" )
  outfile.cd()
  h_fake_count = r.TH1D("weighted_count", "Weighted Count of Events", 4, 0, 4)
  h_fake_count.SetDirectory(0)
  h_fake_count.SetBinContent(2, total_count)
  h_fake_count.SetBinError(2, total_error)
  h_fake_count.Write()
  print("Wrote Total Fake Count: %0.2f +/- %0.2f" % (total_count, total_error) )
  outfile.Close()

def moveIntoCombined(hist_path, SRs):
  """Basically a fancy HADD just for the weighted_counts hist in the list of samples given. 
  Makes a histogram for each sample where the bins are the total counts in each of the signal reigons sent.
  Outputs the final hists into a directory called /Combined/ in the same config"""
  
  #Check if output path for hist exists
  if (not os.path.isdir(hist_path)):
    os.system("mkdir -p %s" % hist_path)

  samples = ["TTBar1l", "WJets", "DY", "Fakes"]
  pretty_SR_names = {"2lepSSEE": "ee", "2lepSSEMu":  "e #mu", "2lepSSMuMu": "#mu #mu", "3lep_0SFOS": "0SFOS", "3lep_1SFOS": "1SFOS", "3lep_2SFOS": "2SFOS"}
  latex_SR_names = {"2lepSSEE": "$ee$", "2lepSSEMu":  "$e \mu$", "2lepSSMuMu": "$\mu \mu$", "3lep_0SFOS": "0SFOS", "3lep_1SFOS": "1SFOS", "3lep_2SFOS": "2SFOS"}
  tex_dict = {}
  
  #Loop over samples, combining the 'weighted count' hist into a new hist in the /Combined/ directory
  for sample in samples:
    outfile = r.TFile("%s/Combined/%s.root" % (hist_path, sample), "RECREATE")
    outfile.cd()
    h_signal_count = r.TH1D("weighted_count", "Weighted Count of Events", 0, 0, 0)
    h_signal_count.SetDirectory(0)
    tex_dict[sample] = {}
    
    for sr in SRs:
      #Fill in an entry in the h_signal_count histogram for each signal region, label the x axis with the name of the SR.
      sample_file_in_sr = r.TFile("%s/%s/%s.root" % (hist_path, sr, sample), "r")
      sample_count = sample_file_in_sr.Get("weighted_count").Clone("weighted_count_%s_%s" % (sr, sample) )
      bin_num = h_signal_count.Fill(pretty_SR_names[sr], sample_count.GetBinContent(2))
      h_signal_count.SetBinError(bin_num, sample_count.GetBinError(2))
      tex_dict[sample][latex_SR_names[sr]] = {"cv": sample_count.GetBinContent(2), "unc": sample_count.GetBinError(1)}
    
    outfile.cd()
    h_signal_count.Write()
    outfile.Close()

  printLatexTable(tex_dict)

def printLatexTable(tex_dict):
  outfile = open("outputs/frhists_%s.tex" % args.study_dir, "w+")
  pretty_sample_names = {"WZ":"WZ", "WW":"WW", "ZZ":"ZZ", "TTBar2l":"TTBar $\\to$ 2lep", "DY":"DY", "WWW":"WWW", "Wh":"WH", "VVV":"VVV", "TTV":"TTV", "SingleTop":"SingleTop", "Other":"ZH", "Data":"Data", "Fakes":"Fakes"}

  if args.signal_regions == "SS":
    SRs = ["$ee$", "$e \mu$", "$\mu \mu$"]
  elif args.signal_regions == "3lep":
    SRs = ["0SFOS", "1SFOS", "2SFOS"]
  else:
    SRs = ["$ee$", "$e \mu$", "$\mu \mu$", "0SFOS", "1SFOS", "2SFOS"]

  line1 = tex_dict[tex_dict.keys()[0]]
  pred_count = {}

  outfile.write("\\begin{table}[ht!]\n")
  outfile.write("\\begin{center}\n")
  outfile.write("\\begin{tabular}{l"+"|c"*len(line1.keys())+"}\n")
  
  header="Sample  "
  for sr in SRs:
    header+="& %s " % sr
    pred_count[sr] = {"cv": 0, "unc": 0}
  header+="\\\\ \\hline\n"

  outfile.write(header)

  ##==============
  ## Predictions
  ##==============
  for sample in tex_dict.keys():
    if sample == "Data":
      continue
    line = "%s" % pretty_sample_names[sample]
    for sr in SRs:
      line+=" & $%0.2f \\pm %0.2f$" % (tex_dict[sample][sr]["cv"], tex_dict[sample][sr]["unc"])
      pred_count[sr]["cv"]+=tex_dict[sample][sr]["cv"]
      pred_count[sr]["unc"]+=tex_dict[sample][sr]["unc"]**2
    line+= "\\\\ \n"
    outfile.write(line)
  
  outfile.write("\\hline \n")

  ##==============
  ## Prediction Sum
  ##==============
  line = "Sum "
  for sr in SRs:
    line+=" & $%0.2f \\pm %0.2f$" % (pred_count[sr]["cv"], math.sqrt(pred_count[sr]["unc"]))
  line+= "\\\\ \n"
  outfile.write(line)

  ##==============
  ## Predictions
  ##==============
  sample="Data"
  line = "%s" % sample
  for sr in SRs:
    line+=" & $%0.2f \\pm %0.2f$" % (tex_dict[sample][sr]["cv"], tex_dict[sample][sr]["unc"])
  line+= "\\\\ \n"
  
  outfile.write(line)

  outfile.write("\\end{tabular} \n")
  outfile.write("\\end{center} \n")
  outfile.write("\\end{table} \n")

def printPlotMakerCommand(prediction_hist_path):
  conf_dir = "configs/%s/Combined" % prediction_hist_path[prediction_hist_path.find("WWWCrossSection_Hists")+22:]
  print("makeAllForDir %s plots" % conf_dir)

def main():
  samples = ["TTBar1l", "WJets"]
  
  if args.signal_regions == "SS":
    SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu"]
  elif args.signal_regions == "3lep":
    SRs = ["3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]
  else:
    SRs = ["2lepSSEE","2lepSSEMu","2lepSSMuMu","3lep_0SFOS","3lep_1SFOS","3lep_2SFOS"]

  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/"

  FR_Hist_path = "auxFiles/fakerate_pt_v_eta.root" 

  pt_bins = [10, 15, 20, 25, 30, 35, 50, 120, 6001]
  eta_bins_m = [0, 1.2, 2.1, 2.4, 3]
  eta_bins_e = [0, 0.8, 1.479, 2.5, 3]
  
  base_hists_path = "/nfs-7/userdata/bobak/WWWCrossSection_Hists/FRClosure/%s" % args.study_dir

  for signal_region in SRs:
    hist_paths= [ "%s/%s/%s.root" % (base_hists_path, signal_region, s) for s in samples]
    print("Signal Region: %s: " % signal_region)
    makeHistos(pt_bins, eta_bins_e, eta_bins_m, hist_paths, FR_Hist_path, signal_region)

  moveIntoCombined(base_hists_path, SRs)

if __name__ == "__main__":
  main()