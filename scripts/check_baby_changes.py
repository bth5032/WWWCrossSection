import ROOT as r
import argparse, os, sys, root_utils

hadoop_baby_path = "/hadoop/cms/store/user/bhashemi/AutoTwopler_babies/"

longnames = {
"ttbar_dilep"      : "TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1",
"ttbar_1lep"       : "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1",
"dy"               : "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2",
"wz"               : "WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1",
"www"              : "WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1"
}

def makePlotsFromTrees(plotdir,t1,t2,name1,name2):
  """Draws all branches from both trees while overlapping branches that exist in both trees so you can see how the distribution changes. Eventually would like to have this also print out event lists for each tree too"""
  print t1, t2

  verbose = False

  info1 = root_utils.get_branch_info(t1)
  info2 = root_utils.get_branch_info(t2)

  keys_1_not_2 = set(info1.keys()) - set(info2.keys())
  keys_2_not_1 = set(info2.keys()) - set(info1.keys())
  keys_both = set(info1.keys()) - keys_1_not_2

  print "Number of branches only in file 1:", len(keys_1_not_2)
  if verbose: print "\t"+"\n\t".join(sorted(map(lambda x: info1[x].get("alias",x), keys_1_not_2)))
  print "Number of branches only in file 2:", len(keys_2_not_1)
  if verbose: print "\t"+"\n\t".join(sorted(map(lambda x: info2[x].get("alias",x), keys_2_not_1)))
  print "Number of branches common to both:", len(keys_both)
  if verbose: print "\t"+"\n\t".join(sorted(map(lambda x: info1[x].get("alias",x), keys_both)))

  max_nevts = 100000
  c1 = r.TCanvas("c1")
  # for b in list(keys_both)[:10]:
  for b in list(keys_both):
    try:
        i1 = info1[b]
        i2 = info2[b]
        is_p4 = "LorentzVector" in i1.get("type","")
        alias = i1["alias"]

        todraw = alias + (".pt()" if is_p4 else "")
        t1.Draw("{0}>>h1_{1}".format(todraw,alias),"","goff",max_nevts)
        h1 = r.gDirectory.Get("h1_{0}".format(alias))

        h2 = h1.Clone("h2_{0}".format(alias))
        h2.Reset()

        t2.Draw("{0}>>h2_{1}".format(todraw,alias),"","goff",max_nevts)

        # print h1.GetEntries(), h2.GetEntries(), h1.GetMean(), h2.GetMean()

        h1.SetLineColor(r.kRed)
        h2.SetLineColor(r.kBlue)

        pchi = h1.Chi2Test(h2)
        print pchi, alias

        pchistr = ("%.2f" % (pchi)).replace(".","p")

        h1.SetTitle(h1.GetTitle() + (" #color[2]{[%s]} #color[4]{[%s]} -- p_{#chi} = %.2f" % (name1,name2,pchi)))


        h1.Draw("hist")
        h2.Draw("samepe")

        c1.SaveAs("{0}/both/h_{1}_{2}.png".format(plotdir,pchistr,alias))
    except:
        print("Error: Caught exception %s" % e)
        pass

  for b in list(keys_1_not_2):
    try:
        i1 = info1[b]
        is_p4 = "LorentzVector" in i1.get("type","")
        alias = i1["alias"]

        todraw = alias + (".pt()" if is_p4 else "")
        t1.Draw("{0}>>h1_{1}_v1".format(todraw,alias),"","goff",max_nevts)
        h1 = r.gDirectory.Get("h1_{0}_v1".format(alias))

        # print h1.GetEntries(), h2.GetEntries(), h1.GetMean(), h2.GetMean()

        h1.SetLineColor(r.kRed)

        h1.SetTitle(h1.GetTitle() + (" #color[2]{[%s]}" % name1))


        h1.Draw("hist")

        c1.SaveAs("{0}/only_{1}/h_{2}.png".format(plotdir,name1,alias))
    except:
        print("Error: Caught exception %s" % e)
        pass

  for b in list(keys_2_not_1):
    try:
        i2 = info2[b]
        is_p4 = "LorentzVector" in i2.get("type","")
        alias = i2["alias"]

        todraw = alias + (".pt()" if is_p4 else "")
        t2.Draw("{0}>>h2_{1}_v2".format(todraw,alias),"","goff",max_nevts)
        h2 = r.gDirectory.Get("h2_{0}_v2".format(alias))

        h2.SetLineColor(r.kBlue)

        h2.SetTitle(h1.GetTitle() + ("#color[4]{[%s]} " % name2))


        h2.Draw("hist")

        c1.SaveAs("{0}/only_{1}/h_{2}.png".format(plotdir,name2,alias))
    except Exception as e:
        print("Error: Caught exception %s" % e)
        pass

def compareForFiles(f1, f2):
  r.gROOT.SetBatch(1)
  f1_basename = f1.split("/")[-1][:-5] # get everything after last '/' and assume it's a .root
  f2_basename = f2.split("/")[-1][:-5] # get everything after last '/' and assume it's a .root

  print("Doing comparison for files %s and %s" % (f1_basename, f2_basename))

  treename = "t"
  plotdir = "~/public_html/WWWComparison/dump/"

  os.system("mkdir -p {0}".format(plotdir))
  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/dump/index.php")

  os.system("mkdir -p {0}/both/".format(plotdir))
  os.system("mkdir -p {0}/only_{1}/".format(plotdir, f1_basename))
  os.system("mkdir -p {0}/only_{1}/".format(plotdir, f2_basename))
  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/dump/both/index.php")
  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/dump/only_%s/index.php" % f1_basename)
  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/dump/only_%s/index.php" % f2_basename)

  t1 = r.TChain(treename)
  t2 = r.TChain(treename)

  t1.Add(f1)
  t2.Add(f2)

  makePlotsFromTrees(plotdir,t1,t2,f1_basename,f2_basename)

def compareForSample(sample_name, v1, v2, baby_path=hadoop_baby_path):
  """Prints a list of which events are different between trees, also prints out all the branches"""
  r.gROOT.SetBatch(1)

  print("Doing comparison for %s" % sample_name)

  treename = "t"
  fnames1 = "%sZMET_WWW_%s/%s/output/output_*" % (baby_path, v1, longnames[sample_name])
  fnames2 = "%sZMET_WWW_%s/%s/output/output_*" % (baby_path, v2, longnames[sample_name])
  plotdir = "~/public_html/WWWComparison/%s_VS_%s/%s/" % (v1,v2, sample_name)

  os.system("mkdir -p {0}".format(plotdir))
  #os.system("cp scripts/index_nick.php ~/public_html/WWWComparison/%s_VS_%s/index.php" % (v1,v2))
  #os.system("cp scripts/index_nick.php ~/public_html/WWWComparison/%s_VS_%s/%s/index.php" % (v1,v2,sample_name))
  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/%s_VS_%s/index.php" % (v1,v2))
  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/%s_VS_%s/%s/index.php" % (v1,v2,sample_name))

  os.system("mkdir -p {0}/both/".format(plotdir))
  os.system("mkdir -p {0}/only_{1}".format(plotdir,v1))
  os.system("mkdir -p {0}/only_{1}".format(plotdir,v2))

  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/%s_VS_%s/%s/both/index.php" % (v1,v2,sample_name))
  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/%s_VS_%s/%s/only_%s/index.php" % (v1,v2,sample_name,v1))
  os.system("cp ~/public_html/index.php ~/public_html/WWWComparison/%s_VS_%s/%s/only_%s/index.php" % (v1,v2,sample_name,v2))

  t1 = r.TChain(treename)
  t2 = r.TChain(treename)

  t1.Add(fnames1)
  t2.Add(fnames2)

  makePlotsFromTrees(plotdir,t1,t2,v1,v2)

def main():
  parser = argparse.ArgumentParser(add_help=False)

  parser.add_argument("old_version", help="Version of older babies to check against", nargs='?', type=str, default="v0.1.4")
  parser.add_argument("new_version", help="Version of newer babies to check", nargs='?', type=str, default="v0.1.5")
  parser.add_argument("--f1", help="Direct path to the first file", type=str)
  parser.add_argument("--f2", help="Direct path to the second file", type=str)
  parser.add_argument("--ttbar_dilep", help="Do check on TTBar to dilepton sample", action="store_true")
  parser.add_argument("--ttbar_1lep", help="Do check on TTBar to single lepton (from top) sample", action="store_true")
  parser.add_argument("--dy", help="Do check on DY m50 sample", action="store_true")
  parser.add_argument("--wz", help="Do check on WZ sample", action="store_true")
  parser.add_argument("--www", help="Do check on WWW sample", action="store_true")
  parser.add_argument("--help", help="Print help message and quit", action="store_true")

  args=parser.parse_args()

  if (args.f1 and args.f2):
    compareForFiles(args.f1,args.f2)

  elif (args.f1 or args.f2):
    print("You must give both file locations --f1 <path_to_old_baby> --f2 <path_to_new_baby>")
    exit(0)

  else:
    if (args.ttbar_dilep):
      compareForSample("ttbar_dilep", args.old_version, args.new_version)
    elif (args.ttbar_1lep):
      compareForSample("ttbar_1lep", args.old_version, args.new_version)
    elif (args.dy):
      compareForSample("dy", args.old_version, args.new_version)
    elif (args.wz):
      compareForSample("wz", args.old_version, args.new_version)
    elif (args.www):
      compareForSample("www", args.old_version, args.new_version)

if __name__ == "__main__":
  main()