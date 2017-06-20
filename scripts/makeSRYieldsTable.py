import ROOT as r
import argparse, os, sys, root_utils

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