#!/usr/bin/python
import sys, re


def copyLines(input_file, output_file):
  """ Copies lines in input_file to output_file unless it comes upon an include statement, 
  in which case it will call itself with the same with the included file as the new input file.
  Through this method, it will recursively replace <include > statments in the original file with
  the contents of the included files in place. """

  for line in input_file:
    line = re.sub( '\s+', ' ', line ).strip()
    #For now the only preprocessing is the <include> command
    if ( ("< include" in line) or ("<include" in line) ):
      print(line)
      inc_path=line[line.find("include")+7:line.find(">")].strip()
      print(inc_path)
      include_file=open(inc_path, 'r')
      copyLines(include_file, output_file)
    else:
      output_file.write(line)
      output_file.write('\n')

def preprocess():
  """This program preprocesses the configs in config parser"""
  if (len(sys.argv) < 2):
    print("There was no config file given\n");
    return

  input_path=sys.argv[1]
  output_path=input_path.replace(".conf", ".conf_tmp")

  input_file = open(input_path, 'r');
  output_file = open(output_path, 'w');

  copyLines(input_file, output_file)

if __name__ == "__main__":
  preprocess();