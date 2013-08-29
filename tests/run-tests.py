#!/usr/bin/env python
import os
import sys
import subprocess

DATASETS = ['artificial', 'artificial-samples', 'real']
IN_EXT  = '.vcf.in'
OUT_EXT = '.csv.out'
ARGS_KEY = '##comment="ARGS='

def main():

  test_dir = os.path.dirname(os.path.relpath(sys.argv[0]))
  print test_dir

  for dataset in DATASETS:
    infile  = test_dir+os.sep+dataset+IN_EXT
    outfile = test_dir+os.sep+dataset+OUT_EXT

    if not os.path.exists(infile):
      sys.stderr.write("Error: file not found: "+infile)
      continue
    if not os.path.exists(outfile):
      sys.stderr.write("Error: file not found: "+outfile)
      continue

    options = read_options(infile)
    script_cmd = 'allele-counts.py '+options+' -i '+infile
    bash_cmd = 'diff '+outfile+' <('+script_cmd+')'
    # print infile+":"
    print script_cmd
    subprocess.call(['bash', '-c', bash_cmd])


def read_options(infile):
  with open(infile, 'r') as infilehandle:
    for line in infilehandle:
      line.strip()
      if ARGS_KEY == line[:len(ARGS_KEY)]:
        return line[len(ARGS_KEY):-2]
  return ''


if __name__ == '__main__':
  main()