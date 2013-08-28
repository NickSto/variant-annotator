#!/usr/bin/env python
import os
import subprocess

DATASETS = ['artificial', 'artificial-samples', 'real']
IN_EXT  = '.vcf.in'
OUT_EXT = '.csv.out'
ARGS_KEY = '##comment="ARGS='

def main():

  for dataset in DATASETS:
    infile  = dataset+IN_EXT
    outfile = dataset+OUT_EXT

    options = read_options(infile)
    bash_cmd = 'diff '+outfile+' <(allele-counts.py -i '+infile+' '+options+')'
    dash_cmd = "bash -c '"+bash_cmd+"'"
    # print infile+":"
    print bash_cmd
    subprocess.call(dash_cmd, shell=True)


def read_options(infile):
  with open(infile, 'r') as infilehandle:
    for line in infilehandle:
      line.strip()
      if ARGS_KEY == line[:len(ARGS_KEY)]:
        return line[len(ARGS_KEY):-2]
  return ''


if __name__ == '__main__':
  main()