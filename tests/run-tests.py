#!/usr/bin/env python3
import os
import sys
import subprocess

SCRIPT_NAME = 'allele-counts.py'
DATASETS = [
  'artificial',
  'artificial-samples',
  'artificial-nofilt',
  'real',
  'real-mit',
  'real-mit-s',
  'real-nofilt',
]
IN_EXT  = '.vcf.in'
OUT_EXT = '.csv.out'
ARGS_KEY = '##comment="ARGS='

XML = {
  'tests_start':'  <tests>',
  'test_start': '    <test>',
  'input':      '      <param name="input" value="tests/%s" />',
  'param':      '      <param name="%s" value="%s" />',
  'output':     '      <output name="output" file="tests/%s" />',
  'test_end':   '    </test>',
  'tests_end':  '  </tests>',
}
PARAMS = {
  '-f':'freq',
  '-c':'covg',
  '-H':'header',
  '-s':'stranded',
  '-n':'nofilt',
  '-r':'seed',
}
PARAM_ARG = {
  '-f':True,
  '-c':True,
  '-H':False,
  '-s':False,
  '-n':False,
  '-r':True,
}

def main():

  do_print_xml = False
  if len(sys.argv) > 1:
    if sys.argv[1] == '-x':
      do_print_xml = True
    else:
      sys.stderr.write("Error: unrecognized option '"+sys.argv[1]+"'\n")
      sys.exit(1)

  test_dir = os.path.dirname(os.path.realpath(__file__))
  script_dir = os.path.relpath(os.path.dirname(test_dir))
  test_dir = os.path.relpath(test_dir)

  if do_print_xml:
    print(XML.get('tests_start'))

  for dataset in DATASETS:
    infile  = os.path.join(test_dir, dataset+IN_EXT)
    outfile = os.path.join(test_dir, dataset+OUT_EXT)

    if not os.path.exists(infile):
      sys.stderr.write("Error: file not found: "+infile+"\n")
      continue
    if not os.path.exists(outfile):
      sys.stderr.write("Error: file not found: "+outfile+"\n")
      continue

    options = read_options(infile)
    if do_print_xml:
      print_xml(infile, outfile, options, XML, PARAMS, PARAM_ARG)
    else:
      run_tests(infile, outfile, options, script_dir)

  if do_print_xml:
    print(XML.get('tests_end'))


def run_tests(infile, outfile, options, script_dir):
  script_cmd = os.path.join(script_dir, SCRIPT_NAME)+' '+options+' -i '+infile
  bash_cmd = 'diff '+outfile+' <('+script_cmd+')'
  print(script_cmd)
  subprocess.call(['bash', '-c', bash_cmd])


def print_xml(infile, outfile, options_str, xml, params, param_arg):
  infile = os.path.basename(infile)
  outfile = os.path.basename(outfile)

  options = options_str.split()  # on whitespace

  print(xml.get('test_start'))
  print(xml.get('input') % infile)

  # read in options one at a time, print <param> line
  i = 0
  while i < len(options):
    opt = options[i]
    if opt not in params or opt not in param_arg:
      sys.stderr.write("Error: unknown option '"+opt+"' in ARGS list in file "+infile+"\n")
      sys.exit(1)
    # takes argument
    if param_arg[opt]:
      i+=1
      arg = options[i]
      print(xml.get('param') % (params[opt], arg))
    # no argument (boolean)
    else:
      print(xml.get('param') % (params[opt], 'true'))
    i+=1

  print(xml.get('output') % outfile)
  print(xml.get('test_end'))


def read_options(infile):
  with open(infile, 'r') as infilehandle:
    for line in infilehandle:
      line.strip()
      if ARGS_KEY == line[:len(ARGS_KEY)]:
        return line[len(ARGS_KEY):-2]
  return ''


if __name__ == '__main__':
  main()