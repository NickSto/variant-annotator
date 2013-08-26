#!/usr/bin/python
# This parses the output of Dan's "Naive Variant Caller" (previously, "BAM
# Coverage"). It was forked from the code of "bam-coverage.py".
#
# TODO:
# - use actual flagged options
# - option to include a header in the output
# - should it technically handle data lines that start with a '#'?
import os
import sys

COLUMNS = ['sample', 'chr', 'pos', 'A', 'C', 'G', 'T', 'coverage', 'alleles',
  'major', 'minor', 'freq'] #, 'bias']
BASES = ['A', 'C', 'G', 'T', 'N']
PRINT_HEADER_DEFAULT = True
FREQ_THRES_DEFAULT = 1 #percent
COVG_THRES_DEFAULT = 100

def main():

  args = len(sys.argv)
  if args == 1 or '-h' in sys.argv[1][0:3]:
    script_name = os.path.basename(sys.argv[0])
    print """USAGE:
  $ """+script_name+""" variants.vcf [frequency] [coverage] > alleles.csv
This will parse the VCF output of Dan's "Naive Variant Caller" (aka "BAM
Coverage") Galaxy tool. For each position reported, it counts the number of
reads of each base, determines the major allele, minor allele (second most
frequent variant), and number of alleles above a threshold.

Requirements:
input: The input VCF must report the variants per strand.
frequency: Alleles are only counted when they are above a default threshold
  of 1% frequency, but a different value can be given as the second argument.
coverage: Positions are required to have at least 100x coverage on each strand
  by default, but a different value can be given as the third argument.
  Positions below this coverage will not be reported.
strand bias: Both strands must show the same bases passing the frequency
  threshold (but not necessarily in the same order of frequency). If the
  position fails this test, the number of alleles is reported as 0."""
    sys.exit(0)

  print_header = PRINT_HEADER_DEFAULT
  freq_thres = FREQ_THRES_DEFAULT
  covg_thres = COVG_THRES_DEFAULT

  # get arguments
  if args > 1:
    filename = sys.argv[1]
  if args > 2:
    try:
      freq_thres = int(sys.argv[2])
    except ValueError, e:
      pass
  if args > 3:
    try:
      covg_thres = int(sys.argv[3])
    except ValueError, e:
      pass

  if print_header:
    print '#'+'\t'.join(COLUMNS)

  # main loop: process and print one line at a time
  sample_names = []
  with open(filename, 'r') as lines:
    for line in lines:
      line = line.rstrip('\r\n')

      # header lines
      if line[0] == '#':
        if line[0:6].upper() == '#CHROM':
          sample_names = line.split('\t')[9:]
        continue

      if not sample_names:
        fail("Error in input VCF: Data line encountered before header line.")

      site_data = read_site(line, sample_names)

      site_summary = summarize_site(site_data, sample_names, BASES, freq_thres,
        covg_thres)

      print_site(site_summary, COLUMNS)

      # return



  print HEADER
  for snp in snps:
    total = As = Ts = Gs = Cs = Ns = 0
    for variant in snp['variants']:
      As += variant.get('a', 0)
      Ts += variant.get('t', 0)
      Gs += variant.get('g', 0)
      Cs += variant.get('c', 0)
      Ns += variant.get('n', 0)
    total = As + Ts + Gs + Cs + Ns
    line = "\t".join(map(str, [total, As, Ts, Gs, Cs, Ns]))
    sys.stdout.write(snp['chr']+':'+snp['pos']+'\t')
    print line


def read_site(line, sample_names):
  """Read in a line, parse the variants into a data structure, and return it.
  The line should be actual site data, not a header line, so check beforehand.
  Notes:
  - The line is assumed to have been chomped."""
  
  site = {}
  fields = line.split('\t')

  if len(fields) < 9:
    fail("Error in input VCF: wrong number of fields in data line.")

  site['chr'] = fields[0]
  site['pos'] = fields[1]
  samples     = fields[9:]

  sample_counts = {}
  for i in range(len(samples)):
    
    variant_counts = {}
    counts = samples[i].split(':')[-1]
    counts = counts.split(',')

    for count in counts:
      if len(count) < 4:
        continue
      if count[0] != '-' and count[0] != '+':
        fail("Error in input VCF: variant data not strand-specific.")
      base = count[0:2]
      reads = count[3:]
      try:
        variant_counts[base.upper()] = int(reads)
      except ValueError, e:
        continue

    if i < len(sample_names):
      sample_counts[sample_names[i]] = variant_counts
    else:
      fail("Error in input VCF: more sample fields in data line than in header.")

  site['samples'] = sample_counts

  return site


def summarize_site(site, sample_names, bases, freq_thres, covg_thres):
  """Take the raw data from the VCF line and transform it into the summary data
  to be printed in the output format."""

  site_summary = []
  for sample_name in sample_names:

    sample = {}
    variants = site['samples'].get(sample_name)
    if not variants:
      fail("Error in input VCF: missing sample fields in data line.")

    sample['sample'] = sample_name
    sample['chr']    = site['chr']
    sample['pos']    = site['pos']

    coverage = 0
    top_base = ''
    top_base_reads = 0
    second_base = ''
    second_base_reads = 0
    for base in bases:
      reads = 0
      reads = variants.get('-'+base, 0)
      reads += variants.get('+'+base, 0)
      sample[base] = reads
      coverage += reads

      if reads > top_base_reads:
        second_base = top_base
        second_base_reads = top_base_reads
        top_base = base
        top_base_reads = reads
      elif reads > second_base_reads:
        second_base = base
        second_base_reads = reads

    if top_base_reads == second_base_reads:
      second_base = 'N'
      second_base_reads = 0

    sample['coverage'] = coverage
    sample['major']    = top_base or '.'
    sample['minor']    = second_base or '.'
    sample['freq']     = second_base_reads / coverage

    sample['alleles']  = count_alleles(sample, bases)

    # Coverage threshold
    if coverage < covg_thres:
      sample['print'] = False
    else:
      sample['print'] = True

    site_summary.append(sample)

  return site_summary


def print_site(site, columns):
  """Print the output lines for one site (one per sample)."""

  for sample in site:
    if sample['print']:
      fields = [str(sample.get(column)) for column in columns]
      print '\t'.join(fields)



def count_alleles(sample, bases):
  """Determine how many alleles to report, based on filtering rules."""
  alleles = 0
  # PLACEHOLDER (no filter)
  for base in bases:
    if sample[base]:
      alleles+=1
  return alleles


def fail(message):
  sys.stderr.write(message+'\n')
  sys.exit(1)

if __name__ == "__main__":
  main()