#!/usr/bin/python
# This parses the output of Dan's "Naive Variant Detector" (previously,
# "BAM Coverage"). It was forked from the code of "bam-coverage.py".
#
# TODO:
# - use actual flagged options
# - option to include a header in the output
# - should it technically handle data lines that start with a '#'?
import os
import sys

debug = True
# debug = False

COLUMNS = ['sample', 'chr', 'pos', 'A', 'C', 'G', 'T', 'coverage', 'alleles',
  'major', 'minor', 'freq'] #, 'bias']
BASES = ['A', 'C', 'G', 'T', 'N']
PRINT_HEADER_DEFAULT = False
FREQ_THRES_DEFAULT = 1/100
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
input: The input VCF must report the variants per strand. Note: the variants
  are case-sensitive.
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
      freq_thres = int(sys.argv[2])/100.0
    except ValueError, e:
      pass
  if args > 3:
    try:
      covg_thres = int(sys.argv[3])
    except ValueError, e:
      pass

  global debug
  if debug:
    if args > 4:
      if ':' in sys.argv[4]:
        (print_chr, print_pos) = sys.argv[4].split(':')
      else:
        print_pos = sys.argv[4]
    else:
      # pass
      debug = False

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

      if debug:
        if site_data['pos'] != print_pos:
          continue
        try:
          if site_data['chr'] != print_chr:
            continue
        except NameError, e:
          "No chr specified. Just go ahead and print the line"


      site_summary = summarize_site(site_data, sample_names, BASES, freq_thres,
        covg_thres)

      if debug and site_summary[0]['print']: print line.split('\t')[9].split(':')[-1]

      print_site(site_summary, COLUMNS)

      # return



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
      if not count:
        continue
      (variant, reads) = count.split('=')
      if variant[0] != '-' and variant[0] != '+':
        fail("Error in input VCF: variant data not strand-specific.")
      try:
        variant_counts[variant] = int(reads)
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

    coverage = sum(variants.values())

    # get stranded coverage
    covg_plus = 0
    covg_minus = 0
    for variant in variants:
      if variant[0] == '+':
        covg_plus += variants[variant]
      elif variant[0] == '-':
        covg_minus += variants[variant]
    # stranded coverage threshold
    if covg_plus < covg_thres or covg_minus < covg_thres:
      sample['print'] = False
      site_summary.append(sample)
      continue
    else:
      sample['print'] = True

    # get an ordered list of read counts for all variants (either strand)
    ranked_bases = get_read_counts(variants, 0, strands='+-')

    # record read counts into dict for this sample
    for base in ranked_bases:
      sample[base[0]] = base[1]
    # fill in any zeros
    for base in bases:
      if not sample.has_key(base):
        sample[base] = 0

    # set minor allele to N if there's a tie for 2nd
    if len(ranked_bases) >= 3 and ranked_bases[1][1] == ranked_bases[2][1]:
      ranked_bases[1] = ('N', 0)

    if debug: print ranked_bases

    sample['coverage'] = coverage
    sample['major']    = ranked_bases[0][0]
    try:
      sample['minor']  = ranked_bases[1][0]
      sample['freq']   = ranked_bases[1][1] / float(coverage)
    except IndexError, e:
      sample['minor']  = '.'
      sample['freq']   = 0.0

    sample['alleles']  = count_alleles(variants, freq_thres, bases)

    site_summary.append(sample)

  return site_summary


def print_site(site, columns):
  """Print the output lines for one site (one per sample)."""
  for sample in site:
    if sample['print']:
      fields = [str(sample.get(column)) for column in columns]
      print '\t'.join(fields)


def get_read_counts(variant_counts, freq_thres, strands='+-', variants=None):
  """Count the number of reads for each base, and create a ranked list of
  alleles passing the frequency threshold.
      Arguments:
  variant_counts: Dict of the stranded variants (keys) and their read counts (values).
  freq_thres: The frequency threshold each allele needs to pass to be included.
  strands: Which strand(s) to count. Can be '+', '-', or '+-' for both (default).
  variants: A list of the variants of interest. Other types of variants will not
    be included in the returned list. If no list is given, all variants found in
    the variant_counts will be used.
      Return value:
  ranked_bases: A list of the alleles and their read counts. The elements are
    tuples (base, read count). The alleles are listed in descending order of
    frequency, and only those passing the threshold are included."""

  # If variant list is not given, create one from the variant_counts list
  if not variants:
   variants = [variant[1:] for variant in variant_counts]
  # deduplicate via a dict
  variant_dict = dict((variant, 1) for variant in variants)
  variants = variant_dict.keys()

  ranked_bases = []
  for variant in variants:
    reads = 0
    for strand in strands:
      reads += variant_counts.get(strand+variant, 0)
    ranked_bases.append((variant, reads))

  # get coverage for the specified strands
  coverage = 0
  for variant in variant_counts:
    if variant[0] in strands:
      coverage += variant_counts.get(variant, 0)

  # sort the list of alleles by read count
  ranked_bases.sort(reverse=True, key=lambda base: base[1])

  if debug:
    print 'coverage: '+str(coverage)+', freq_thres: '+str(freq_thres)
    for base in ranked_bases:
      print (base[0]+': '+str(base[1])+'/'+str(float(coverage))+' = '+
        str(base[1]/float(coverage)))

  # remove bases below the frequency threshold
  ranked_bases = [base for base in ranked_bases
    if base[1]/float(coverage) >= freq_thres]

  return ranked_bases


def count_alleles(variant_counts, freq_thres, bases):
  """Determine how many alleles to report, based on filtering rules.
  The current rule determines which bases pass the frequency threshold on each
  strand individually, then compares the two sets of bases. If they are the same
  (regardless of order), the allele count is the number of bases. Otherwise it
  is zero."""
  allele_count = 0

  alleles_plus  = get_read_counts(variant_counts, freq_thres, variants=bases,
    strands='+')
  alleles_minus = get_read_counts(variant_counts, freq_thres, variants=bases,
    strands='-')
  alleles_plus_sorted  = sorted([base[0] for base in alleles_plus if base[1]])
  alleles_minus_sorted = sorted([base[0] for base in alleles_minus if base[1]])

  if alleles_plus_sorted == alleles_minus_sorted:
    allele_count = len(alleles_plus)

  return allele_count


def fail(message):
  sys.stderr.write(message+'\n')
  sys.exit(1)

if __name__ == "__main__":
  main()