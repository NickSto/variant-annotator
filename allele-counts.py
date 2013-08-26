#!/usr/bin/python
# This parses the output of Dan's "Naive Variant Detector" (previously,
# "BAM Coverage"). It was forked from the code of "bam-coverage.py".
#
# New in this version:
#   Made header line customizable
#     - separate from internal column labels, which are used as dict keys
import os
import sys
from optparse import OptionParser

COLUMNS = ['sample', 'chr', 'pos', 'A', 'C', 'G', 'T', 'coverage', 'alleles', 'major', 'minor', 'freq'] #, 'bias']
COLUMN_LABELS = ['SAMPLE', 'CHR',  'POS', 'A', 'C', 'G', 'T', 'CVRG', 'ALLELES', 'MAJOR', 'MINOR', 'MINOR.FREQ.PERC.'] #, 'STRAND.BIAS']
CANONICAL_VARIANTS = ['A', 'C', 'G', 'T']
USAGE = """Usage: cat variants.vcf | %prog [options] > alleles.csv
       %prog [options] -i variants.vcf -o alleles.csv"""
OPT_DEFAULTS = {'infile':'-', 'outfile':'-', 'freq_thres':1.0, 'covg_thres':100,
  'print_header':False, 'stdin':False}
DESCRIPTION = """This will parse the VCF output of Dan's "Naive Variant Caller" (aka "BAM Coverage") Galaxy tool. For each position reported, it counts the number of reads of each base, determines the major allele, minor allele (second most frequent variant), and number of alleles above a threshold. So currently it only considers SNVs (ACGT), including in the coverage figure. By default it reads from stdin and prints to stdout."""
EPILOG = """Requirements:
The input VCF must report the variants for each strand.
The variants should be case-sensitive (e.g. all capital base letters).
Strand bias: Both strands must show the same bases passing the frequency threshold (but not necessarily in the same order). If the site fails this test, the number of alleles is reported as 0."""


def get_options(defaults, usage, description='', epilog=''):
  """Get options, print usage text."""

  parser = OptionParser(usage=usage, description=description, epilog=epilog)

  parser.add_option('-i', '--infile', dest='infile',
    default=defaults.get('infile'),
    help='Read input VCF data from this file instead of stdin.')
  parser.add_option('-o', '--outfile', dest='outfile',
    default=defaults.get('outfile'),
    help='Print output data to this file instead of stdout.')
  parser.add_option('-f', '--freq-thres', dest='freq_thres', type='float',
    default=defaults.get('freq_thres'),
    help='Frequency threshold for counting alleles, given in percentage: -f 1 = 1% frequency. Default is %default%.')
  parser.add_option('-c', '--covg-thres', dest='covg_thres', type='int',
    default=defaults.get('covg_thres'),
    help='Coverage threshold. Each site must be supported by at least this many reads on each strand. Otherwise the site will not be printed in the output. The default is %default reads per strand.')
  parser.add_option('-H', '--header', dest='print_header', action='store_const',
    const=not(defaults.get('print_header')), default=defaults.get('print_header'),
    help='Print header line. This is a #-commented line with the column labels. Off by default.')
  parser.add_option('-d', '--debug', dest='debug', action='store_true',
    default=False,
    help='Turn on debug mode. You must also specify a single site to process in a final argument using UCSC coordinate format.')

  (options, args) = parser.parse_args()

  # read in positional arguments
  arguments = {}
  if options.debug:
    if len(args) >= 1:
      arguments['print_loc'] = args[0]
      args.remove(args[0])

  return (options, arguments)


def main():

  (options, args) = get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)

  infile = options.infile
  outfile = options.outfile
  print_header = options.print_header
  freq_thres = options.freq_thres / 100.0
  covg_thres = options.covg_thres
  debug = options.debug

  if debug:
    print_loc = args.get('print_loc')
    if print_loc:
      if ':' in print_loc:
        (print_chr, print_pos) = print_loc.split(':')
      else:
        print_pos = print_loc
    else:
      sys.stderr.write("Warning: No site coordinate found in arguments. "
        +"Turning off debug mode.\n")
      debug = False

  # set infile_handle to either stdin or the input file
  if infile == OPT_DEFAULTS.get('infile'):
    infile_handle = sys.stdin
    sys.stderr.write("Reading from standard input..\n")
  else:
    if os.path.exists(infile):
      infile_handle = open(infile, 'r')
    else:
      fail('Error: Input VCF file '+infile+' not found.')

  # set outfile_handle to either stdout or the output file
  if outfile == OPT_DEFAULTS.get('outfile'):
    outfile_handle = sys.stdout
  else:
    try:
      outfile_handle = open(outfile, 'w')
    except IOError, e:
      fail('Error: The given output filename '+outfile+' could not be opened.')

  if len(COLUMNS) != len(COLUMN_LABELS):
    fail('Error: Internal column names do not match column labels.')
  if print_header:
    outfile_handle.write('#'+'\t'.join(COLUMN_LABELS)+"\n")

  # main loop: process and print one line at a time
  sample_names = []
  for line in infile_handle:
    line = line.rstrip('\r\n')

    # header lines
    if line[0] == '#':
      if line[0:6].upper() == '#CHROM':
        sample_names = line.split('\t')[9:]
      continue

    if not sample_names:
      fail("Error in input VCF: Data line encountered before header line. "
        +"Failed on line:\n"+line)

    site_data = read_site(line, sample_names, CANONICAL_VARIANTS)

    if debug:
      if site_data['pos'] != print_pos:
        continue
      try:
        if site_data['chr'] != print_chr:
          continue
      except NameError, e:
        pass  # No chr specified. Just go ahead and print the line.

    site_summary = summarize_site(site_data, sample_names, CANONICAL_VARIANTS,
      freq_thres, covg_thres, debug=debug)

    if debug and site_summary[0]['print']:
      print line.split('\t')[9].split(':')[-1]

    print_site(outfile_handle, site_summary, COLUMNS)

  # close any open filehandles
  if infile_handle is not sys.stdin:
    infile_handle.close()
  if outfile_handle is not sys.stdout:
    outfile_handle.close()

  # keeps Galaxy from giving an error if there were messages on stderr
  sys.exit(0)



def read_site(line, sample_names, canonical):
  """Read in a line, parse the variants into a data structure, and return it.
  The line should be actual site data, not a header line, so check beforehand.
  Notes:
  - The line is assumed to have been chomped."""
  
  site = {}
  fields = line.split('\t')

  if len(fields) < 9:
    fail("Error in input VCF: wrong number of fields in data line. "
          +"Failed on line:\n"+line)

  site['chr'] = fields[0]
  site['pos'] = fields[1]
  samples     = fields[9:]

  if len(samples) < len(sample_names):
    fail("Error in input VCF: missing sample fields in data line. "
          +"Failed on line:\n"+line)
  elif len(samples) > len(sample_names):
    fail("Error in input VCF: more sample fields in data line than in header. "
          +"Failed on line:\n"+line)

  sample_counts = {}
  for i in range(len(samples)):
    
    variant_counts = {}
    counts = samples[i].split(':')[-1]
    counts = counts.split(',')

    for count in counts:
      if not count:
        continue
      fields = count.split('=')
      if len(fields) != 2:
        fail("Error in input VCF: Incorrect variant data format (must contain "
          +"a single '='). Failed on line:\n"+line)
      (variant, reads) = fields
      if variant[1:] not in canonical:
        continue
      if variant[0] != '-' and variant[0] != '+':
        fail("Error in input VCF: variant data not strand-specific. "
          +"Failed on line:\n"+line)
      try:
        variant_counts[variant] = int(float(reads))
      except ValueError, e:
        fail("Error in input VCF: Variant count not a valid number. Failed on variant count string '"+reads+"'\nIn the following line:\n"+line)

    sample_counts[sample_names[i]] = variant_counts

  site['samples'] = sample_counts

  return site


def summarize_site(site, sample_names, canonical, freq_thres, covg_thres,
  debug=False):
  """Take the raw data from the VCF line and transform it into the summary data
  to be printed in the output format."""

  site_summary = []
  for sample_name in sample_names:

    sample = {'print':False}
    variants = site['samples'].get(sample_name)
    if not variants:
      site_summary.append(sample)
      continue

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
    if coverage <= 0 or covg_plus < covg_thres or covg_minus < covg_thres:
      site_summary.append(sample)
      continue
    else:
      sample['print'] = True

    # get an ordered list of read counts for all variants (either strand)
    ranked_bases = get_read_counts(variants, 0, strands='+-', debug=debug)

    # record read counts into dict for this sample
    for base in ranked_bases:
      sample[base[0]] = base[1]
    # fill in any zeros
    for variant in canonical:
      if not sample.has_key(variant):
        sample[variant] = 0

    sample['alleles']  = count_alleles(variants, freq_thres, debug=debug)

    # set minor allele to N if there's a tie for 2nd
    if len(ranked_bases) >= 3 and ranked_bases[1][1] == ranked_bases[2][1]:
      ranked_bases[1] = ('N', 0)
      sample['alleles'] = 1 if sample['alleles'] else 0

    if debug: print ranked_bases

    sample['coverage'] = coverage
    try:
      sample['major']  = ranked_bases[0][0]
    except IndexError, e:
      sample['major']  = '.'
    try:
      sample['minor']  = ranked_bases[1][0]
      sample['freq']   = round(ranked_bases[1][1]/float(coverage), 5)
    except IndexError, e:
      sample['minor']  = '.'
      sample['freq']   = 0.0

    site_summary.append(sample)

  return site_summary


def print_site(filehandle, site, columns):
  """Print the output lines for one site (one per sample).
  filehandle must be open."""
  for sample in site:
    if sample['print']:
      fields = [str(sample.get(column)) for column in columns]
      filehandle.write('\t'.join(fields)+"\n")


def get_read_counts(variant_counts, freq_thres, strands='+-', debug=False):
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

  # Get list of all variants from variant_counts list
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
  # if debug: print "strands: "+strands+', covg: '+str(coverage)

  if coverage < 1:
    return []

  # sort the list of alleles by read count
  ranked_bases.sort(reverse=True, key=lambda base: base[1])

  if debug:
    print strands+' coverage: '+str(coverage)+', freq_thres: '+str(freq_thres)
    for base in ranked_bases:
      print (base[0]+': '+str(base[1])+'/'+str(float(coverage))+' = '+
        str(base[1]/float(coverage)))

  # remove bases below the frequency threshold
  ranked_bases = [base for base in ranked_bases
    if base[1]/float(coverage) >= freq_thres]

  return ranked_bases


def count_alleles(variant_counts, freq_thres, debug=False):
  """Determine how many alleles to report, based on filtering rules.
  The current rule determines which bases pass the frequency threshold on each
  strand individually, then compares the two sets of bases. If they are the same
  (regardless of order), the allele count is the number of bases. Otherwise it
  is zero."""
  allele_count = 0

  alleles_plus  = get_read_counts(variant_counts, freq_thres, debug=debug,
    strands='+')
  alleles_minus = get_read_counts(variant_counts, freq_thres, debug=debug,
    strands='-')

  if debug:
    print '+ '+str(alleles_plus)
    print '- '+str(alleles_minus)

  # check if each strand reports the same set of alleles
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