#!/usr/bin/python
# This parses the output of Dan's "Naive Variant Detector" (previously,
# "BAM Coverage"). It was forked from the code of "bam-coverage.py".
# Or, to put it briefly,
# cat variants.vcf | grep -v '^#' | cut -f 10 | cut -d ':' -f 4 | tr ',=' '\t:'
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
USAGE = """Usage: %prog [options] -i variants.vcf -o alleles.csv
       cat variants.vcf | %prog [options] > alleles.csv"""
OPT_DEFAULTS = {'infile':'-', 'outfile':'-', 'freq_thres':1.0, 'covg_thres':100,
  'print_header':False, 'stdin':False, 'stranded':False}
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
    help=('Frequency threshold for counting alleles, given in percentage: -f 1 '
      +'= 1% frequency. Default is %default%.'))
  parser.add_option('-c', '--covg-thres', dest='covg_thres', type='int',
    default=defaults.get('covg_thres'),
    help=('Coverage threshold. Each site must be supported by at least this '
      +'many reads on each strand. Otherwise the site will not be printed in '
      +'the output. The default is %default reads per strand.'))
  parser.add_option('-H', '--header', dest='print_header', action='store_const',
    const=not(defaults.get('print_header')), default=defaults.get('print_header'),
    help=('Print header line. This is a #-commented line with the column '
      +'labels. Off by default.'))
  parser.add_option('-s', '--stranded', dest='stranded', action='store_const',
    const=not(defaults.get('stranded')), default=defaults.get('stranded'),
    help='Report variant counts by strand, in separate columns. Off by default.')
  parser.add_option('-d', '--debug', dest='debug', action='store_true', default=False,
    help=('Turn on debug mode. You must also specify a single site to process '
      +'in a final argument using UCSC coordinate format. Also, you can add a '
      +'sample ID after another ":" to restrict it further.'))

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
  stranded = options.stranded
  debug = options.debug

  print_sample = ''
  if debug:
    print_loc = args.get('print_loc')
    if print_loc:
      coords = print_loc.split(':')
      print_chr = coords[0]
      print_pos = ''
      if len(coords) > 1: print_pos = coords[1]
      if len(coords) > 2: print_sample = coords[2].lower()
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

  # Take care of column names, print header
  if len(COLUMNS) != len(COLUMN_LABELS):
    fail('Error: Internal column names list do not match column labels list.')
  if stranded:
    COLUMNS[3:7]       = ['+A', '+C', '+G', '+T', '-A', '-C', '-G', '-T']
    COLUMN_LABELS[3:7] = ['+A', '+C', '+G', '+T', '-A', '-C', '-G', '-T']
  if print_header:
    outfile_handle.write('#'+'\t'.join(COLUMN_LABELS)+"\n")

  # main loop
  # each iteration processes one VCF line and prints one or more output lines
  # one VCF line    = one site, one or more samples
  # one output line = one site, one sample
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
      if print_pos != site_data['pos']:
        continue
      if print_chr != site_data['chr'] and print_chr != '':
        continue
      if print_sample != '':
        for sample in site_data['samples'].keys():
          if sample.lower() != print_sample:
            site_data['samples'].pop(sample, None)


    site_summary = summarize_site(site_data, sample_names, CANONICAL_VARIANTS,
      freq_thres, covg_thres, stranded, debug=debug)

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
  Only the variants in 'canonical' will be read; all others are ignored.
  Note: the line is assumed to have been chomped.
  The returned data is stored in a dict, with the following structure:
  {
    'chr': 'chr1',
    'pos': '2617',
    'samples': {
      'THYROID': {
        '+A': 32,
        '-A': 45,
        '-G': 1,
      },
      'BLOOD': {
        '+A': 2,
        '-C': 1,
        '+G': 37,
        '-G': 42,
      },
    },
  }
  """
  
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
  stranded=False, debug=False):
  """Take the raw data from the VCF line and transform it into the summary data
  to be printed in the output format."""

  site_summary = []
  for sample_name in sample_names:

    sample = {'print':False}
    variants = site['samples'].get(sample_name)

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
      site_summary.append(sample)
      continue
    else:
      sample['print'] = True

    # get an ordered list of read counts for all variants (both strands)
    bases = get_read_counts(variants, '+-')
    ranked_bases = process_read_counts(bases, sort=True, debug=debug)

    # prepare stranded or unstranded lists of base counts
    base_count_lists = []
    if stranded:
      strands = ['+', '-']
      base_count_lists.append(get_read_counts(variants, '+'))
      base_count_lists.append(get_read_counts(variants, '-'))
    else:
      strands = ['']
      base_count_lists.append(ranked_bases)

    # record read counts into output dict
    # If stranded, this will loop twice, once for each strand, and prepend '+'
    # or '-' to the base name. If not stranded, it will loop once, and prepend
    # nothing ('').
    for (strand, base_count_list) in zip(strands, base_count_lists):
      for base_count in base_count_list:
        sample[strand+base_count[0]] = base_count[1]
      # fill in any zeros
      for base in canonical:
        if not sample.has_key(strand+base):
          sample[strand+base] = 0

    sample['alleles'] = count_alleles(variants, freq_thres, debug=debug)

    if debug: print "ranked +-: "+str(ranked_bases)

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

    # set minor allele to N if there's a tie for 2nd
    if len(ranked_bases) >= 3 and ranked_bases[1][1] == ranked_bases[2][1]:
      sample['minor'] = 'N'
      sample['freq'] = 0.0

    site_summary.append(sample)

  return site_summary


def get_read_counts(stranded_counts, strands='+-'):
  """Do a simple sum of the read counts per variant, on the specified strands.
      Arguments:
  stranded_counts: Dict of the stranded variants (keys) and their read counts
    (values).
  strands: Which strand(s) to count. Can be '+', '-', or '+-' for both (default).
      Return value:
  summed_counts: A list of the alleles and their read counts. The elements are
    tuples (variant, read count)."""

  variants = stranded_counts.keys()

  summed_counts = {}
  for variant in variants:
    strand = variant[0]
    base = variant[1:]
    if strand in strands:
      summed_counts[base] = stranded_counts[variant] + summed_counts.get(base, 0)

  return summed_counts.items()


def process_read_counts(variant_counts, freq_thres=0, sort=False, debug=False):
  """Process a list of read counts by frequency filtering and/or sorting.
      Arguments:
  variant_counts: List of the non-stranded variants and their read counts. The
    elements are tuples (variant, read count).
  freq_thres: The frequency threshold each allele needs to pass to be included.
  sort: Whether to sort the list in descending order of read counts.
      Return value:
  variant_counts: A list of the processed alleles and their read counts. The
    elements are tuples (variant, read count)."""

  # get coverage for the specified strands
  coverage = 0
  for variant in variant_counts:
    coverage += variant[1]

  if coverage <= 0:
    return []

  # sort the list of alleles by read count
  if sort:
    variant_counts.sort(reverse=True, key=lambda variant: variant[1])

  if debug:
    print 'coverage: '+str(coverage)+', freq_thres: '+str(freq_thres)
    for variant in variant_counts:
      print (variant[0]+': '+str(variant[1])+'/'+str(float(coverage))+' = '+
        str(variant[1]/float(coverage)))

  # remove bases below the frequency threshold
  if freq_thres > 0:
    variant_counts = [variant for variant in variant_counts
      if variant[1]/float(coverage) >= freq_thres]

  return variant_counts


def count_alleles(variant_counts, freq_thres, debug=False):
  """Determine how many alleles to report, based on filtering rules.
  The current rule determines which bases pass the frequency threshold on each
  strand individually, then compares the two sets of bases. If they are the same
  (regardless of order), the allele count is the number of bases. Otherwise it
  is zero."""
  allele_count = 0

  alleles_plus  = get_read_counts(variant_counts, '+')
  alleles_plus  = process_read_counts(alleles_plus, freq_thres=freq_thres,
    sort=False, debug=debug)
  alleles_minus = get_read_counts(variant_counts, '-')
  alleles_minus = process_read_counts(alleles_minus, freq_thres=freq_thres,
    sort=False, debug=debug)

  if debug:
    print '+ '+str(alleles_plus)
    print '- '+str(alleles_minus)

  # Check if each strand reports the same set of alleles.
  # Sorting by base is to compare lists without regard to order (as sets).
  alleles_plus_sorted  = sorted([base[0] for base in alleles_plus if base[1]])
  alleles_minus_sorted = sorted([base[0] for base in alleles_minus if base[1]])
  if alleles_plus_sorted == alleles_minus_sorted:
    allele_count = len(alleles_plus)

  return allele_count


def print_site(filehandle, site, columns):
  """Print the output lines for one site (one per sample).
  filehandle must be open."""
  for sample in site:
    if sample['print']:
      fields = [str(sample.get(column)) for column in columns]
      filehandle.write('\t'.join(fields)+"\n")


def fail(message):
  sys.stderr.write(message+'\n')
  sys.exit(1)


if __name__ == "__main__":
  main()