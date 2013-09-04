#!/usr/bin/python
import os
import sys
from optparse import OptionParser

OPT_DEFAULTS = {'freq_thres':0, 'covg_thres':0}
BASES = ['A', 'C', 'G', 'T']
USAGE = ("Usage: %prog (options) 'variant_str' ('variant_str2' (etc))\n"
    +"       cat variants.txt > %prog (options)")

def main():

  parser = OptionParser(usage=USAGE)
  parser.add_option('-f', '--freq-thres', dest='freq_thres', type='float',
    default=OPT_DEFAULTS.get('freq_thres'),
    help=('Frequency threshold for counting alleles, given in percentage: -f 1 '
      +'= 1% frequency. Default is %default%.'))
  parser.add_option('-c', '--covg-thres', dest='covg_thres', type='int',
    default=OPT_DEFAULTS.get('covg_thres'),
    help=('Coverage threshold. Each site must be supported by at least this '
      +'many reads on each strand. Otherwise the site will not be printed in '
      +'the output. The default is %default reads per strand.'))
  (options, args) = parser.parse_args()
  freq_thres = options.freq_thres
  covg_thres = options.covg_thres

  if len(sys.argv) > 1 and '-h' in sys.argv[1][0:3]:
    script_name = os.path.basename(sys.argv[0])
    print """USAGE:
  $ """+script_name+""" [sample column text]
  $ """+script_name+""" '+A=10,+G=2,-A=20,-G=41,'
  $ """+script_name+""" '0/1:10,1:0.25,0.025:-A=29,-AG=1,-G=10,'
  $ """+script_name+""" '+A=10,+G=2,-A=20,-G=41,' '0/1:10,1:0.25,0.025:-A=29,-AG=1,-G=10,'
Or invoke with no arguments to use interactively. It will read from stdin, so
just paste one sample per line."""
    sys.exit(0)
  
  if len(args) > 0:
    stdin = False
    samples = args
  else:
    stdin = True
    samples = sys.stdin
    print "Reading from standard input.."

  for sample in samples:
    print ''
    sample = sample.split(':')[-1]
    if not sample:
      continue
    print sample
    counts_dict = parse_counts(sample)
    compute_stats(counts_dict, freq_thres, covg_thres)
  

def parse_counts(sample_str):
  counts_dict = {}
  counts = sample_str.split(',')
  for count in counts:
    if '=' not in count:
      continue
    (var, reads) = count.split('=')
    if var[1:] in BASES:
      counts_dict[var] = int(reads)
  return counts_dict


def compute_stats(counts_dict, freq_thres, covg_thres):
  # totals for A, C, G, T
  counts_unstranded = {}
  for base in BASES:
    counts_unstranded[base] = 0
    for strand in '+-':
      counts_unstranded[base] += counts_dict.get(strand+base, 0)
  print '+- '+str(counts_unstranded)

  # get counts for each strand
  plus_counts  = get_stranded_counts(counts_dict, '+')
  minus_counts = get_stranded_counts(counts_dict, '-')
  counts_lists = {'+':plus_counts, '-':minus_counts}
  print ' + '+str(plus_counts)
  print ' - '+str(minus_counts)

  # stranded coverage threshold
  coverages = {}
  failing_strands = {}
  for strand in '+-':
    coverages[strand] = 0
    for count in counts_lists[strand].values():
      coverages[strand] += count
    if coverages[strand] < covg_thres:
      failing_strands[strand] = coverages[strand]
    sys.stdout.write(strand+'coverage: '+str(coverages[strand])+"\t")
  coverages['+-'] = coverages['+'] + coverages['-']
  sys.stdout.write("+-coverage: "+str(coverages['+-'])+"\n")

  if failing_strands:
    for strand in failing_strands:
      print ('coverage on '+strand+' strand too low ('
        +str(failing_strands[strand])+' < '+str(covg_thres)+")")
    return

  # apply frequency threshold
  for strand in counts_lists:
    strand_counts = counts_lists[strand]
    for variant in strand_counts.keys():
      # print (variant+" freq: "+str(strand_counts[variant])+"/"
      #   +str(coverages[strand])+" = "
      #   +str(strand_counts[variant]/float(coverages[strand])))
      if strand_counts[variant]/float(coverages[strand]) < freq_thres:
        strand_counts.pop(variant)
  plus_variants  = sorted(plus_counts.keys())
  minus_variants = sorted(minus_counts.keys())
  if plus_variants == minus_variants:
    strand_bias = False
    print "no strand bias: +"+str(plus_variants)+" == -"+str(minus_variants)
    sys.stdout.write("alleles: "+str(len(plus_variants))+"\t")
  else:
    strand_bias = True
    print "   strand bias: +"+str(plus_variants)+" != -"+str(minus_variants)
    sys.stdout.write("alleles: 0\t")

  variants_sorted = sort_variants(counts_unstranded)
  if len(variants_sorted) >= 1:
    sys.stdout.write("major: "+variants_sorted[0]+"\t")
  minor = '.'
  if len(variants_sorted) == 2:
    minor = variants_sorted[1]
  elif len(variants_sorted) > 2:
    if (counts_unstranded.get(variants_sorted[1]) ==
        counts_unstranded.get(variants_sorted[2])):
      minor = 'N'
    else:
      minor = variants_sorted[1]

  sys.stdout.write("minor: "+minor+"\tfreq: ")
  print round(float(counts_unstranded.get(minor, 0))/coverages['+-'], 5)


def get_stranded_counts(unstranded_counts, strand):
  stranded_counts = {}
  for variant in unstranded_counts:
    if variant[0] == strand:
      stranded_counts[variant[1:]] = unstranded_counts[variant]
  return stranded_counts

def sort_variants(variant_counts):
  """Sort the list of variants based on their counts. Returns a list of just
  the variants, no counts."""
  variants = variant_counts.keys()
  var_del = []
  for variant in variants:
    if variant_counts.get(variant) == 0:
      var_del.append(variant)
  for variant in var_del:
    variants.remove(variant)
  variants.sort(reverse=True, key=lambda variant: variant_counts.get(variant,0))
  return variants

if __name__ == "__main__":
  main()
