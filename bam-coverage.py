#!/usr/bin/python
# For parsing the output of Dan's Bam Coverage tool
import sys

HEADER = "#location\tdepth\tA's\tT's\tG's\tC's\tN's"

if len(sys.argv) > 1:
  filename = sys.argv[1]
else:
  sys.stderr.write("Error: must give an input VCF file as first argument.\n")
  sys.exit(1)

names = []
snps = []
with open(filename, 'r') as lines:
  for line in lines:
    snp = {}
    variants = []
    line = line.rstrip('\r\n')
    if line[0:2] == '##':
      continue

    fields = line.split('\t')
    if len(fields) < 9:
      sys.stderr.write("Error: wrong number of fields. Skipping line.\n")
      continue
    snp['chr'] = fields[0]
    snp['loc'] = fields[1]
    snp['id']  = fields[2]
    snp['ref'] = fields[3]
    samples    = fields[9:]
    if line[0:6].lower() == '#chrom':
      names = samples
      print "#" + "\t".join(names)
      continue

    for sample in samples:
      variant = {}
      counts = sample.split(':')[1]
      counts = counts.split(',')
      for count in counts:
        if count == '': continue
        (base, reads) = count.split('-')
        variant[base.lower()] = int(reads)

      variants.append(variant)

    snp['variants'] = variants
    snps.append(snp)

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
  sys.stdout.write(snp['chr']+':'+snp['loc']+'\t')
  print line
