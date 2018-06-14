#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Written by Wu, Yang <wuyang@vt.edu>

'''
Usage:
python stringtie_FPKM.py <ballgown_dir>|<sample_list>

For example:
  python stringtie_FPKM.py /foo/stringtie_output/
Or (it is better if the groups are ordered, i.e. Group1, Group2, Group3...):
  python stringtie_FPKM.py sample_list.tsv

sample_list.tsv content (tabluar format):
/foo/Sample1/\tGroup1
/foo/Sample2/\tGroup1
/foo/Sample3/\tGroup2
/foo/Sample4/\tGroup2
'''

import sys, subprocess, os, re
from copy import deepcopy

try:
	arg = sys.argv[1]
except IndexError:
	print __doc__
	sys.exit(1)
if os.path.isfile(arg):
	with open(arg, 'r') as inf:
		filelist = [tuple(line.strip().split('\t')) for line in inf]
elif os.path.isdir(arg):
	filelist = [(f, os.path.split(f)[0]. split('/')[-1]) \
	            for f in subprocess.check_output('find %s -name "*.gtf" -type f | sort -u' % arg,
	                                             shell=True).strip().split('\n')]
	print filelist

FPKM_dict = {}
groups = []
i = 0
max_len = 0
total = len(filelist)
gene_id_p = r'gene_id "(.+?)";'
trx_id_p = r'transcript_id "(.+?)";'
ref_gene_name_p = r'ref_gene_name "(.+?)";'
FPKM_p = r'FPKM "(.+?)";'
#TPM_p = r'TPM "(.+?)";'

for item in filelist:
	i += 1
	filename = item[0]
	if os.path.isdir(filename):
		#filename += '/t_data.ctab'
		if filename[-1] != '/':
			filename += '/'
		filename += subprocess.check_output('ls %s|grep "\.gtf$"' % filename, shell=True).strip()
	group = item[-1]
	info = "(%s/%s) %s" % (i, total, group)
	if len(info) < max_len:
		blank = max_len - len(info)
	else:
		max_len = len(info)
		blank = 0
	print(info + ' '*blank + '\r'),
	#sys.stdout.write('(%s/%s) %s\r' % (i, total, group))
	sys.stdout.flush()
	groups.append(group)
	with open(filename, 'r') as inf:
		for line in inf:
			line = line.strip().split('\t')
			try:
				if line[2] == "transcript":
					info = line[-1]
					# May vary:
					#gene_id = re.findall(ref_gene_name_p, info)
					gene_id = re.findall(gene_id_p, info)
					FPKM = re.findall(FPKM_p, info)
					if len(gene_id) == 1 and len(FPKM) == 1:
						gene_id = gene_id[0]
						FPKM = FPKM[0]
						FPKM_stats = deepcopy(FPKM_dict.get(gene_id, {}))
						FPKM_stats[group] = FPKM
						FPKM_dict[gene_id] = deepcopy(FPKM_stats)
			except IndexError:
				pass

with open("FPKM_stats.tsv", 'w') as outf:
	print "\nSummarizing..."
	n = 0
	for item in sorted(FPKM_dict.iteritems()):
		n += 1
		if n == 1:
			outf.write("Gene_id\t" + '\t'.join(groups) + "\n")
		FPKMs = [item[-1].get(g, "NA") for g in groups]
		outf.write(item[0] + '\t' + "\t".join(FPKMs) + '\n')
print "All set."