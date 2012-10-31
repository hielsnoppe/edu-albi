import sys
import os
import re

def read_alignment(filename):
	with open(filename, 'r') as seqfile:
		seq = {}
		regex = r"^>(\S+)\s"
		name = ""
		for line in seqfile:
			c = re.match(regex, line)
			if (c):
				name = c.group(1)
				seq[name] = ""
			else:
				seq[name] += line.strip()
	return seq

def alignment_to_distancematrix(alignment):
	return

def distancematrix_to_phylogenetictree():
	return

seq = read_alignment("data/test.aln")
print seq
