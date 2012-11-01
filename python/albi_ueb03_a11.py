import sys
import os
import re
import itertools
import heapq
import argparse

def dist(seqa, seqb):
	""" returns normalized hamming distance of the sequences seqa and seqb """
	return sum(0.0 if a == b else 1.0 for a, b in zip(seqa, seqb)) / len(seqa)

def jcest(dist, len):
	""" returns estimated evolutionary distance for proportion dist/len of differences """
	p = float(dist) / float(len)
	critical = 4.0 * p / 3.0
	if critical < 1.0:
		return (-3.0) * math.log(1.0 - critical) / 4.0
	else:
		return 0.0

def read_alignment(filename):
	with open(filename, 'r') as f:
		seq = {}
		regex = r"^>(\S+)\s"
		name = ""
		for line in f:
			c = re.match(regex, line)
			if (c):
				name = c.group(1)			# get name for next sequence
				seq[name] = ""				# create empty entry for next sequence
			else:
				seq[name] += line.strip()	# add line to current sequence
	return seq

def create_distance_matrix(alignment):
	return [[dist(seqa, seqb) for seqb in alignment] for seqa in alignment]

def upgma(alignment, d):
	clusters = dict([(i, i) for i in range(len(d))])
	next = len(clusters)
	while len(clusters) > 1:
		candidates = [(d[a][b], a, b) for a, b in itertools.combinations(clusters.keys(), 2)]
		heapq.heapify(candidates)
		dist, a, b = heapq.heappop(candidates)
		for row in d:
			row.append((row[a] + row[b]) / 2)		# compute entry for new column
		d.append([row[next] for row in d] + [0])	# copy last column to new row
		clusters[next] = (clusters[a], clusters[b])
		del(clusters[a])
		del(clusters[b])
		next += 1
	return clusters[next-1]



parser = argparse.ArgumentParser(description="Uebung 3 - Aufgabe 11")

parser.add_argument('-i', '--input',
type=str,
required=True,
help='Eingabedatei')

parser.add_argument('-o', '--output',
type=str,
required=False,
help='Ausgabedatei')

args = parser.parse_args()
infile = args.input
outfile = args.output

a = read_alignment(infile)
d = create_distance_matrix(a)
tree = upgma(a, d)

if outfile:
	with open(outfile, 'w') as f:
		f.write(str(tree))
else:
	print(tree)
