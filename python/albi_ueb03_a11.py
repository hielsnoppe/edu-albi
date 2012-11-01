import sys
import os
import re
import itertools
import heapq

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

def create_distance_matrix(alignment):
	return [[dist(seqa, seqb) for seqb in alignment] for seqa in alignment]

def upgma(alignment, d):
	indices = [i for i in range(len(d))]
	clusters = dict([(i, i) for i in range(len(d))])
	next = len(indices)
	while len(indices) > 1:
		candidates = [(d[a][b], a, b) for a, b in itertools.combinations(indices, 2)]
		heapq.heapify(candidates)
		dist, a, b = heapq.heappop(candidates)
		for row in d:
			row.append((row[a] + row[b]) / 2)		# compute entry for new column
		d.append([row[next] for row in d] + [0])	# copy last column to new row
		indices.remove(a)
		indices.remove(b)
		indices.append(next)
		clusters[next] = (clusters[a], clusters[b])
		next += 1
	return clusters[next-1]

alignment = read_alignment("../data/test.aln")
#alignment = {"1":"ACG--TA", "2":"ACAGGTA", "3":"ACGA-TA", "4":"CC-GGTA"}

d = create_distance_matrix(alignment)
tree = upgma(alignment, d)
print (tree)
