import sys
import os
import re

def dist(seqa, seqb):
	""" returns hamming distance of the sequences seqa and seqb """
	return sum(0 if a == b else 1 for a, b in zip(seqa, seqb))

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

# https://github.com/SergioFierens/ai4r/blob/master/lib/ai4r/clusterers/single_linkage.rb
#def create_distance_matrix(data_set)
#	@distance_matrix = Array.new(data_set.data_items.length-1) {|index| Array.new(index+1)}
#	data_set.data_items.each_with_index do |a, i|
#		i.times do |j|
#			b = data_set.data_items[j]
#			@distance_matrix[i-1][j] = @distance_function.call(a, b)
#		end
#	end
#end

def create_distance_matrix(alignment):
	return [[dist(seqa, seqb) for seqb in alignment] for seqa in alignment]

def upgma(alignment, d):
	for seq in alignment:
		skip
	return

alnm = read_alignment("data/test.aln")
for row in create_distance_matrix(alnm):
	print row
