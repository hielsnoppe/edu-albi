import sys
import os
import argparse
import math
import random

def rseq(len):
	""" returns random sequence of length len """
	return [random.choice("ACGT") for x in xrange(len)]

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

def mutate(nuc, aE):
	""" returns mutation of given nucleotide nuc with mutation rate 10^aE """
	""" TODO: Ratenmatrix verwenden! """
	r = random.random() * pow(10, (-1) * aE)
	if r < 1:
		return "A" if nuc == "T" else ("C" if nuc == "A" else ("G" if nuc == "C" else "T"))
	elif r < 2:
		return "C" if nuc == "T" else ("G" if nuc == "A" else ("T" if nuc == "C" else "A"))
	elif r < 3:
		return "G" if nuc == "T" else ("T" if nuc == "A" else ("A" if nuc == "C" else "G"))
	else:
		return nuc

def simulation(L, TE, dtE, aE, sE):
	"""
	simulation(100, 9, 5, -7, 6)
	"""
	seq = rseq(L)
	mut = seq
	data = []
	for step in xrange(pow(10, TE - sE)):
		for m in xrange(pow(10, sE - dtE)):
			mut = map(lambda s: mutate(s, aE), mut)
		d = dist(seq, mut)
		e = jcest(d, L)
		data.append((step, d, e));
	return data

def latex(data):
	xscale = 1
	yscale = 5
	print "\\begin{picture}(" + str(xscale * len(data)) + "," + str(yscale * L) + ")"
	print "\\put(0,0){\\vector(1,0){" + str(xscale * len(data)) + "}}"
	print "\\put(0,0){\\vector(0,1){" + str(yscale * L) + "}}"
	for step, dist, jcest in data:
		print "\\put(" + str(xscale * step) + "," + str(yscale * dist) + "){\\circle*{0.1}}"
		print "\\put(" + str(xscale * step) + "," + str(yscale * 10 * jcest) + "){\\circle*{0.1}}"
	print "\\end{picture}"



parser = argparse.ArgumentParser(description="Uebung 2 - Aufgabe 7")

parser.add_argument('-L', '--length',
type=int,
required='True',
help='Laenge der DNA-Sequenz')

parser.add_argument('-TE', '--time',
type=int,
required='True',
help='zu simulierender Zeitraum in 10^TE Jahren')

parser.add_argument('-dtE', '--delta',
type=int,
required='True',
help='Zeitraum fuer eine Mutation in 10^dtE Jahren')

parser.add_argument('-aE', '--alpha',
type=int,
required='True',
help='Mutationsrate in 10^aE Jahren')

parser.add_argument('-sE', '--step',
type=int,
required='True',
help='Schrittweite fuer Zwischenergebnisse in 10^stepE Jahren')

args = parser.parse_args()

L = args.length
TE = args.time
dtE = args.delta
aE = args.alpha
sE = args.step

data = simulation(L, TE, dtE, aE, sE)
latex(data)



"""
sys.stdout.write(", ".join(str(dist) for step, mut, dist in data))
sys.stdout.write(os.linesep)
"""















