import string
import math
import itertools
import sys, os

def fastaread(fname):
	seqs = []
	seq = ""
	with open(fname, "r") as fh:
		for line in fh:
			if line.startswith('>'):
				seqs.append(seq)
				seq = ""
			else:
				seq += string.strip(line)
	return seqs

def train(seqs):
	nucs = ["a", "c", "g", "t"]
	binucs = list(itertools.product(nucs, repeat=2))
	trinucs = list(itertools.product(nucs, repeat=3))

	C = [dict([("".join(trinuc), 1.0) for trinuc in trinucs]) for i in range(3)]
	A = [dict([("".join(binuc), dict([("".join(binuc), 0.0) for binuc in binucs])) for binuc in binucs]) for i in range(3)]

	for seq in seqs:
		step = 0
		cod = seq[:3]			# current codon

		for nuc in seq:
			rfn = step % 3		# reading frame number is 0, 1 or 2
			C[rfn][cod] = C[rfn][cod] + 1
			step = step + 1;
			cod = cod[1:] + nuc

	for i in range(len(C)):
		for trinuc in C[i]:
			rs, st = trinuc[:2], trinuc[1:]
			A[i][rs][st] = C[i][trinuc] / sum([C[i][rs+l] for l in nucs])

	return (A[0], A[1], A[2])

def Pr(seq, A):
	return sum([math.log(A[seq[n:n+2]][seq[n+1:n+3]]) for n in range(0, len(seq)-3)])

def analyze(seq, Gmod, NCmod, winsize):
	data = []
	S = [0.0, 0.0, 0.0]

	step = 0
	win = seq[:winsize]		# current window

	for n in seq[winsize:]:
		rfn = step % 3		# reading frame number is 0, 1 or 2
		S[rfn] = Pr(win, Gmod) - Pr(win, NCmod)

		if rfn == 2:
			data.append((step, S[0], S[1], S[2]))

		step = step + 1;
		win = win[1:] + n
	return data

(G1, G2, G3) = train(fastaread("../data/Genes.fa"))
(NC1, NC2, NC3) = train(fastaread("../data/Noncoding.fa"))
sequence = fastaread("../data/Test.fa")[0]
data = analyze(sequence, G1, NC1, 100)
data = analyze(sequence, G2, NC2, 100)
data = analyze(sequence, G3, NC3, 100)

"""
def latex_dict(d):
	print "\\begin{tabular}{c|" + "".join(["c" for j in range(16)]) + "}"
	sys.stdout.write("   & " + " & ".join([h for h in d.keys()]) + "\\\\\\hline" + os.linesep)
	for row in d.iterkeys():
		sys.stdout.write(row + " & ")
		sys.stdout.write(" & ".join(["%3.0f" % (f*1000) for f in d[row].values()]) + "\\\\" + os.linesep)
	print "\\end{tabular}"

def latex_plot(data):
	xscale = 1.0/3.0
	yscale = 5
	ymax = 15
	print "\\begin{picture}(" + str(xscale * len(data)) + "," + str(yscale * 2 * ymax) + ")(0,-" + str(yscale * ymax) + ")"
	print "\\put(0,0){\\vector(1,0){" + str(xscale * len(data) * 3) + "}}"
	print "\\put(0,0){\\vector(0,1){" + str(yscale * ymax) + "}}"
	print "\\put(0,0){\\vector(0,-1){" + str(yscale * ymax) + "}}"
	for step, S1, S2, S3 in data:
		print "\color{red}"
		print "\\put(" + str(xscale * step) + "," + str(yscale * S1) + "){\\circle*{0.1}}"
		print "\color{green}"
		print "\\put(" + str(xscale * step) + "," + str(yscale * S2) + "){\\circle*{0.1}}"
		print "\color{blue}"
		print "\\put(" + str(xscale * step) + "," + str(yscale * S3) + "){\\circle*{0.1}}"
	print "\\end{picture}"
"""
