import math

sequence = "ATGAATACTTTACGTATTGGCTTAGTTTCCATCTCATACTAATAAGAAGGACAAGACGGGGGTAAGGAGA\
CCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCC\
TATTGGTCTATTTTCCCACCCTTAGATGATCGCGCATCCAGCGGCGTTTATCAGGATAAAGGCATCCCTG\
CGCTGGAAGAATGGCTGACATCGGCGCTAACCACGCCGTTTGAACTGGAAACCCGCTTAATCCCCGATGA\
GCAGGCGATCATCGAGCAAACGTTGTGTGAGCTGGTGGATGAAATGAGTTGCCATCTGGTGCTCACCACG\
GGCGGAACTGGCCCGGCGCGTCGTGACGTAACGCCCGATGCGACGCTGGCAGTAGCGGACCGCGAGATGC\
CTGGCTTTGGTGAACAGATGCGCCAGATCAGCCTGCATTTTGTACCAACTGCGATCCTTTCGCGTCAGGT\
GGGCGTGATTCGCAAACAGGCGCTGATCCTTAACTTACCCGGTCAGCCGAAGTCTATTAAAGAGACGCTG\
GAAGGTGTGAAGGACGCTGAGGGTAACGTTGTGGTACACGGTATTTTTGCCAGCGTACCGTACTGCATTC\
AGTTGCTGGAAGGGCCATACGTTGAAACGGCACCGGAAGTGGTTGCAGCATTCAGACCGAAGAGTGCAAG\
ACGCGACGTTAGCGAATAA";

FF = {
"UUU": 19.7,
"UCU": 5.7,
"UAU": 16.8,
"UGU": 5.9,
"UUC": 15.0,
"UCC": 5.5,
"UAC": 14.6,
"UGC": 8.0,
"UUA": 15.2,
"UCA": 7.8,
"UAA": 1.8,
"UGA": 1.0,
"UUG": 11.9,
"UCG": 8.0,
"UAG": 0.0,
"UGG": 10.7,
"CUU": 11.9,
"CCU": 8.4,
"CAU": 15.8,
"CGU": 21.1,
"CUC": 10.5,
"CCC": 6.4,
"CAC": 13.1,
"CGC": 26.0,
"CUA": 5.3,
"CCA": 6.6,
"CAA": 12.1,
"CGA": 4.3,
"CUG": 46.9,
"CCG": 26.7,
"CAG": 27.7,
"CGG": 4.1,
"AUU": 30.5,
"ACU": 8.0,
"AAU": 21.9,
"AGU": 7.2,
"AUC": 18.2,
"ACC": 22.8,
"AAC": 24.4,
"AGC": 16.6,
"AUA": 3.7,
"ACA": 6.4,
"AAA": 33.2,
"AGA": 1.4,
"AUG": 24.8,
"ACG": 11.5,
"AAG": 12.1,
"AGG": 1.6,
"GUU": 16.8,
"GCU": 10.7,
"GAU": 37.9,
"GGU": 21.3,
"GUC": 11.7,
"GCC": 31.6,
"GAC": 20.5,
"GGC": 33.4,
"GUA": 11.5,
"GCA": 21.1,
"GAA": 43.7,
"GGA": 9.2,
"GUG": 26.4,
"GCG": 38.5,
"GAG": 18.4,
"GGG": 8.6
}

f = {
"TTT": 19.7,
"TCT": 5.7,
"TAT": 16.8,
"TGT": 5.9,
"TTC": 15.0,
"TCC": 5.5,
"TAC": 14.6,
"TGC": 8.0,
"TTA": 15.2,
"TCA": 7.8,
"TAA": 1.8,
"TGA": 1.0,
"TTG": 11.9,
"TCG": 8.0,
"TAG": 0.000001,	# must not be zero
"TGG": 10.7,
"CTT": 11.9,
"CCT": 8.4,
"CAT": 15.8,
"CGT": 21.1,
"CTC": 10.5,
"CCC": 6.4,
"CAC": 13.1,
"CGC": 26.0,
"CTA": 5.3,
"CCA": 6.6,
"CAA": 12.1,
"CGA": 4.3,
"CTG": 46.9,
"CCG": 26.7,
"CAG": 27.7,
"CGG": 4.1,
"ATT": 30.5,
"ACT": 8.0,
"AAT": 21.9,
"AGT": 7.2,
"ATC": 18.2,
"ACC": 22.8,
"AAC": 24.4,
"AGC": 16.6,
"ATA": 3.7,
"ACA": 6.4,
"AAA": 33.2,
"AGA": 1.4,
"ATG": 24.8,
"ACG": 11.5,
"AAG": 12.1,
"AGG": 1.6,
"GTT": 16.8,
"GCT": 10.7,
"GAT": 37.9,
"GGT": 21.3,
"GTC": 11.7,
"GCC": 31.6,
"GAC": 20.5,
"GGC": 33.4,
"GTA": 11.5,
"GCA": 21.1,
"GAA": 43.7,
"GGA": 9.2,
"GTG": 26.4,
"GCG": 38.5,
"GAG": 18.4,
"GGG": 8.6
}

def getH(seq):
	return sum([math.log(f[seq[n:n+3]]) for n in range(0, len(seq), 3)])

def main(seq, winsize):
	fsize = 3 * winsize
	H, P = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
	orf = seq[:fsize]
	step = 0
	
	for n in seq[:-fsize]:
		rfn = step % 3
		H[rfn] = H[rfn] + getH(orf)

		if rfn == 2:
			A = sum(H) / float(len(H))
			H = map(lambda Hi: math.exp(Hi - A), H)
			print(step, rfn, A, H)
			P = map(lambda Hi: Hi / sum(H), H)
			#print(math.log10(P / 1.0 - P))

		step = step + 1;
		orf = orf[1:] + n
	return

def main2(seq, winsize):
	H1, H2, H3 = 0, 0, 0
	for step in range(winsize):
		take_start, take_end = step * 3, (step + 1) * 3
		print(step, take_start, take_end)
		H1 = H1 + F[seq[take_start:take_end]]
		H2 = H2 + F[seq[take_start+1:take_end+1]]
		H3 = H3 + F[seq[take_start+2:take_end+2]]

	for step in range(winsize, len(seq) / 3):
		take_start, take_end = step * 3, (step + 1) * 3
		drop_start, drop_end = take_start - 3 * winsize, take_end - 3 * winsize
		print(step, drop_start, drop_end, take_start, take_end)

		H1 = H1 - F[seq[drop_start:drop_end]] + F[seq[take_start:take_end]]
		H2 = H2 - F[seq[drop_start+1:drop_end+1]] + F[seq[take_start+1:take_end+1]]
		H3 = H3 - F[seq[drop_start+2:drop_end+2]] + F[seq[take_start+2:take_end+2]]

		A = (H1 + H2 + H3) / 3.0
		h1, h2, h3 = math.exp(H1 - A), math.exp(H2 - A), math.exp(H3 - A)
		dvsr = h1 + h2 + h3
		P1, P2, P3 = h1 / dvsr, h2 / dvsr, h3 / dvsr

		print(P1, P2, P3)
		print(math.log10(P1 / 1.0 - P1))
		print(math.log10(P2 / 1.0 - P2))
		print(math.log10(P3 / 1.0 - P3))
	return

main(sequence, 30);













