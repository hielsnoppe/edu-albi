#This software is a free software. 
#Thus, it is licensed under GNU General Public License.
#Python implementation to Nussinov Folding Algorithm 
#for Homework 7 of Bioinformatics class.
#Forrest Bao, Nov. 29 <http://fsbao.net> <forrest.bao aT gmail.com>

import sys,string
from numpy import *
from matplotlib import *

#read sequences
#sequences are stored in many lines
f=open(sys.argv[1], 'r')
seq=[]
stringSeq = ""

for line in f:
	# meta daten des fasta formats skippen
	if not line.startswith('>'):
		stringSeq += string.strip(line)
	# else: TODO: neue Sequenz zu seq hinzufuegen

seq.append(stringSeq)

seq = map(lambda x: str.replace(x, "T", "U"), seq)

def delta(l,m):
	delta=0;
	if l=='A' and m=='U':
		return 2;
	elif l=='U' and m=='A':
		return 2;
	elif l=='G' and m=='C':
		return 3;
	elif l=='C' and m=='G':
		return 3;
	#elif l=='G' and m=='U':
	#	return 1;
	#elif l=='U' and m=='G':
	#	return 1;	
	else:
		return 0;

def buildDP(seq):
	L=len(seq);
	s=zeros((L,L));
	for n in xrange(1,L):
		for j in xrange(n,L):
			i=j-n;
			case1=s[i+1,j-1]+delta(seq[i],seq[j]);
			case2=s[i+1,j];
			case3=s[i,j-1];
			if i+3<=j:
				tmp=[];
				for k in xrange(i+1,j):
					tmp.append(s[i,k]+s[k+1,j]);
				case4=max(tmp);
				s[i,j]=max(case1,case2,case3,case4);
			else:
				s[i,j]=max(case1,case2,case3);
	return s;

def traceback(s,seq,i,j,pair):
	if i<j:
		if s[i,j]==s[i+1,j]:
			traceback(s,seq,i+1,j,pair);
		elif s[i,j]==s[i,j-1]:
			traceback(s,seq,i,j-1,pair);
		# auf hairpin loops > 3 achten!
		elif s[i,j]==s[i+1,j-1]+delta(seq[i],seq[j]):# and j-i > 3:
			pair.append([i,j,str(seq[i]),str(seq[j])]);
			traceback(s,seq,i+1,j-1,pair);
		else:
			for k in xrange(i+1,j):
				if s[i,j]==s[i,k]+s[k+1,j]:
					traceback(s,seq,i,k,pair);
					traceback(s,seq,k+1,j,pair);
					break;
	return pair;

for q in xrange(0,len(seq)):
	pair=traceback(buildDP(seq[q]),seq[q],0,len(seq[q])-1,[])
 	print "max # of folding pairs: ",len(pair);
	#for x in xrange(0,len(pair)):
	#	print '%d %d %s==%s' % (pair[x][0],pair[x][1],pair[x][2],pair[x][3]);
	print "---";
	
	# Vienna Darstellung
	braketStr = ''
	for seqPos in xrange(0,len(stringSeq)):
		ch = '.'
		for x in xrange(0,len(pair)):
			if seqPos == pair[x][0]:
				ch = '('
			elif seqPos == pair[x][1]:
				ch = ')'
		braketStr += ch
	
	for seqPos in xrange(0,len(stringSeq), 30):
		print '%d\t%s' % (seqPos, stringSeq[seqPos:(seqPos+30)])
		print '\t%s' % (braketStr[seqPos:(seqPos+30)])
