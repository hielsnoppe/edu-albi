#!/usr/bin/env python

import math
import csv
import sys

data = []
reader = csv.reader(open(sys.argv[1], "rb")) # TATA Box einlesen
for row in reader:
  row = list(int(elem) for elem in row)
  data.append(row)
  
multi_align = []
for x in range(389):
  multi_align.append("")

nukleotid = ['A','C','G','T']

for pos in range(len(data[0])):
  seq_index = 0
  for nuk in range(4): 
    val = data[nuk][pos]
    for i in range(val):
      multi_align[seq_index] += nukleotid[nuk]
      seq_index += 1
      
for seq in multi_align:
  print seq
