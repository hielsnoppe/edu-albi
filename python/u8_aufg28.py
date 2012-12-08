#!/usr/bin/env python

import math
import csv
import sys

def H(data):
  if not data:
    return 0
  entropy = []
  
  total = []
  for i in range(len(data)):
  	total.append(sum(data[i]))

  for k in range(len(data[0])):
  	entropy_k = 0
  	for x in range(len(data)):
  		p_x = float(data[x][k])/float(total[x])
    	if p_x > 0.0:
      		entropy_k += - p_x*math.log(p_x, 2)
  	entropy.append(entropy_k)  	
  return entropy

data = []
reader = csv.reader(open(sys.argv[1], "rb")) 
for row in reader:
  row = list(int(elem) for elem in row)
  data.append(row)


entropy = H(data)
print entropy
print sum(entropy)
