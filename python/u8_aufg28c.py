#!/usr/bin/env python
import math

# Spalte 1,2,3
f_i = [	[2,2,2], # A
        [2,2,1], # C
        [2,2,4], # G
        [2,2,1]] # U
f_ij = [[0,1,0], # AA
        [0,0,1], # AC
        [0,1,1], # AG
        [2,0,0], # AU
        [0,0,1], # CA
        [0,0,0], # ...
        [2,2,0],
        [0,0,1],
        [0,1,1], # GA
        [2,0,0],
        [0,0,2],
        [0,1,0],
        [2,0,1], # UA
        [0,1,0],
        [0,1,1],
        [0,0,0]] # UU

H = 0

for i in range(4):
  for j in range(4): 
    for spalte in range(3):
      p = f_ij[i*4+j][spalte]
      a = f_i[i][spalte]
      b = f_i[j][spalte]
      
      if p == 0:
        p = float(1)/16 # pseudocount
      else:
        p = float(p)/8
      
      a = float(a)/8
      b = float(b)/8
      
      H += p * math.log((p/(a*b)), 2)
      
print (-H) # -4.9375
