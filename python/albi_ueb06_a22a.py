import random
import os

N = 10000
LEN = 1000

f = open('library.fa', 'w')

nukl = ['A','C','G','T']


for i in range(N):
    f.write('> random'+str(i)+os.linesep)
    for j in range(LEN):
        f.write(random.choice(nukl))
    f.write(os.linesep)

f.close()
        
