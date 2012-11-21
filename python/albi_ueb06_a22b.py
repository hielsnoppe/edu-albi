import random
import os

K = 100

f = open('query.fa', 'w')

nukl = ['A','C','G','T']

f.write('> myrandomnukl'+os.linesep)

for i in range(K):
    f.write(random.choice(nukl))
    
f.close()
        
