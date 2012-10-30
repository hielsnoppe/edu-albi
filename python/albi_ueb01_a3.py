import re
import argparse

parser = argparse.ArgumentParser(description="Uebung 1 - Aufgabe 3")
parser.add_argument('-i',
'--inputfile',
type=str,
required='True',
help='Datei mit Homologen von FAL1')

parser.add_argument('-p',
'--properties',
type=str,
required='True',
help='Datei mit Proteineigenschaften')

args = parser.parse_args()

in_file = open(args.inputfile, "rt")
prop_file = open(args.properties, "rt")

DALTONS = 40000
homo_list = []


for line in in_file:
    line = line[:-1]
    homo_list.append(line)

for line in prop_file:
    line_array = re.split(r'\t', line)
    prot = line_array[0].strip()
    weight = int(line_array[2].strip())
    
    if(prot in homo_list and weight > DALTONS):
        print (prot + ": " + str(weight))
        
in_file.close()
prop_file.close()
