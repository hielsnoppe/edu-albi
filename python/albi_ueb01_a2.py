import sys
import re
import os

regex = r"(^>\S*\s)"

in_file = open(sys.argv[1], "rt")
out_file = open(sys.argv[2], "wt")

for line in in_file:
    c = re.match(regex, line)
    if (c):
        out_file.write(c.group(1))
        out_file.write(os.linesep)
    else:
        out_file.write(line)
		
in_file.close()
out_file.close()
