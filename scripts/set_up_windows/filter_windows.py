import os
import sys

chr_names = {str(i):1 for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,"X","Y","M","MT"]}
with open(sys.argv[1], "r") as input:
    line = input.readline()
    print line.strip()
    line = input.readline()
    while(line):
        row = line.split("\t")
        if row[2] in chr_names:
            print line.strip()
        line = input.readline()