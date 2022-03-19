#!/usr/bin/env python
import sys

strand_info={}
for line in open(sys.argv[1]):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    strand_info[a[1]]=a[3]

output=open(sys.argv[3], 'w')
for line in open(sys.argv[2]):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    strand=strand_info[a[3]]
    output.write(line+'\t'+strand+'\n')
    
