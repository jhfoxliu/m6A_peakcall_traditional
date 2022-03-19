
import sys

"""
argv:
    1, input file:  uscs gene annotation. for example ""
    2, output file: all non redundant exon regions. If exons have overlaps, I merge it into one exons.
"""

output=open(sys.argv[2],'w')
exon={}
input=open(sys.argv[1])
for line in input:
    if (line[0]=='#'):
        continue
    line=line[:-1]
    a=line.split('\t')
    start=a[9].split(',')[:-1]
    end=a[10].split(',')[:-1]
    chr=a[2]
    strand=a[3]
    ensg=a[12]
    for i in range(len(start)):
        if (not exon.has_key(ensg)):
            exon[ensg]=[[int(start[i]), int(end[i]), a[1], chr, strand]]
        else:
            exon[ensg].append([int(start[i]), int(end[i]), a[1], chr, strand])    # a[2] not useful for this programme

for key in exon.keys():
    chr=exon[key][0][3]
    strand=exon[key][0][4]
    exon[key].sort(key=lambda x:x[0])
    i=0
    while (i<len(exon[key])-1):
        if (exon[key][i][1]>exon[key][i+1][0] and exon[key][i+1][1]>exon[key][i][0] ):
            exon[key][i][0]=min(exon[key][i][0], exon[key][i+1][0])
            exon[key][i][1]=max(exon[key][i][1], exon[key][i+1][1])
            exon[key].pop(i+1)
        else:
            i+=1
    start=[]
    end=[]
    for i in range(len(exon[key])):
        start.append(str(exon[key][i][0]))
        end.append(str(exon[key][i][1]))
    output.write(key+'\t'+chr+'\t'+','.join(start+[''])+'\t'+','.join(end+[''])+'\t'+strand+'\n')
    
    
    
    
    


    
    
    
    
    
