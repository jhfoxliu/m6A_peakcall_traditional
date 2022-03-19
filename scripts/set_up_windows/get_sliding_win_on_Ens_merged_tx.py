#!/usr/bin/env python
import sys

window_width = 50 
half_window = 25

if len(sys.argv) < 5:
	print "Usage: python script.py anno size noredudant output"

gene_strand={}

for line in open(sys.argv[1]):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    genename=a[12]
    gene_strand[genename]=a[3]
    
hg19_size={}
for line in open(sys.argv[2]):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    hg19_size[a[0]]=int(a[1])

output=open(sys.argv[4], 'w')
output.write('\t'.join(["#bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds",
                        "score","name2","cdsStartStat","cdsEndStat","exonFrames"])+'\n')
num=1
for line in open(sys.argv[3],'r'):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    genename=a[0]
    strand=gene_strand[genename]
    chr=a[1]
    exonstart=a[2].split(',')[:-1]
    exonend=a[3].split(',')[:-1]
    for i in range(len(exonstart)):
        window=1
        exonid=i+1
        start=int(exonstart[i])
        end=int(exonend[i])
        if chr not in hg19_size:
            continue
        if (end+window_width> hg19_size[chr]):
            continue
        exonlen=end-start
        for j in range(exonlen/half_window):
            #output.write('\t'.join([chr, str(start+50*j), str(start+50*j+100), genename+'_'+str(window),  '1000', strand])+'\n')
            output.write('\t'.join([str(num), genename+'_exon_'+str(exonid)+'_'+str(window), chr,  strand, str(start+half_window*j), str(start+half_window*j+window_width), str(start+half_window*j), str(start+half_window*j+window_width),
                                    '1', str(start+half_window*j), str(start+half_window*j+window_width), '0', genename+'_exon_'+str(exonid)+'_'+str(window), 'cmpl', 'cmpl', '0'])+'\n')
            window+=1
            num+=1
            


        

            
