# hg38
import readline
import sys

genelist = "/share/public1/data/liujh/database/human/ensembl_release104_GRCh38/sliding_window/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.noheader.exons.w50s25.chr.txt"

data = {}
with open(genelist, "r") as input:
    for line in input.readlines():
        line = line.strip().split("\t")
        gene = line[1]
        strand = line[3]
        data[gene] = strand
with open(sys.argv[1], "r") as input:
    for line in input.readlines():
        if not line:
            continue    
        line = line.strip().split("\t")
        win = line[3]
        chr = line[0]
        start = line[1]
        end = str(int(line[2]) + 1)
        peak = line[4]
        strand = data[win]
        print "{}\t{}\t{}\t{}\t{}\t{}".format(chr, start, end, win, peak, strand)
        
