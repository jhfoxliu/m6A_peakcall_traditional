import sys


data = {}
with open(sys.argv[1],'r') as input:
	for line in input.readlines():
		line = line.strip().split("\t")
		chr = line[0]
		pos = line[1]
		strand = line[2]
		cov = line[3]
		level = line[4]
		peak_start = line[5]
		peak_end = line[6]
		peak = float(line[7])
		
		key = "@".join([chr,pos,strand,cov,level])
		if key not in data:
			data[key] = [peak_start,peak_end,peak]
		else:
			if data[key][2] < peak:
				data[key] = [peak_start,peak_end,peak]

for key,values in data.items():
	print key.replace("@","\t")+ "\t" + "\t".join([str(i) for i in values])