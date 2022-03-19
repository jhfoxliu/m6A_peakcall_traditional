import sys

fn = sys.argv[1] #make sure this file is sorted by chr position
step = 25

data = {}
with open(fn,'r') as input:
	for line in input.readlines():
		line = line.strip().split("\t")
		chr,start,end,window,peak = line
		start = int(start)
		half = start + step
		end = int(end)
		peak = float(peak)
		window = "_".join([i for i in window.split("_")[0:3]])
		key1 = "@".join([chr,window,str(start),str(half)])
		key2 = "@".join([chr,window,str(half),str(end)])
		data[key1] = peak if key1 not in data else (data[key1]+peak)/2
		data[key2] = peak if key2 not in data else (data[key2]+peak)/2

for key,values in sorted(data.items()):
	key= key.split("@")
	print "\t".join([key[0],key[2],key[3],key[1],str(values)])
