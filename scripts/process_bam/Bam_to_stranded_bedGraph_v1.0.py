import sys,os
import pysam
import argparse
import collections
from collections import defaultdict
import numpy as np
import time
from time import gmtime, strftime

def base_to_hist(data,chr,output):
	''' This function convert the single base BED to a histogram BED to save space '''
	for pos,coverage in data.iteritems(): #input pos should be 0-based
		try:
			if pos == last_pos_1:
				if coverage == last_coverage:
					last_pos_1 = pos + 1
				else:
					output.write("\t".join([str(chr),str(last_pos),str(last_pos_1),str(last_coverage)]))
					output.write("\n")
					last_pos = pos  #0-based
					last_pos_1 = pos + 1
					last_coverage = coverage
			else:
				output.write("\t".join([str(chr),str(last_pos),str(last_pos_1),str(last_coverage)]))
				output.write("\n")
				last_pos = pos #0-based
				last_pos_1 = pos + 1
				last_coverage = coverage
		except NameError: #first one
			last_pos = pos  #0-based
			last_pos_1 = pos + 1
			last_coverage = coverage
	output.write("\t".join([str(chr),str(last_pos),str(last_pos_1),str(last_coverage)])) #last one
	output.write("\n")
	
def write_output(chr,plus,minus):
	not_empty = False
	dir_use = collections.OrderedDict(sorted(blocks_dict['+'].items()))
	if dir_use:
		not_empty = True
		#Output each base
		# for pos,coverage in dir_use.iteritems():
			# plus.write("\t".join([str(chr),str(pos),str(pos+1),str(coverage)]))
			# plus.write("\n")
		#Output as hist	
		base_to_hist(dir_use,chr,plus)
		
	dir_use = collections.OrderedDict(sorted(blocks_dict['-'].items()))
	if dir_use:
		not_empty = True
		#Output each base
		# for pos,coverage in dir_use.iteritems():
			# minus.write("\t".join([str(chr),str(pos),str(pos+1),str(coverage)]))
			# minus.write("\n")
		#Output as hist	
		base_to_hist(dir_use,chr,minus)
		
	if not_empty:
		sys.stderr.write("[%s] Finished: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()), str(chr)))
		
def blocks_to_dict(chr,strand,blocks):
	for block in blocks:
		for i in xrange(block[0],block[1]+1):
			blocks_dict[strand][i] += 1

def read_windows():
	''' read sliding windows file '''
	win_dict = {} #chr --> strand --> interval --> windows
	with open(options.windows) as windows:
		line = windows.readline()
		if line.startswith("#"):
			line = windows.readline()
		while(line):
			line = line.strip().split("\t")
			name = line[1]
			chr = line[2]
			strand = linep[3]
			start = line[9]
			end = line[10]		
			line = windows.readline()
	return win_dict
	
def union_interval(input):
	tmp = []
	for begin,end in sorted(input):
		if tmp and tmp[-1][1] >= begin - 1:
			tmp[-1][1] = max(tmp[-1][1], end)
		else:
			tmp.append([begin, end])
	return tmp

def map_length(input):
	sum = 0
	for i in input:
		sum += abs(i[1]-i[0])
	return sum
	
def get_read_stack_single(read1):
	global read_lengths
	blocks = read.get_blocks() #in genomic coordinates
	if read1.is_reverse:
		strand = "-"
	else:
		strand = "+"
	length = map_length(blocks)
	read_lengths.append(length)
	#print 1,blocks,map_length(blocks)
	chr = read1.reference_name
	blocks_to_dict(chr,strand,blocks)
	
def get_read_stack_mate(read1,read2):
	global read_lengths
	if read1.is_read2:
		read1,read2 = read2,read1
	if read1.reference_name != read2.reference_name: #check point, read1 and read2 should align to same conting
		return None
	blocks_1 = read1.get_blocks()
	blocks_2 = read2.get_blocks()
	blocks = union_interval(blocks_1+blocks_2)
	if read1.is_reverse:
		strand = "-"
	else:
		strand = "+"
	length = map_length(blocks)
	read_lengths.append(length)
	#print 2,blocks,map_length(blocks)	
	chr = read1.reference_name
	blocks_to_dict(chr,strand,blocks)
	#DEBUG OUTPUT
	#print "#1",blocks_1,"#2",blocks_2,"#merge",blocks,"#length",map_length(blocks)
	
	
def get_rpkm():
	pass

if __name__ == "__main__":
	description = """
	To reduce memory usage, please use query name coordinates;
	Make sure that read1 and read2 have same names;
	For fwd-rev stranded library only;
	Note: Input bam should be sorted, and have an index file. (For tophat2 output, you just need to index the bam)
	"""
	parser = argparse.ArgumentParser(prog="Bam_to_stranded_bedGraph",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_require = parser.add_argument_group("Required")
	group_require.add_argument("-i","--input",dest="input",required=True,help="")
	group_require.add_argument("-n","--name",dest="name",required=True,help="Name prefix. File name = [name].plus.bedGraph & [name].minus.bedGraph")
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("-r","--reverse",dest="reverse",action="store_true",default=False,help="Reverse read1 and read2. If use dUTP library")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	sys.stderr.write("[%s] Work begins...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	
	#windows = read_windows()
	
	TMP = {}
	read_lengths = []
	
	with pysam.AlignmentFile(options.input,"rb") as input, file(options.name+".plus.bedGraph",'w') as plus, file(options.name+".minus.bedGraph",'w') as minus:
		#Set default values
		contings = input.header['SQ']
		for conting in contings:
			blocks_dict = {'+':{},'-':{}}
			blocks_dict['+'] = defaultdict(int)
			blocks_dict['-'] = defaultdict(int)

			for read in input.fetch(conting['SN']):
				if read.is_secondary == False: #primary alignamnet for a pair only
					if read.query_name not in TMP:
						if read.mate_is_unmapped or not read.is_paired:
							get_read_stack_single(read)
						TMP[read.query_name] = read #place the read into tmp (one mate)
					else:
						read2 = TMP.get(read.query_name)
						get_read_stack_mate(read,read2)
						TMP.pop(read.query_name)
			if options.reverse:
				write_output(conting['SN'],minus,plus)
			else:
				write_output(conting['SN'],plus,minus)
	sys.stderr.write("[%s] All finished!\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	
	#Print mean map lengths
	read_lengths = np.array(read_lengths)
	#print read_lengths.mean(),read_lengths.std(),read_lengths.max(),read_lengths.min()
	sys.stdout.write("\n")
	sys.stdout.write("Mean map length is: %.2f +/- %.2f bp\nMedian length: %.2f bp\nMax length: %.2f bp\nMin length: %.2f bp\n" % (read_lengths.mean(),read_lengths.std(),np.median(read_lengths),read_lengths.max(),read_lengths.min()))
	
	
