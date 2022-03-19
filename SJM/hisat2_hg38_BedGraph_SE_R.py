#!/usr/bin/env python
#Jianheng Liu, 20180709
#Usage: create SJM files

import sys,os
from tabnanny import check
from sjm_tools import job,check_env

if not sys.argv[1]:
	print "Usage: python %s [directory] [submit(optional)]" % (sys.argv[0])
	
#metadata
genome_index = check_env('/share/public1/data/liujh/database/human/ensembl_release104_GRCh38/hisat2_index/Homo_sapiens.GRCh38.dna_sm.primary_assembly.format',is_prefix=True)
# rRNA_index = check_env('/share/public1/data/liujh/database/human/rRNA_silva/silva_human_rRNA',is_prefix=True)
gtf = check_env('/share/public1/data/liujh/database/human/ensembl_release104_GRCh38/metadata/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.noheader.gtf')
# codon = check_env('/share/public1/data/liujh/database/fly/dm6/Drosophila_melanogaster.BDGP6.93.no_pseudogene.codon')
size = check_env('/share/public1/data/liujh/database/human/ensembl_release104_GRCh38/genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.size')
#environment
python = '/share/public/apps/bin/python'
samtools = '/share/public1/data/liujh/software/bin/samtools'
hisat2 = '/share/public/apps/hisat2/2.1.0/hisat2'
stringtie = "/share/public/apps/Stringtie/stringtie-2.0.3/stringtie"
cutadapt = '/share/public/apps/bin/cutadapt3'
trimmomatic = '/share/public1/data/liujh/software/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar'
java = '/share/public1/data/liujh/software/java/jre1.8.0_131/bin/java'
sortmerna = '/share/public/apps/sortmerna-2.1/sortmerna'
htseq = '/share/public/apps/bin/htseq-count'
samtools = '/share/public/apps/bin/samtools'
bowtie2 = '/share/public/apps/bowtie2/2.3.4.2/bowtie2'
cufflinks = '/share/public/apps/bin/cufflinks'
bedtools = '/share/public/apps/bedtools/2.29.2/bin/bedtools'
#namespace

PATH = "./"#sys.argv[1]
if os.path.isdir(PATH) == False:
	os.mkdir(PATH)

with open(sys.argv[1],"r") as input:
	for line in input.readlines():
		line = line.strip().split("\t")
		name, read1 = line
		
		SJM = name + ".IP.sjm"
		workpath = PATH + "/" + name + "/"
		JOB = job(workpath,SJM)
		
		if os.path.isdir(workpath) == False:
			os.mkdir(workpath)
				
		JOB.step_start(step_name="QC",memory="20G")
		JOB.add_process("{cutadapt} -j 6 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --max-n 0 -e 0.1 -q 20 -m 20 --trim-n -o read1.cutadapt.fastq ../{read1}".format(cutadapt=cutadapt,read1=read1))
		JOB.step_end()

		JOB.step_start(step_name="hisat2",memory="100G")
		JOB.add_process("{hisat2} -p 12 --no-unal  --rna-strandness R -x {genome_index} -U read1.cutadapt.fastq -S hisat2.sam".format(hisat2=hisat2,genome_index=genome_index))
		# unique only
		JOB.add_process("{samtools} view -q 60 -bS -@ 6 hisat2.sam > hisat2.bam ;\n".format(samtools=samtools))
		JOB.add_process("{samtools} sort -@ 4 -m 4G -o hisat2.sorted.bam hisat2.bam;\n".format(samtools=samtools))
		JOB.add_process("{samtools} index hisat2.sorted.bam ;\n".format(samtools=samtools))
		JOB.add_process("{samtools} flagstat hisat2.sorted.bam ;\n".format(samtools=samtools))
		JOB.step_end()
		
		JOB.step_start(step_name="BedGraph",memory="100G")
		JOB.add_process("{bedtools} genomecov -bg -split -strand + -ibam hisat2.sorted.bam -g {size} > {name}.plus.bedGraph".format(bedtools=bedtools, size=size, name=name))
		JOB.add_process("{bedtools} genomecov -bg -split -strand - -ibam hisat2.sorted.bam -g {size} > {name}.minus.bedGraph".format(bedtools=bedtools, size=size, name=name))

		JOB.step_end()
		
		JOB.job_finish()
		JOB.submit()
