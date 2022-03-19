# Traditional way to call modification RIP-seq peaks via sliding windows

Can be used in m6A-seq, m6Am-seq, m5C-seq, etc.

Usage:

If you are using SJM:

firstly, prepare a Sample Sheet (`Samples.txt`) like this (reads can be fastq or fastq.gz):

Single-end:

`Sample[\t]read.fastq`

Paired-end:

`Sample[\t]read1.fastq[\t]read2.fastq`

Then run the first pipeline:

`python hisat2_hg38_BedGraph[xxx].py Samples.txt`

You will get the bedGraph files. Please check the logs to get the median lengths.

When this is finished, prepare another Sample Sheet for peak call (`Peakcall.txt`), should be like this:

`Experiment[\t]Input name[\t]Input median length[\t]IP name[\t]IP median length`

Then run 

`python m6A_peak_call_hg38.py Peakcall.txt`

to get the winscores and peaks.