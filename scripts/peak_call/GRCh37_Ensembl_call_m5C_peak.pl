#!/usr/bin/perl
############################################################
#Author: Jianheng Liu
#Usage: call m5C-seq peaks, for pair-end reads
use warnings;
use strict;

die "Usage: perl script.pl input.bam IP.bam name readlength" if (!$ARGV[3]);

my $input = $ARGV[0];
$input=~ s/\.bam//i;
my $IP = $ARGV[1];
$IP =~ s/\.bam//i;
my $name = $ARGV[2];
my $readlength = $ARGV[3];

#Directory of some softwares and scripts
#my $myPerl = '/share/public1/data/liujh/software/perlbrew/perls/perl-5.24.1/bin/perl';
my $bedtools = '/share/public/apps/bin/bedtools';
my $python = '/share/public/apps/bin/python';
#Directory of the files:
my $windows = "/share/public1/data/liujh/database/m5C_artificial/sliding_window_on_merged_ensembl.artificial.txt";
my $size = "/share/public1/data/liujh/database/m5C_artificial/chr_size.txt";
#Validate directories:

my $script = "$name.m5C_seq_call_peak.sh";

open (my $OUTPUT, ">", $script);
	print $OUTPUT "$bedtools genomecov -bg -split -pc -strand + -ibam $input.bam > $input.plus.bedGraph;\n";
	print $OUTPUT "$bedtools genomecov -bg -split -pc -strand + -ibam $IP.bam > $IP.plus.bedGraph;\n";
	print $OUTPUT "$bedtools genomecov -bg -split -pc -strand - -ibam $input.bam > $input.minus.bedGraph;\n";
	print $OUTPUT "$bedtools genomecov -bg -split -pc -strand - -ibam $IP.bam > $IP.minus.bedGraph;\n";
	print $OUTPUT "perl /share/public1/data/liujh/script/m5C-seq/rnaexp_rpkm_strand_ZH.pl $windows $size $input.minus.bedGraph $input.plus.bedGraph $input","_window_rpkm.txt $readlength;\n";
	print $OUTPUT "perl /share/public1/data/liujh/script/m5C-seq/rnaexp_rpkm_strand_ZH.pl $windows $size $IP.minus.bedGraph $IP.plus.bedGraph $IP","_window_rpkm.txt $readlength;\n";
	print $OUTPUT "$python /share/public1/data/liujh/script/m5C-seq/merge_input_eluate_win_RPKM.py $input","_window_rpkm.txt $IP","_window_rpkm.txt $name","_eluate_win_rpkm.txt;\n";
	print $OUTPUT "perl /share/public1/data/liujh/script/m5C-seq/calculate_winscore.pl $name","_eluate_win_rpkm.txt;\n";
	print $OUTPUT "perl /share/public1/data/liujh/script/m5C-seq/remove_data_RPKM1.pl $name","_eluate_win_rpkm_winscore.txt $name","_winscore_filtered.txt;\n";
	print $OUTPUT "perl /share/public1/data/liujh/script/m5C-seq/match_two_files.pl $name","_winscore_filtered.txt $windows $name","_winscore_filtered_matched.txt;\n";
	print $OUTPUT "cut -f 3,5,6,17,18 $name","_winscore_filtered_matched.txt > $name","_winscore_filtered_matched_formated.txt;\n";
	print $OUTPUT "perl /share/public1/data/liujh/script/m5C-seq/combine_peak.0.5.pl $name","_winscore_filtered_matched_formated.txt;\n";
	print $OUTPUT "perl /share/public1/data/liujh/script/m5C-seq/combine_peak2.pl"," ","name","_winscore_filtered_matched_formated.txt;\n";	
	print $OUTPUT "grep -P \'Y\$\' $name","_winscore_filtered_matched_formated.txt.co.pk | cut -f 1-5 > $name","_winscore_filtered_matched_formated_greped.txt;\n";
	print $OUTPUT "cat $name","_winscore_filtered_matched_formated_greped.txt $name","_winscore_filtered_matched_formated.txt.co.pk2  > $name","_peaklist.txt;\n";
	print $OUTPUT "$python /share/public1/data/liujh/script/m5C-seq/add_strand_to_final_peaks.py $windows $name","_peaklist.txt $name","_final_peaks_with_strand.bed;\n";