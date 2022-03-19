import sys,os
from sjm_tools import job,check_env

if not sys.argv[1]:
	print "Usage: python %s [Sample Sheet]" % (sys.argv[0])

windows_file = check_env("/share/public1/data/liujh/database/human/ensembl_release104_GRCh38/sliding_window/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.noheader.exons.w50s25.chr.txt")
python = check_env('/share/public/apps/bin/python')
bedtools = check_env('/share/public/apps/bedtools/2.29.2/bin/bedtools')
perl = check_env('/share/public1/data/liujh/software/perlbrew/perls/perl-5.24.1/bin/perl')
chr_size = "/share/public1/data/liujh/database/human/ensembl_release104_GRCh38/genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chr.size"
script_prefix = '/share/public1/data/liujh/script/m5C-seq/'

PATH = "./"
with open(sys.argv[1], "r") as input:
    for line in input.readlines():
        if not line:
            continue
        name, input_prefix, input_winsize, ip_prefix, ip_winsize = line.strip().split("\t")
        
        SJM = name + "_call_peak.sjm"
        workpath = PATH + "/" + name + "_call_peak/"

        if os.path.isdir(workpath) == False:
            os.mkdir(workpath)
            
        JOB = job(workpath,SJM) 
        
        JOB.step_start(step_name="CallPeak",memory="50G")
        JOB.add_process("{perl} {script_prefix}/rnaexp_rpkm_strand_ZH.pl {windows_file} {chr_size} ../{prefix}/{prefix}.minus.bedGraph ../{prefix}/{prefix}.plus.bedGraph {prefix}_window_rpkm.txt {winsize}".format(script_prefix=script_prefix,perl=perl, chr_size=chr_size, windows_file=windows_file, winsize=input_winsize, prefix=input_prefix))
        JOB.add_process("{perl} {script_prefix}/rnaexp_rpkm_strand_ZH.pl {windows_file} {chr_size} ../{prefix}/{prefix}.minus.bedGraph ../{prefix}/{prefix}.plus.bedGraph {prefix}_window_rpkm.txt {winsize}".format(script_prefix=script_prefix,perl=perl, chr_size=chr_size, windows_file=windows_file, winsize=ip_winsize, prefix=ip_prefix))

        JOB.add_process("{python} {script_prefix}/merge_input_eluate_win_RPKM.py {prefix_input}_window_rpkm.txt {prefix_ip}_window_rpkm.txt {name}_eluate_win_rpkm.txt".format(script_prefix=script_prefix,python=python, prefix_input=input_prefix, prefix_ip=ip_prefix, name=name))

        JOB.add_process("{perl} {script_prefix}/calculate_winscore.pl {name}_eluate_win_rpkm.txt".format(script_prefix=script_prefix,perl=perl, name=name))

        JOB.add_process("{perl} {script_prefix}/remove_data_RPKM1.pl {name}_eluate_win_rpkm_winscore.txt {name}_winscore_filtered.txt".format(script_prefix=script_prefix,perl=perl, name=name))

        JOB.add_process("{perl} {script_prefix}/match_two_files.pl {name}_winscore_filtered.txt {windows_file} {name}_winscore_filtered_matched.txt".format(script_prefix=script_prefix,perl=perl, name=name, windows_file=windows_file))

        JOB.add_process("cut -f 3,5,6,17,18 {name}_winscore_filtered_matched.txt > {name}_winscore_filtered_matched_formated.txt".format(script_prefix=script_prefix,name=name))

        JOB.add_process("{python} {script_prefix}/m6A_add_strand.py {name}_winscore_filtered_matched_formated.txt > {name}_winscore_filtered_matched_formated.bed".format(script_prefix=script_prefix,python=python, name=name))

        JOB.add_process("{perl} {script_prefix}/combine_peak.0.5.pl {name}_winscore_filtered_matched_formated.txt".format(script_prefix=script_prefix,perl=perl, name=name))

        JOB.add_process("{perl} {script_prefix}/combine_peak2.pl {name}_winscore_filtered_matched_formated.txt".format(script_prefix=script_prefix,perl=perl, name=name))

        JOB.add_process("grep -P 'Y$' {name}_winscore_filtered_matched_formated.txt.co.pk | cut -f 1-5 > {name}_winscore_filtered_matched_formated_greped.txt".format(script_prefix=script_prefix,name=name))

        JOB.add_process("cat {name}_winscore_filtered_matched_formated_greped.txt {name}_winscore_filtered_matched_formated.txt.co.pk2 > {name}_peaklist.txt".format(script_prefix=script_prefix,name=name))

        JOB.add_process("{python} {script_prefix}/add_strand_to_final_peaks.py {windows_file} {name}_peaklist.txt {name}_final_peaks_with_strand.bed".format(script_prefix=script_prefix,python=python, windows_file=windows_file, name=name))
        
        JOB.step_end()
        JOB.job_finish()
        JOB.submit()
