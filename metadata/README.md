For hg38 (Ensembl GRCh38 release 104).

Sliding windows: 50 nt (25 nt overlaps)

Please download the genome and annotation files from Ensembl:

http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/

http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/

and build the index and window file by yourself.

To set up new window file (scripts can be found in /scripts/set_up_windows/):

`python get_non_redundent_region_from_ucsc_format_Ens_gene_with_strand.py <UCSC table-like> <non-redundent exons>`

where `<UCSC table-like>` is a table as UCSC table browser output. If you only have `GTF`, you can try [this script](https://github.com/jhfoxliu/RNA-m5C/blob/master/0_m5C_step-by-step_metadata/gtf2anno.py) or [this script](https://github.com/jhfoxliu/RNA-m5C/blob/master/0_m5C_step-by-step_metadata/gtf2anno_plus.py) to convert your file. `<non-redundent exons>` is the output name of this step.

`python get_sliding_win_on_Ens_merged_tx.py <UCSC table-like> <chr sizes> <non-redundent exons> <output windows>`

where `<chr sizes>` should in format as `Homo_sapiens.GRCh38.dna_sm.primary_assembly.chr.size`. `<output windows>` is the name of sliding window output.

To change the window length, please modify `window_width` and `half_window` in `get_sliding_win_on_Ens_merged_tx.py`.


