bowtie -v 3 -p 8 --sam /compbio/data/bowtieIndex/hg19 lncap.foax1.export.txt.gz.fq >foxa1.sam
epd_python get_flagstat_summary.py foxa1.sam >summary.txt
