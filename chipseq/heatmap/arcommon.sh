for i in `cat AR.txt`
do
	epd_python hm_matrix.pe.py '/scratch/bioinfo2/ryang24/irfan/mapping/'$i.nodup.bed common.bed.summit.bed 2000 40 /home/ryang24/Genome_file/hg19.chromend 'read'
	epd_python hm_matrix.pe.py '/scratch/bioinfo2/ryang24/irfan/mapping/'$i.nodup.bed common.bed.summit.bed 10000 200 /home/ryang24/Genome_file/hg19.chromend 'read'
done
