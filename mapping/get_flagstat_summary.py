import sys,os
os.system('samtools view -bS '+sys.argv[1]+' >'+sys.argv[1]+'.bam')
os.system('samtools sort '+sys.argv[1]+'.bam '+sys.argv[1]+'.sorted')
os.system('samtools index '+sys.argv[1]+'.sorted.bam')
os.system('samtools flagstat '+sys.argv[1]+'.sorted.bam')
