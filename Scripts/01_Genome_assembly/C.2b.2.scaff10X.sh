export PATH="/usr/local/software/sScaff10X/bin/:$PATH"
scaff_reads -nodes 10 file.dat reads_1.fastq.gz reads_2.fastq.gz
scaff10x -nodes 30 -align bwa -score 20 -matrix 2000 -reads 10 -longread 0 -gap 100 -edge 50000 -link 8 Acromyrmex_octospinosus_10Xgenomics.fasta reads_1.fastq.gz reads_2.fastq.gz Acromyrmex_octospinosus.genome.scaff10x.fasta
echo Ok
