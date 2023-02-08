# Add reads and genome in fasta
READ1=/home/projects/ku_00039/people/dinghe/data/GAGA/Raw_genome_reads/ncbi_sra/fq/NCBI-0006_1.fq.gz
READ2=/home/projects/ku_00039/people/dinghe/data/GAGA/Raw_genome_reads/ncbi_sra/fq/NCBI-0006_2.fq.gz
FASTA=/home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies/NCBI-0006_Cyphomyrmex_costatus_dupsrm_filt.fasta

meryl k=18 count output read1.meryl $READ1
meryl k=18 count output read2.meryl $READ2
meryl union-sum output db.meryl read*.meryl
merqury.sh db.meryl $FASTA output
sh ./merqury-1.3/eval/spectra-cn.sh db.meryl $FASTA output

