###### EDIT INPUT ######
path="/path/to/Gene-annotation-pipeline"
table="/path/to/gene_families.xlsx"
gene_fam_db="/path/to/Gene_families_db"
genome_name="Genome-name"
output_directory="${path}/Data/Genomes/Genome-name"
proteome="/path/to/Proteome-in-fasta-file"
gff="/path/to/gff3-file"
genome="/path/to/Genome-in-fasta-file"
interpro="/path/to/predicted-domains-from-InterPro-in-TSV-format"
run_bitacora="${path}/bitacora-master/runBITACORA_command_line.sh" # Make sure to edit the script from bitacora and add the path of bitacora scripts and gemoma in there.
num_threads=4
############################

python3 Scripts/run_analysis.py $path $table $gene_fam_db $proteome $interpro $gff $genome $run_bitacora $genome_name $output_directory $num_threads
