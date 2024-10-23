# execute this file before running the pipeline. It will assess the gene families databases files, and generate files called GENE-FAMILY-NAME_rare_db.fasta if it find strange features in the sequences, such as very long, short sequences. 
# this files can be manually analyzed to decide to exclude or not the rare sequences from the db

python3 Scripts/db_analysis.py
