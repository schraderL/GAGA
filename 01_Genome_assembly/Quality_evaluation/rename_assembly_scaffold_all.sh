#!/bin/bash

# Usage: bash rename_assembly_scaffold_all.sh folder_with_assemblies

FILES=$1/*fasta

for f in $FILES
do
	perl ~/scripts/rename_assembly_scaffold.pl $f
done

