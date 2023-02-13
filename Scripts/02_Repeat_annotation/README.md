# GAGA repeat annotation
Folder with the scripts used in GAGA to conduct the repeat annotation in the final genome assemblies. 

### Repeat annotation (Note Zijun upload the scripts and add numbers if you think that it is relevant)
We conducted an initial general repeat annotation that we used to generate a soft-masked genome assembly to proceed with gene annotation. 
Homology-based and de novo methods were conducted in combination to identify the transposable elements (TEs). The genome sequences were aligned against the Rebase TE library (v25.03) and TE protein database using RepeatMasker and RepeatProteinMask (version 4.1.2). In addition, RepeatModeler v2.0.2 was used to build de novo L. neglectus repeat library, which was subsequently used to annotate repeats using RepeatMasker. TRF v4.10.0 was used to find tandem repeats with parameters: "Match=2, Mismatch=7, Delta=7, PM=80, PI=10, Minscore=50". Finally, we combined all evidence producing a general repeat annotation, see the results [here](../01_Genome_assembly/GAGA_genome_stats.xlsx).

Run the main script with following commands, then run the output shells step-by step
```
perl Run_repeatAnnotation.pl  -dcpu 50 -dcutf 30 -LTR_FINDER -RepeatModeler -Trf -RepeatMasker -ProteinMask -kcpu 50 -kcutf 50 -queue st.q -pro_code P18Z10200N0102 assembled_genome.fasta
```



### Repeat annotation v2 -- LUKAS
