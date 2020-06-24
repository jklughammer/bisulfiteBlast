## bisulfiteBlast

These scripts allow to use blast with bisulfite converted sequences (i.e. sequences in which most Cs have been converted to Ts).
To this end, a blast db is created in which all Cs have been replaced by Ts. 
The creation of this "converted" blast db is descibed in createDB.

The inteneded application is to confirm the species annotation of samples by identifying the "best hits" in the blast nucleotide database for a number of randomly sample reads from reduced representation bisulfite sequencing (RRBS) data. It is intended to work well with the [RefFreeDMA pipeline](https://github.com/jklughammer/RefFreeDMA). 

Requirements:
1. R/ Rscript with packages data.table and ggplot2
2. [blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
3. [samtools](http://www.htslib.org/download/)
4. BLASTDB environment variable: `export BLASTDB=$BLASTDB:"<you path>/blastDB/nt":"<you path>/blastDB/taxdb"`
 



