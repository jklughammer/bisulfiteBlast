#! /bin/bash

sample_dir=$1  #results_pipeline directory of the auto RefFreeDMA pipeline
subset=$2 #species id
topn=1000
db="nt_conv" #"bacterial_conv" #"16SMicrobial_conv" #specify blast db to use (needs to be created first --> see createDB)
subfolder="bisulfiteBlast_$db"
convert="TRUE"
resMotifs="CGG|CGA" #for RRBS select restriction binding sites (i.e. bases with which valid reads can start)

scripts=$(cd "$(dirname $0)"; pwd)

log_dir=$sample_dir/00_RefFreeDMA_log/bisulfiteBlast
mkdir -p $log_dir


#check all the samples
shopt -s extglob

for sample_path in `ls -d $sample_dir/*/fastq/$subset*.bam`; do
sample=$(basename $sample_path "_trimmed_bwaMeth.bam")
echo $sample

if [ ! -s $sample_dir/*/bisulfiteBlast_nt_conv/${sample}_unmapped-blast.stats ]; then
sbatch --export=ALL --get-user-env=L --job-name=bisulfiteBlast_$sample --ntasks=1 --cpus-per-task=4 --mem-per-cpu=4000 --partition=mediumq --time=1-00:00:00 -o $log_dir/${sample}_bisulfiteBlast.log $scripts/bisulfiteBlast.sh $sample_path $topn $subfolder $db $convert $resMotifs $sample
else
echo "done"
fi
sleep 2s
while [  `sacct --parsable --long|grep '|PENDING|'|wc -l` -gt 10 ]; do echo Too many pending!; sleep 5m; done

done


#read -p "Press [Enter] key to continue after all samples are done (check with squeue)."
#Rscript analyse_blast_unmappable.R $sample_dir $topn
