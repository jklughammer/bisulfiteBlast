#! /bin/bash

in_bam=$1
topn=$2 #number of random reads to check
subfolder=$3
db=$4 #nt, nt_conv, nt_amb #blast db to use. Needs to be created first --> see createDB
convert=$5 #TRUE #FALSE #should input sequences be converted? In case of bisulfite converted data coose TRUE.
resMotifs=$6 #"CGG|CGA" # restriction sides in RRBS data.
sample=$7

out=$(dirname $(dirname $in_bam))/$subfolder/$sample
temp="$(dirname $(dirname $in_bam))/TEMP"
mkdir -p $(dirname $(dirname $in_bam))/$subfolder/
mkdir -p $temp

convMotifs=`echo $resMotifs| awk '{orig=$1;conv=$1;gsub(/C/,"T",conv);gsub("[|]","|\\\t",conv);gsub("[|]","|\\\t",orig);print "\\\t"conv"|\\\t"orig}'`

if [ $convert == "TRUE" ]; then
echo "converting"
samtools view -f 4 $in_bam|grep -P $convMotifs|shuf -n $topn|awk '{gsub(/C/,"T",$10);print ">"$1"\n"$10}'> ${out}_unmapped-conv.txt
else
echo "not converting using sequences as they are"
samtools view -f 4 $in_bam|grep -P $convMotifs|shuf -n $topn|awk '{print ">"$1"\n"$10}' > ${out}_unmapped-conv.txt
fi

echo "blast"
blastn -db $db -query ${out}_unmapped-conv.txt -max_target_seqs 100 -num_threads 4 -word_size 15 -evalue 0.00000001 -outfmt "6 qseqid sseqid sscinames scomnames qlen slen sstart send  pident length evalue bitscore qseq"|sort -t $'\t' -k1,1 -k12,12nr >  ${out}_unmapped-blast_raw.txt

echo "filter"
cat ${out}_unmapped-blast_raw.txt|sort -T $temp -u -k1,1 --merge >  ${out}_unmapped-blast.txt

echo "summarize"
cut -f3 ${out}_unmapped-blast.txt|sort -T $temp|uniq -c|sort -T $temp -k1,1nr> ${out}_unmapped-blast.stats


rm ${out}_unmapped-conv.txt ${out}_unmapped-blast_raw.txt


#output alignments (takes quite a long time)
#echo "output alignments"
#blastn -db $db -query ${out}_unmapped-conv.txt -num_alignments 5 -num_threads 4 -word_size 15 -evalue 10 -outfmt "0" >  ${out}_unmapped-blast_alignment.txt

