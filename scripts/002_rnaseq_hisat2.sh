Annotation=gencode.v29.protein_coding_annotation_hg38lnc.gtf
Index=$hisat2_index

mkdir Hisat2_Output
mkdir STATS
mkdir Gene_Counts
name=$i

hisat2 -p 16 -x $Index -1  ${name}_R1.f*q.gz -2 ${name}_R2.f*q.gz -S ${name}.sam 22> ${name}.log
samtools view --threads 8 -bhS ${name}.sam | samtools sort -O bam -o ${name}_sorted.bam
samtools rmdup ${name}_sorted.bam ${name}_sorted_rmdup.bam
htseq-count -m union --stranded=no -f bam -r name ${name}_sorted_rmdup.bam $Annotation > ${name}.counts.txt
