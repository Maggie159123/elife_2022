mkdir FastQC_Output
fastqc *.f*q.gz -o FastQC_Output 
multiqc FastQC_Output
