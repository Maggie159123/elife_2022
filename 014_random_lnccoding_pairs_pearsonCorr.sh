#!/bin/bash

lnc_file=$1
coding_file=$2
six_types=$3
cell_stage=$4

mkdir random_pairs

perl test.pl Female_Meiotic_lncRNA_SC_expression_2017_genes.txt1 NONCODEv5_human_XH.txt > NONCODEv5_human_XH_Female_Meiotic_lncRNA.txt

perl test1.pl Female_Meiotic_XH_pearson_corr.txt NONCODEv5_human_XH_Female_Meiotic_lncRNA.txt > Female_Meiotic_XH_pearson_corr.bed

cat Female_Meiotic_XH_pearson_corr.bed | cut -f1 | uniq -c | sort -k2,2V  |  awk '{print $1,$2}' OFS='\t' > Female_Meiotic_XH_pearson_corr_num.txt 

bash random_extract.sh
cat random_Female_Meiotic_XH_same_num.txt | awk '{print $3,$5,"Random"}' OFS='\t' > random_Female_Meiotic_XH_same_num_types.txt

sed -i "1iID\tID\tType" random_Female_Meiotic_XH_same_num_types.txt
perl test2.pl random_Female_Meiotic_XH_same_num_types.txt Female_Meiotic_lncRNA_SC_expression_2017.txt > random_Female_Meiotic_XH_SC_expression_2017_108cells.txt

bash extract_lnc_expressed_sc.sh random_Female_Meiotic_XH_SC_expression_2017_108cells.txt Female_Meiotic_scRNA_cells_coding_expression.txt random_Female_Meiotic_XH_same_num_types.txt Female_Meiotic

perl two_sets_gene_compare.pl random_Female_Meiotic_XH_lnc_coding_pearson_correlation.txt random_Female_Meiotic_XH_same_num_types.txt > random_Female_Meiotic_XH_lnc_coding_pearson_correlation.txt1
rm -f random_Female_Meiotic_XH_lnc_coding_pearson_correlation.txt
mv random_Female_Meiotic_XH_lnc_coding_pearson_correlation.txt1	random_Female_Meiotic_XH_lnc_coding_pearson_correlation.txt

cat Female_Meiotic_XH_pearson_corr.txt random_Female_Meiotic_XH_lnc_coding_pearson_correlation.txt > Female_Meiotic_XH_vs_Random_sam_numandchrom_lnc_coding_pearson_correlation.txt 

sed -i "1iID\tID\tType\tCorr" random_Female_Meiotic_XH_lnc_coding_pearson_correlation.txt
Rscript violin_plot.R Female_Meiotic_XH_vs_Random_sam_numandchrom_lnc_coding_pearson_correlation.txt Female_Meiotic_XH_vs_Random_sam_numandchrom_lnc_coding_pearson_correlation.pdf
