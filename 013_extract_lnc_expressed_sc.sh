#!/bin/bash

lnc_file=$1
coding_file=$2
six_types=$3
cell_stage=$4

mkdir tmp
lnc_num=$(cat $lnc_file | wc -l)

head -n1 $1 > tmp/header_${lnc_file}

head -n1 $coding_file > tmp/header_${coding_file}

for lnc in $(seq 2 $lnc_num)
do
	name=$(cat $lnc_file | tail -n+$lnc | head -n1 | cut -f1)
	#echo $name
	cat $lnc_file | tail -n+$lnc | head -n1 | awk '{for(i=2;i<=NF;i++){if($i>0) print i}}' OFS='\t' >> tmp/${name}_lnc_sc.tmp
	sed -i "1i1" tmp/${name}_lnc_sc.tmp
#	sed ":a;N;s/\n/,/g;ta" tmp/${name}_lnc_sc.tmp > tmp/${name}_lnc_sc.txt

        declare -i num1 ##计数器
        num1=0;

	for sc in `cat tmp/${name}_lnc_sc.tmp`
	do
		cat $lnc_file | cut -f$sc > tmp/${name}_lnc_sc_${num1}.tmp	
		num1=$num1+1
	done
	paste tmp/${name}_lnc_sc_*.tmp | grep -E "^ID|^$name" > tmp/${name}_lnc_sc_expression.txt
	
	for loci in `cat tmp/${name}_lnc_sc.tmp`
	do
		cat tmp/header_${lnc_file} | cut -f$loci >> tmp/${name}_lnc_sc.name
	done

for loci2 in `cat tmp/${name}_lnc_sc.name`
	do
		cat tmp/header_${coding_file} | awk '{for(i=1;i<=NF;i++) {if($i==loci2) print i}}' loci2=$loci2 OFS='\t' >> tmp/${name}_coding_sc.tmp
	done
        declare -i num2 ##计数器
        num2=0;

	for sc2 in `cat tmp/${name}_coding_sc.tmp`
	do
		cat ${coding_file} | cut -f$sc2 > tmp/${name}_coding_sc_${num2}.tmp
		num2=$num2+1
	done
	paste tmp/${name}_coding_sc_*.tmp > tmp/${name}_all_coding_sc_expression.txt

        #rm -f tmp/${name}_lnc_sc_*.tmp
	#rm -f tmp/${name}_coding_sc_*.tmp
	#rm -f tmp/${name}_lnc_sc.tmp 
	#rm -f tmp/${name}_coding_sc.tmp
	#rm -f tmp/*.name

	cat $six_types | grep -E "^ID|^${name}" > tmp/${name}_lnc_coding_pair.tmp

	cat tmp/${name}_lnc_coding_pair.tmp | while read lnc coding type
	do
		cat tmp/${name}_all_coding_sc_expression.txt | awk '{if($1==coding) print $_}' coding=$coding OFS='\t' >> tmp/${name}_coding_sc_expression.tmp	
	done

	wait
	cat -n tmp/${name}_coding_sc_expression.tmp | sort -k2,2 -k1,1n | uniq -f1 | sort -k1,1n | cut -f2- > tmp/${name}_coding_sc_expression.txt
	
	#rm -f tmp/${name}_coding_sc_expression.tmp
	#rm -f tmp/${name}_lnc_coding_pair.tmp
	#rm -f tmp/${name}_all_coding_sc_expression.txt

	Rscript /Share2/home/jijk/users/hejing/Project/LncRNA-seq/WN/extraction/cells_expression/FemaleGC_and_both_Somatic/for_plot/pearson_correlation/expressed_lnc_corresponding_sc/scripts/pearson_correlation.R tmp/${name}_lnc_sc_expression.txt tmp/${name}_coding_sc_expression.txt tmp/${name}_lnc_coding_pearson_correlation.txt
	sed -i "1d" tmp/${name}_lnc_coding_pearson_correlation.txt
	
done

#### merge
cat tmp/*_lnc_coding_pearson_correlation.txt > ./${cell_stage}_six_types_lnc_coding_pearson_correlation.txt

