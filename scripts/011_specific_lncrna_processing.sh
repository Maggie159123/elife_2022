## extract lncrnas that expressed in more than 2 cell in per celltype && max expression more than 3
#tail -n +2 Female_Oogenesis_lncRNA_SC_expression_2017.txt | awk '{for(i=2; i<=NF; i++) {a[$1]+=$i;if($i>0){b[$1]++}}}END{for(j in a) {if(b[j]>=2 && a[j]>=3) print j"\t"a[j]"\t"b[j]}}' > tmp/Female_Oogenesis_lncRNA_SC_expression_2017_genes.txt1
for i in *_lncRNA_SC_expression_2017.txt
do
  name=${i%_lncRNA_SC_expression_2017.txt}
  tail -n +2 $i | awk '{for(i=2; i<=NF; i++) {a[$1]+=$i;if($i>0){b[$1]++};if($i>c[$1]){c[$1]=$i}}}END{for(j in a) {if(b[j]>=2 && c[j]>=3) print j"\t"a[j]"\t"b[j]"\t"c[j]"\t"a[j]/b[j]}}' OFS='\t' > specific_lncrna2/${name}_lncRNA_SC_expression_2017_genes.txt1
  sed -i "1iID\tSum\tCount\tMax\tMean" specific_lncrna2/${name}_lncRNA_SC_expression_2017_genes.txt1 
done

### other cell type lncrna exclusion
for i in *_lncRNA_SC_expression_2017.txt
do
  name=${i%_lncRNA_SC_expression_2017.txt}
  cat $i | cut -f2- > ./tmp/${name}_expression_2017.txt
done

cat Female_lncRNA_SC_expression_2017_name.txt | cut -f1 > ./tmp/lncRNA_SC_expression_2017_name.txt 
cd tmp

for i in *_expression_2017.txt
do
  name=${i%_expression_2017.txt}
  mv $i ${name}.txt
  paste *_expression_2017.txt > all_except.tmp
  paste lncRNA_SC_expression_2017_name.txt all_except.tmp >${name}_all_except.txt
  tail -n +2 ${name}_all_except.txt | awk '{for(i=2; i<=NF; i++) {a[$1]+=$i;if($i>0){b[$1]++};if($i>c[$1]){c[$1]=$i}}}END{for(j in a) {if(b[j]<2 && c[j]<3) print j"\t"a[j]"\t"b[j]"\t"c[j]}}' >${name}_othertype_no.txt
  mv ${name}.txt ${name}_expression_2017.txt
  rm -f all_except.tmp
done

for i in *_lncRNA_SC_expression_2017_genes.txt1
do
  name=${i%_lncRNA_SC_expression_2017_genes.txt1}
  perl two_sets_gene_compare_same.pl ./tmp/${name}_othertype_no.txt $i > ./tmp/${name}_none_other.txt
done

for i in *_none_other.txt
do
  name=${i%_none_other.txt}
  mv $i ${name}.txt
  cat *_none_other.txt  > all_except.txt
  perl ../two_sets_gene_compare_diff.pl all_except.txt ${name}.txt > ${name}_specific_final.txt
  mv ${name}.txt ${name}_none_other.txt
  rm -f all_except.txt
done
