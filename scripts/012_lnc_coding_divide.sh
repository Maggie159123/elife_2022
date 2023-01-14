mkdir tmp
mkdir -p types/Intergenic types/Antisense types/Sense

for i in *_lncRNA_transform.gtf
do
	name=${i%_lncRNA_transform.gtf}
	cat $i | sort -k1,1V -k2,2V -k3,3V > tmp/$i

# intergenic
	bedtools closest -a tmp/$i -b human_gencode.v29.basic.annotation.bed -d | grep -w "protein_coding" | awk '{if($16>5000) print}' | sort -u | tac > tmp/${name}_intergenic.tmp
	perl unique_pairing.pl tmp/${name}_intergenic.tmp > types/Intergenic/${name}_intergenic.txt

# antisense
	bedtools window -a tmp/$i -b human_gencode.v29.basic.annotation_protein_coding.bed -w 5000 | sort -u | sort -k1,1V -k2,2V -k3,3V -k9,9V -k10,10V -k11,11V | tac > tmp/${name}_less_5kb.tmp

	perl unique_pairing.pl tmp/${name}_less_5kb.tmp > tmp/${name}_less_5kb_nearby_genes.tmp
	les tmp/${name}_less_5kb_nearby_genes.tmp | awk '{if($4=="+" && $12=="+") print}' > tmp/${name}_11.txt
	les tmp/${name}_less_5kb_nearby_genes.tmp | awk '{if($4=="+" && $12=="-") print}' > tmp/${name}_12.txt
	les tmp/${name}_less_5kb_nearby_genes.tmp | awk '{if($4=="-" && $12=="-") print}' > tmp/${name}_22.txt
	les tmp/${name}_less_5kb_nearby_genes.tmp | awk '{if($4=="-" && $12=="+") print}' > tmp/${name}_21.txt
# XH
	les tmp/${name}_12.txt | awk '{if($2>=$11) print}' > types/Antisense/${name}_XH1.txt
	les tmp/${name}_21.txt | awk '{if($3<=$10) print}' > types/Antisense/${name}_XH2.txt
# XT
	les tmp/${name}_12.txt | awk '{if($3<=$10) print}' > types/Antisense/${name}_XT1.txt
	les tmp/${name}_21.txt | awk '{if($11<=$2) print}' > types/Antisense/${name}_XT2.txt 
# XI
	les tmp/${name}_12.txt | awk '{if($2>=$10 && $3<=$11) print}' > types/Antisense/${name}_XI1.txt
	les tmp/${name}_21.txt | awk '{if($10>=$2 && $11<=$3) print}' > types/Antisense/${name}_XI2.txt
# XO
	les tmp/${name}_12.txt | awk '{if($2<=$10 && $3>=$11) print'} > types/Antisense/${name}_XO1.txt
	les tmp/${name}_21.txt | awk '{if($10<=$2 && $11>=$3) print}' > types/Antisense/${name}_XO2.txt

# sense
# SD
	les tmp/${name}_11.txt | awk '{if($2>=$11) print}' > types/Sense/${name}_SD1.txt
	les tmp/${name}_22.txt | awk '{if($2>=$11) print}' > types/Sense/${name}_SD2.txt
# SU
	les tmp/${name}_11.txt | awk '{if($10>=$3) print}' > types/Sense/${name}_SU1.txt
	les tmp/${name}_22.txt | awk '{if($10>=$3) print}' > types/Sense/${name}_SU2.txt
done
