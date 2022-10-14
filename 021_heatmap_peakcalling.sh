name=$1

echo Convert bam to bigwig
  bamCoverage -b ${name}_sorted_rmdup.bam -o ${name}_rpkm.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --effectiveGenomeSize 2150570000 \
    --extendReads

 # plot profile and heatmap of TSS  
  computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 -R $RefSeq -S ${name}.bw --skipZeros -o matrix1_${name}_TSS.gz --outFileSortedRegions regions1_${name}_genes.bed
  plotHeatmap --dpi 720 -m matrix1_${name}_TSS.gz -out ${name}_TSS.png
  plotProfile --dpi 720 -m matrix1_${name}_TSS.gz -out ${name}_TSS.pdf --plotFileFormat pdf
  echo Finished!

echo Start Peak calling
  macs2 callpeak -c $Control -t ${name}_sorted_rmdup.bam -m 10 30 -p 1e-5 -f BAM -g mm -n $name 2>${name}.masc2.log
