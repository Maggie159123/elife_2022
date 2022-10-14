awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' sample_peaks.bed >sample_homer.bed
findMotifsGenome.pl sample_homer.bed hg19 motifDir -len 8,10,12
