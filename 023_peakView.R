## usage: Rscript  peakView.R  peaks.bed  IP.sorted.bam  input.sorted.bam  10

#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
if(length(args) <2){
	print(" usage: Rscript  peakView.R  peaks.bed  IP.sorted.bam  input.sorted.bam  10")
}

bedFile   		= args[1] ;
IpBamFile 		= args[2] ;
InputBamFile 	= args[3] ;
depthThreshold	= ifelse(is.na(args[4]),10,args[4]);  

peakRegion=read.table(bedFile,sep='\t',comment.char = "#",header = FALSE,stringsAsFactors=F)
## samtools depth -r chr1:1607777-1608177   SRR1042593.sorted.bam 
#record = as.character(peakRegion[1,])
trimws <-function(char=' abc '){
  tmp <- sub('^\\s+','',char)
  char <- sub('\\s+$','',tmp)
  return(char)
}
#trimws('   knkvnd    ')  ##it's a small test for the function trimws I just defined.
apply(peakRegion,1,function(record){
	cmd=paste0('samtools depth -r ',trimws(record[1]),":",trimws(record[2]),"-",trimws(record[3]),"  ",IpBamFile)
	print(cmd)
	depthTable=system(cmd,intern=TRUE)
	if ( is.na(InputBamFile) ) {
		if ( length(depthTable) >0 ){
			depthList= as.numeric(unlist(lapply(strsplit(depthTable,"\t"),function(x)x[3])) )
			if( max(depthList) > depthThreshold){
				print('Draw figures!')
				png(paste0(record[4],'.png'))
				plot(depthList,type = 'l',col='red',xlab = 'position',ylab = 'depth',main=paste(record,collapse=':'))
				dev.off()
			}
		}
	}else{
		cmd2=paste0('samtools depth -r ',trimws(record[1]),":",trimws(record[2]),"-",trimws(record[3]),"  ",InputBamFile)
		print(cmd2)
		depthTable2=system(cmd2,intern=TRUE)
		if ( length(depthTable) >0  && length(depthTable2) >0  ){
			depthList  = as.numeric(unlist(lapply(strsplit(depthTable ,"\t"),function(x)x[3])) )
			if( max(depthList) > depthThreshold){
				print('Draw figures!')
				depthList2 = as.numeric(unlist(lapply(strsplit(depthTable2,"\t"),function(x)x[3])) )
				png(paste0(record[4],'.png'))
				plot(depthList,type = 'l',col='red',xlab = 'position',ylab = 'depth',main=paste(record,collapse=':'))
				lines(depthList2,col='blue')
				legend("topright",c("IP","Input"),col=c('red','blue'), lty=1)
				dev.off()
			}
		}
	}
})


