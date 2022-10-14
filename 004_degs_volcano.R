
# Usage: Plot volcano plot for RNA-seq DEGs

args <- commandArgs(T)
DEGs_file=read.table(args[1],head=T,sep='\t')

## For DESeq2 DEGs result
my_volcano(DEGs_file[,c(8,4)],0.01,0,0) ## padj<0.01, abs(logFC)>=0
#my_volcano(DEGs_file[,c(8,4)],0.01,1,0) ## padj<0.01, abs(logFC)>=1


######################## DO NOT EDIT #######################

my_volcano <- function(DEG,t_p,t_FC=0,t_marker){
  ## example: print(my_volcano(limma_results[,c(5,1)],0.01,0.6,0))
  library(ggplot2)
  DEG=na.omit(DEG)
  colnames(DEG)=c('p','logFC')
  DEG$gene <- rownames(DEG)
  DEG[DEG$p < 1e-10,'p']=1e-10 # change pvalue limit

#  if (t_FC == 0) {
#    t_FC <- with(DEG, mean(abs(logFC)) + 2 * sd(abs(logFC)))
#  }

  DEG$change= as.factor(ifelse(DEG$p<t_p & abs(DEG$logFC) > t_FC, 
                                    ifelse(DEG$logFC > t_FC,"UP", "DOWN"), "NOT"))
  this_tile <- paste0("Cutoff for logFC is ", round(t_FC, 3),
                      "\nCutoff for pvalue is ", round(t_p, 3), 
                      "\nUp regulated genes ", 
                      nrow(DEG[DEG$change == "UP", ]), 
                      "\nDown regulated genes ", 
                      nrow(DEG[DEG$change == "DOWN", ]))

  p <- ggplot(data=DEG, aes(x=logFC, y =-log10(p),color =change)) +
    geom_point() +
    scale_color_manual(values =c('blue',"black","red"))+
    geom_hline(yintercept = -log10(t_p),lty=4,lwd=0.6,alpha=0.8)+
    geom_vline(xintercept = c(t_FC,-t_FC),lty=4,lwd=0.6,alpha=0.8)+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),   
          axis.line = element_line(colour = "black"))+ggtitle(this_tile) +
    labs( x="log2 (fold change)",y="-log10 (p-value)")+
    theme(plot.title = element_text(hjust = 0.5))
  if (t_marker !=0 ) p = p+geom_text(DEG=subset(DEG, abs(logFC) > t_marker), aes(label=gene),col="green",alpha = 0.5)
  return(p)
}


