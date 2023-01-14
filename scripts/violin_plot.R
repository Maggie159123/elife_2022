library(ggplot2)
library(ggpubr)

args <- commandArgs(T)
data <- read.table(args[1], header=T, sep='\t')
data$Type <- as.factor(data$Type)

pdf(args[2])
p <- ggplot(data, aes(x=Type, y=Corr, fill=Type)) + geom_violin()+ geom_boxplot(width=0.05, fill="black") + stat_summary(fun.y=median, geom="point", size=2, color="white")+scale_fill_brewer(palette="Dark2")+ theme(legend.position="none") + scale_x_discrete(limits=c("fGC_Meiotic","Random")) + scale_y_continuous(limits=c(-1,1), breaks=seq(-1,1,0.5)) + stat_compare_means(label.y = 0.5, label.x = 1.3, method="t.test")

p + theme(panel.grid = element_blank()) + theme_bw() + theme(axis.line = element_line(colour = "black",size=0.3))

dev.off()
