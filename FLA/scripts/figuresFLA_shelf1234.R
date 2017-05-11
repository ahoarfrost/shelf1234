require(ggplot2)

maxes <- read.csv('FlaRatesMaxes_shelf1234.csv',row.names=1,colClasses=c("character","numeric","numeric",rep("factor",6)))
factors <- read.csv('FlaRatesWithFactors_shelf1234.csv',row.names=1)

substrateColors <- c("#FFFFFF","#01FFFF","#017F01","#FFFE41","#0000FE","#D50002")
#adjust factor labels
maxes$stn = paste('stn',maxes$stn,sep="")
maxes$depthid <- factor(maxes$depthid,levels=c("d1","d2"),labels=c("surface","bottom"))
factors$stn = paste('stn',factors$stn,sep="")
factors$depthid <- factor(factors$depthid,levels=c("d1","d2"),labels=c("surface","bottom"))

png('../figures/FLA_barplot.png',width=1000,height=600)
a <- ggplot(maxes,aes(x=substrate,y=mean.kcrate.nM.hr)) + geom_bar(aes(fill=substrate),color="gray30",position="dodge",stat="identity") + facet_grid(depthid~stn) + geom_errorbar(aes(ymin=mean.kcrate.nM.hr-sd.kcrate.nM.hr,ymax=mean.kcrate.nM.hr+sd.kcrate.nM.hr),color="grey30",width=0.5)
barplot <- a + theme_bw() + coord_cartesian(ylim=c(0,30)) + theme(strip.text=element_text(size=16),axis.text.x = element_text(angle = 70, hjust = 1,size=16),axis.text.y=element_text(size=16),axis.title.y=element_text(size=20,vjust=2),axis.title.x=element_text(size=20),legend.position="none",strip.background=element_blank(),strip.text.x=element_text(face="bold"),strip.text.y=element_text(face="bold")) + scale_fill_manual(values=substrateColors) + labs(y='Kill-corrected hydrolysis rate (nM/hr)')
print(barplot)
dev.off()

png('../figures/FLA_timeplot.png',width=1000,height=600)
a <- ggplot(factors,aes(x=timepoint,y=mean.kcrate.nM.hr,group=substrate)) + geom_point(aes(fill=substrate),color="gray30",pch=21,size=5,alpha=0.6) + geom_line(aes(color=substrate),size=1,alpha=0.6) + facet_grid(depthid~stn) + geom_errorbar(aes(ymin=mean.kcrate.nM.hr-sd.kcrate.nM.hr,ymax=mean.kcrate.nM.hr+sd.kcrate.nM.hr, color=substrate),alpha=0.3,width=0.5)
timeplot <- a + theme_bw() + coord_cartesian(ylim=c(0,30)) + theme(strip.text=element_text(size=16),title=element_text(size=20),axis.text.x = element_text(angle = 70, hjust = 1,size=16),axis.text.y=element_text(size=16),axis.title.y=element_text(size=20,vjust=2),axis.title.x=element_text(size=20),legend.position="none",strip.background=element_blank(),strip.text.x=element_text(face="bold"),strip.text.y=element_text(face="bold")) + scale_color_manual(values=c("black","#01FFFF","#017F01","#FFFE41","#0000FE","#D50002")) + scale_fill_manual(values=c("black","#01FFFF","#017F01","#FFFE41","#0000FE","#D50002")) + labs(y='Kill-corrected hydrolysis rate (nM/hr)')
print(timeplot)
dev.off()

