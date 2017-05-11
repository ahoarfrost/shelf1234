master <- read.csv("PlateRatesFinal_shelf1234.csv",row.names=1,header=TRUE)

library(RColorBrewer)
library(ggplot2)

substrateColors <- brewer.pal(n=7,name="Set2")
master$depthid <- factor(master$depthid,levels=c("d1","d2"),labels=c("surface","bottom"))

#barplot
png('../figures/Plate_barplot.png',width=1000,height=600)
plot <- ggplot(master,aes(x=substrate,y=avg_potential_rate)) + geom_bar(aes(fill=substrate),color="gray30",stat="identity") + facet_grid(depthid~stn) + geom_errorbar(aes(ymin=avg_potential_rate-avg_sd,ymax=avg_potential_rate+avg_sd),color="grey30",width=0.5) 
plot1 <- plot + coord_cartesian(ylim=c(0,85)) + scale_fill_manual(name=substrate,values=substrateColors) + theme_bw() + theme(strip.text=element_text(size=16),axis.text.x = element_text(angle = 70, hjust = 1,size=16),axis.text.y=element_text(size=16),axis.title.y=element_text(size=20,vjust=2),axis.title.x=element_text(size=20),legend.position="none",strip.background=element_blank(),strip.text.x=element_text(face="bold"),strip.text.y=element_text(face="bold"))  + labs(y="potential rate (nmol/L/hr)")
print(plot1)
dev.off()
