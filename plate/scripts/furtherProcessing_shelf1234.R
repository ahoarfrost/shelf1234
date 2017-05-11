#further processing for plate reader rates
#any rates below zero, change to zero; 

######################## bulk ##########################
bulk <- read.csv("PlateRates_shelf1234.csv",header=TRUE,row.names=1)
#subset impt columns
bulk_factors <- data.frame(row.names=rownames(bulk),"avg_potential_rate"=bulk$avg_potential_rate,"avg_sd"=bulk$potential_sd)
#if rate is below 0, change to zero
for (row in 1:nrow(bulk_factors)) { 
    if (bulk_factors[row,"avg_potential_rate"]<0) {
        bulk_factors[row,"avg_potential_rate"] = 0
        bulk_factors[row,"avg_sd"] = 0
    }
}
#add factors
bulk_factors$stn <- gsub(pattern="stn([0-9])-d([0-9])-bulk-([a-zA-Z0-9]+)",replacement="stn\\1",x=rownames(bulk))
bulk_factors$depthid <- gsub(pattern="stn([0-9])-d([0-9])-bulk-([a-zA-Z0-9]+)",replacement="d\\2",x=rownames(bulk))
bulk_factors$substrate <- gsub(pattern="stn([0-9])-d([0-9])-bulk-([a-zA-Z0-9]+)",replacement="\\3",x=rownames(bulk))
bulk_factors$site <- gsub(pattern="stn([0-9])-d([0-9])-bulk-([a-zA-Z0-9]+)",replacement="stn\\1.d\\2",x=rownames(bulk))
#save as factorsbulk
write.csv(bulk_factors, "PlateRatesFinal_shelf1234.csv",row.names=TRUE)


