#FlaRates_MaxesWithFactors_shelf1234

#take final adjusted rates, extract just kill-corrected mean and sd rate (and round to two decimal points),
#add some factor labels, find timepoint with max rate for each incubation set, 
#save as FlaRatesWithFactors_shelf1234.csv

master.bulk <- read.csv("FlaRatesAdjusted_shelf1234.csv",row.names=1)
#add stn, depth id, site (stn.depthid), depth sampled (m), expt type (bulk), substrate, timepoint factor columns

#define new df with just mean kill-corrected rate and its sd
factors.bulk <- data.frame(mean.kcrate.nM.hr=round(master.bulk$mean.kcrate.nM.hr,2), sd.kcrate.nM.hr=round(master.bulk$sd.kcrate.nM.hr,2),row.names=row.names(master.bulk))

#stn
factors.bulk$stn <- factor(sub(pattern="stn([0-9]+)-d([0-9])-bulk-([a-z]+)-t([0-9])",replacement="\\1",row.names(factors.bulk)))
#depthid
factors.bulk$depthid <- factor(sub(pattern="stn([0-9]+)-d([0-9])-bulk-([a-z]+)-t([0-9])",replacement="d\\2",row.names(factors.bulk)))
#site (stn.depth id)
factors.bulk$site <- factor(sub(pattern="stn([0-9]+)-d([0-9])-bulk-([a-z]+)-t([0-9])",replacement="\\1.d\\2",row.names(factors.bulk)))
#depth sampled (m) (retrieved from cruise notebook)
#expt type
factors.bulk$expt <- factor(sub(pattern="stn([0-9]+)-d([0-9])-bulk-([a-z]+)-t([0-9])",replacement="bulk",row.names(factors.bulk)))
#substrate
factors.bulk$substrate <- factor(sub(pattern="stn([0-9]+)-d([0-9])-bulk-([a-z]+)-t([0-9])",replacement="\\3",row.names(factors.bulk)))
#timepoint
factors.bulk$timepoint <- factor(sub(pattern="stn([0-9]+)-d([0-9])-bulk-([a-z]+)-t([0-9])",replacement="t\\4",row.names(factors.bulk)))

#save factors.bulk as FlaRatesFinalWithFactorsEN556_bulk.csv
write.csv(factors.bulk,"FlaRatesWithFactors_shelf1234.csv",row.names=TRUE)

#take factors from factors.bulk, extract only timepoint with maximum hydrolysis rate for each incubation set; save as FlaRatesMaxes_shelf1234.csv
maxes.bulk <- matrix(ncol=ncol(factors.bulk),dimnames=list(NULL,colnames(factors.bulk)))

#loop through each stn, then depth, then substrate, find max activity timepoint for that substrate, and rbind to maxes.bulk matrix
#for one stn...
for (stn in levels(factors.bulk$stn)) {
    #and one depth...
    for (dep in levels(factors.bulk$depthid)) {
        #and one substrate...
        for (sub in levels(factors.bulk$substrate)) {
            subset <- factors.bulk[factors.bulk$stn==stn&factors.bulk$depthid==dep&factors.bulk$substrate==sub,]
            #...find row with max activity
            m <- max(subset$mean.kcrate.nM.hr)
            maxsub <- subset[subset$mean.kcrate.nM.hr==m,]
            #if >1 row is max (e.g. rates are 0 at every timepoint), use first row
            if(nrow(maxsub)>1) {
                maxsub <- maxsub[1,]
            }
            
            #Add max to FLAmax
            maxes.bulk <- rbind(maxes.bulk,maxsub)
        }
    }
}

#remove row with NA (from creating matrix); (keep only rows where no missing data)
maxes.bulk <- maxes.bulk[complete.cases(maxes.bulk)==TRUE,]
#save as FlaRatesMaxes_shelf1234.csv
write.csv(maxes.bulk,"FlaRatesMaxes_shelf1234.csv",row.names=TRUE)

