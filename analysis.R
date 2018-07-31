library(beadarray)
library(dplyr)
library(parallel)
library(rlist)
library(grid)
library(ggplot2)
library(scater)
library(stringr)
library(Biobase)


setwd("~/Documents/IFN_array/")

#Calculate SNR for each array                                                                     DONE
arraylist <- c(7656782095,7656766032,7509352086,7509352070,7509352065,7509352063,7509352061,7509352001)
metrics_vec <- vector(mode = "character", length = 0)
metrics_list <- list(name = "list of metrics")

for(a in arraylist){
    metrics_list <- list.append(metrics_list, read.table(paste("~/Documents/IFN_array/",a,"/Metrics.txt", sep = ""), sep="\t", header =TRUE, as.is = TRUE))
    metrics_vec <- c(metrics_vec, paste("ht12metrics", a, sep = "_"))
}
metrics_list <- list.remove(metrics_list, 'name')
names(metrics_list) <- metrics_vec

metrics_list <- lapply(metrics_list, function(x){
  dplyr::mutate(x, snr = P95Grn/P05Grn)
})

#Produce SNR plot for each array (ggplot2)                                                     DONE
plots <- list()
for (a in 1:c(length(metrics_list))) {
  plots[[a]] <- qplot(1:12,snr, data = metrics_list[[a]], main = paste("snr", a, sep = "_"))
}
layout <- matrix(c(1:8), nrow = 2, byrow = TRUE)
multiplot(plotlist = plots, layout = layout)


#Read bead level data                                                                         DONE
targets=read.table("targets.txt",header=TRUE,sep="\t",as.is=TRUE)

#generate section names from targets                                                          DONE
sn = paste ( targets [, 3], targets [, 4], sep = "_")



#read and annotate data                                                                       DONE
data = readIllumina(sectionNames = sn, useImages=TRUE, sampleSheet ="sampleSheet.csv")        
suggestAnnotation(data,verbose=TRUE)
data <- setAnnotation (data ,"Humanv4")
save(data, file="BLData.RData")
save.image()
unlink(".RData")

###QA                                                                                         DONE
sectionNames(data)
boxplot (data, transFun = logGreenChannelTransform, col = "green", ylab=expression(log[2](intensity)),las=2,outline=FALSE, main = "HT-12 IFN Experiment ")

#plot array false colour images                                                               DONE
pdf("imageplot.pdf",paper="a4r", width=15)
for (n in c(1:96)) {
  print (paste("Array", n, sep=" "))
  print(beadarray::imageplot(data, array=n, transFun = logGreenChannelTransform, horizontal=TRUE ,low="lightgreen", high="darkgreen"))
}
dev.off()

#Outlier plot                                                                                 DONE
pdf("outlierplot", paper = "a4r", width = 15)
for (i in c(1: 96)) {
  print(outlierplot(data , array = i, main = paste ( sectionNames ( data )[i]," outliers ")))
}
dev.off()

#plot combined control plots for individual arrays                                           NOT DONE
pdf("combined_controls.pdf", paper="a4r", width=15)
for(n in c(1:96)){
  combinedControlPlot(data, array=n)
}
dev.off()

combinedControlPlot(data, array=1)

#BASH correction                                                                             DONE
BASH_output <- mclapply(c(1:96), function(x) BASH(data, array = x, useLocs = FALSE))
wts_list <- list()
for (i in 1:96) {
  a <- BASH_output[[i]]$wts
  wts_list[[i]] <- a
}
data_BASH = setWeights(data , wts = wts_list , array = c(1:96), wtName="wts")
save(data_BASH, file="BLData_BASH.RData")
table(getBeadData(data_BASH2, array=1, what="wts"))
head (data_BASH[[8]])


# Plot masked region                                                                         DONE 
pdf("masked_region.pdf", paper="a4r", width=15)
for(i in c(1:96)){
  showArrayMask(data_BASH, array=i, override=FALSE)
}
dev.off()


#HULK Normalize probe intensities                                                            NOT SURE
for(j in c(1:96)){
  data_HULK <- HULK(data_BASH2, array=j, neighbours = NULL, invasions = 20, useLocs = TRUE, weightName = "wts", transFun = logGreenChannelTransform, outlierFun = illuminaOutlierMethod)
}

#QC plots and HTML summary                                                                   PROBLEM
library(illuminaHumanv4.db)
expressionQCPipeline(data, transFun = logGreenChannelTransform, qcDir = "QC", plotType = ".jpeg")
dev.off()

###Summarise                                                                                 DONE
myMean = function (x) mean (x, na.rm = TRUE )
mySd = function (x) sd(x, na.rm = TRUE )
greenchannel = new("illuminaChannel", logGreenChannelTransform ,illuminaOutlierMethod , myMean , mySd , "G")
data_sum <- beadarray::summarize(data_BASH, channelList = list(greenchannel))


#save datasumm object
save(data_sum, file="data_sum.RData")

###QC                                                                                        DONE
library (arrayQualityMetrics)
arrayQualityMetrics (data_sum, force = TRUE)


# normalization                                                                              DONE
detectionP <- calculateDetection(data_sum, status=fData(data_sum)$Status, negativeLabel="negative")
Detection(data_sum) <- detectionP
data_sum_log2 <- channel(data_sum, "G")
data_sum_log2 <- addFeatureData(data_sum_log2, toAdd = c("SYMBOL", "PROBEQUALITY", "CODINGZONE", "PROBESEQUENCE", "GENOMICLOCATION"))

data_sum_norm <- normaliseIllumina(data_sum_log2, method="quantile", transform="none" )
data_sum_norm_neqc <- normaliseIllumina(data_sum_log2, method="neqc", transform="none" )

save(data_sum_norm, file="data_sum_norm.RData")
save(data_sum_norm_neqc, file="data_sum_norm_neqc.RData")

pdf("normalization2.pdf", paper="a4r", width=15)
par(mai=c(1.5,1,0.2,0.1), mfrow=c(1,2))
boxplot(exprs(data_sum), ylab=expression(log[2](intensity)), las=2, outline=FALSE)
boxplot(nObservations(data_sum), ylab="number of beads", las=2, outline=FALSE)

par(mai=c(1.5,1,0.2,0.1), mfrow=c(1,2))
boxplot(exprs(data_sum_norm), ylab=expression(log[2](intensity)), las=2, outline=FALSE)
boxplot(nObservations(data_sum_norm), ylab="number of beads", las=2, outline=FALSE)

par(mai=c(1.5,1,0.2,0.1), mfrow=c(1,2))
boxplot(exprs(data_sum_norm_neqc), ylab=expression(log[2](intensity)), las=2, outline=FALSE)
boxplot(nObservations(data_sum_norm_neqc), ylab="number of beads", las=2, outline=FALSE)
dev.off()

#####DE analysis###                                                                          DONE
library (limma)

#Design matrix                                                                               DONE
pd <-pData(data_sum_norm)
sample_name <- pd$Sample_Name
sample_group <- toupper(str_extract(pd$Sample_Name,  "^[A-z][0-9]{1,2}"))


design1 <- model.matrix(~0 + sample_group)
colnames(design1) <- str_remove(colnames(design1), "sample_group")

#Fit linear model                                                                            DONE
fit<-limma::lmFit(exprs(data_sum_norm),design1)

#Define contrast matrix (Timepoints vs 0)                                                    DONE
contrasts_group_1 <- as.vector(read.delim("contrasts1.txt",header = FALSE)[[1]])  
cont.matrix_1 <- makeContrasts(contrasts = contrasts_group_1, levels = design1)
 
#Define contrast matrix (Same time diff treatment)                                   
contrasts_group_2 <- as.vector(read.delim("contrasts2.txt",header = FALSE)[[1]])  
cont.matrix_2 <- makeContrasts(contrasts = contrasts_group_2, levels = design1)

#Define contrast matrix (Between time points)                                   
contrasts_group_3 <- as.vector(read.delim("contrasts3.txt",header = FALSE)[[1]])  
cont.matrix_3 <- makeContrasts(contrasts = contrasts_group_3, levels = design1)

#Extract linear model fit for the contrasts                                                  DONE
fit_1  <- contrasts.fit(fit, cont.matrix_1) 
fit_1  <- eBayes(fit_1)
fit_2  <- contrasts.fit(fit, cont.matrix_2) 
fit_2  <- eBayes(fit_2)
fit_3  <- contrasts.fit(fit, cont.matrix_3) 
fit_3  <- eBayes(fit_3)

#Attach annotation information                                                               DONE
anno <- fData(data_sum_norm)
colnames(anno)
anno <- anno[,c("SYMBOL", "IlluminaID")]
fit_1$genes <- anno
topTable(fit_1)

fit_2$genes <- anno
topTable(fit_2)

fit_3$genes <- anno
topTable(fit_3)

#Decide total DEGs                                                                           DONE
decideTests(fit_1)
table(decideTests(fit_1))
sum(abs(decideTests(fit_1))==1, na.rm = TRUE)

decideTests(fit_2)
table(decideTests(fit_2))
sum(abs(decideTests(fit_2))==1, na.rm = TRUE)

decideTests(fit_3)
table(decideTests(fit_3))
sum(abs(decideTests(fit_3))==1, na.rm = TRUE)

#Visualize the DE analysis results                                                           DONE
volcanoplot(fit_1, highlight = 10, names = fit_1$genes$SYMBOL)
volcanoplot(fit_2, highlight = 10, names = fit_2$genes$SYMBOL)
volcanoplot(fit_3, highlight = 10, names = fit_3$genes$SYMBOL)

#save DE contrast file
write.fit(fit_1 , file = "Result_fit_1.csv", adjust="BH")
write.fit(fit_2 , file = "Result_fit_2.csv", adjust="BH")
write.fit(fit_3 , file = "Result_fit_3.csv", adjust="BH")

save.image(file = "analysis.RData")
unlink(".RData")