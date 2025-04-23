rm(list = ls())
cat("\014")

#### this code is to remove miRNA features that have low mu and low sigma ###
library(matrixStats)
library(Yamm)
library(dplyr)
library(ggplot2)
#loading datasets
Data=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/miRNA_features_labels.csv")
x=Data[,-1]
col1=x[,1]
mir_means=colMeans(x)
mir_vars=colVars(as.matrix(x[sapply(x, is.numeric)]))
mir_variance=x %>% summarise_if(is.numeric, var)

mir_MuVar=data.frame('variance'=rowVars(t(x)),'mu'=rowMeans(t(x)),row.names = colnames(x))
highlight_df=mir_MuVar[c("miR.3162","miR.7",'miR.363','miR.133b','miR.196a','miR.18b','miR.505','miR.3613'),]

#heatmap(cor(x[,c("miR.3162","miR.7",'miR.363','miR.133b','miR.196a','miR.18b','miR.505','miR.3613')]))

muqie(mir_MuVar, dm=c(1,2), probs=0.25, nsegs=1,
      nprojs=2000, reltol=0.001, plot.it=FALSE, 
      full.return=FALSE, xlab=NULL, ylab=NULL)
meadianMu <- quantile(mir_MuVar$mu,probs = 0.5)
filteredmus  <- mir_MuVar %>% filter(mu<meadianMu)
medianVar= quantile(filteredmus$variance,probs = 0.5)
filteredVals  <- filteredmus %>% filter(variance<medianVar)
filtered_check= mir_MuVar[!(rownames(mir_MuVar)%in%rownames(filteredVals)),]

mir_MuVar %>%ggplot(aes(x=variance,y=mu))+ geom_point(alpha=0.3) +geom_point(data=highlight_df, aes(x=variance,y=mu),color='red') + geom_point(data=filteredVals, aes(x=variance,y=mu),color='orange')
filtered_check %>%ggplot(aes(x=variance,y=mu))+ geom_point(alpha=0.3) +geom_point(data=highlight_df, aes(x=variance,y=mu),color='red')

##### this is the X after filtering out the features with low mu(abundance) and low variance [i.e teh bottom 25% of the distribution]
filt_X= x[,!(rownames(mir_MuVar)%in%rownames(filteredVals))]
write.csv(filt_X,file="X_maxEXP.muvariance.filt.miRNA.csv",row.names=TRUE)
