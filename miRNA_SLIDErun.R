rm(list = ls())
cat("\014")

library(EssReg)
library(doParallel)
library(dplyr)
library(pROC)
library(ROCR)
library(SLIDE)
library(matrixStats)
library(ggplot2)
#library(Yamm)
registerDoParallel(detectCores())

#x=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/3way_revisedLabels.processed.csv", row.names=1)
x=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/unsupV2/class2_AbScAg.MFI50_30May23_V2_X.csv", row.names=1)
y=read.csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/unsupV2/class2_AbScAg.MFI50_30May23_V2_Y.csv", row.names = 1)
x_std <- scale(x, T, T)
#dataset=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/X_maxExpression_miRNA.csv", row.names=1)
##applying filters on dataset
# subtrating the min of median across all columns =3.6 from every element of dataframe
#sub_minmedian <- function(x){return (x-3.6)}
#df_filter1=data.frame(lapply(dataset,sub_minmedian))
#x_1=df_filter1
##loading er result
er_result=readRDS("final_delta_0.1_lambda_1.rds")
##getting zs
z <- predZ(x_std, er_result)
colnames(z)=paste0("Z",seq(1,ncol(z)))
#setwd("/ix/djishnu/Trirupa/ABomics.Prj/EssReg-main/ER_kidneyTransplant/miRNA_ABMR/maxExpression_filtered/fdr1.0/del0.2_lam1.0/")

#write.csv(z,file="zmatrix_c2.MFI50.d0.05_lam1.auc.csv",row.names=TRUE)
z=read.csv("zmatrix_c2.MFI50.d0.1_lam1.auc.csv",row.names=1)
#z=read.csv("zmatrix.COV.maxEXPmuvar.miRNA.d0.01_lam1.corr.csv",row.names=1)

## running SLIDE (without interaction terms) ##significant latent factors
SLIDE_marginal_res =SLIDE(z,y,niter=1000,do_interacts = F, spec=0.2, f_size = 39,fdr = 0.1)
SLIDE_marginal_res
#saveRDS(SLIDE_marginal_res, file="hierER.3waySchisto.del0.001_lam1.rds")
#saveRDS(SLIDE_marginal_res, file="hierER.c2.MFI50_SchistoAg.d0.1_lam1_spec0.1.rds")
readRDS(file="hierER.c2.MFI50_SchistoAg.d0.1_lam1_spec0.2.rds")
A_abs <-  apply(er_result$A,c(1,2),abs)
A_k <- as.data.frame(A_abs) %>% select(10)
var_name <- rownames(subset(A_k,A_k > 0.08)) 
var_name

sum(er_result$A[, ] > 0.1)


#### getting marginals at different loading thresholds


summary(lm(as.matrix(y)~z[,21]+z[,75]))

#"z1" "z9" "z23" "z33" "z45” 

### Wilcox test on on features in significant latent factors 
result_df=data.frame(name=NULL,pval=NULL, corrXiY=NULL)

for (i in var_name){res=wilcox.test(x[,i]~as.matrix(y)); result_df<-rbind(result_df,data.frame(name=i,pval=res$p.value, corrXiY=cor(x[,i],y)))}
ii <- order(as.numeric(result_df$pval),decreasing = F)
result_df_sorted<- result_df[ii,]
row.names(result_df) <- result_df$name
mirII  <- intersect(result_df$name,row.names(A_k))
#corr_XiY =cor(x[,i])
final_df <- cbind(result_df[mirII,],A_k[mirII,])
colnames(final_df)<- c("name","pval","corrXiY","A")
final_df_pvalsort=final_df %>% arrange(desc(A))

#write.csv(final_df_pvalsort,file="EvsNE.pval_A_corr.csv")
#z= solve(t(er_result$A)%*%er_result$A)%*%t(er_result$A)%*%t(x)

### finding which Zs have our univariate features in.
z_sigunivariate <- data.frame('miR.3162' = A_abs['miR.3162',], 'miR.7'=A_abs['miR.7',],'miR.363'=A_abs['miR.363',],'miR.133b'=A_abs['miR.133b',],'miR.196a'=A_abs['miR.196a',],'miR.18b'=A_abs['miR.18b',],'miR.505'=A_abs['miR.505',],'miR.3613'=A_abs['miR.3613',] )

dataset=read.csv("/ix/djishnu/Trirupa/ABomics.Prj/ERinputs_kidneyTransplant/miRNA_data/X_maxExpression_miRNA.csv", row.names=1)
##applying filters on dataset
 #subtrating the min of median across all columns =3.6 from every element of dataframe
sub_minmedian <- function(x){return (x-3.6)}
df_filter1=data.frame(lapply(dataset,sub_minmedian))
x_1=df_filter1
col1=x_1[,1]
mir_means=colMeans(x_1)
mir_vars=colVars(as.matrix(x_1[sapply(x_1, is.numeric)]))
mir_variance=x_1 %>% summarise_if(is.numeric, var)

mir_MuVar=data.frame('variance'=rowVars(t(x_1)),'mu'=rowMeans(t(x_1)),row.names = colnames(x_1))
highlight_df=mir_MuVar[c("miR.3162","miR.7",'miR.363','miR.133b','miR.196a','miR.18b','miR.505','miR.3613'),]
mir_MuVar %>%ggplot(aes(x=variance,y=mu))+ geom_point(alpha=0.3) +geom_point(data=highlight_df, aes(x=variance,y=mu),color='red') + geom_point(data=filteredVals, aes(x=variance,y=mu),color='orange')
heatmap(cor(x[,c("miR.3162","miR.7",'miR.363','miR.133b','miR.196a','miR.18b','miR.505','miR.3613')]))
  
muqie(mir_MuVar, dm=c(1,2), probs=0.25, nsegs=1,
      nprojs=2000, reltol=0.001, plot.it=FALSE, 
      full.return=FALSE, xlab=NULL, ylab=NULL)
#meadianMu <- quantile(mir_MuVar$mu,probs = 1.00)
#filteredmus  <- mir_MuVar %>% filter(mu<meadianMu)
medianVar= quantile(mir_MuVar$variance,probs = 0.25)
filteredVals  <- mir_MuVar %>% filter(variance<medianVar)
filtered_check= mir_MuVar[!c(rownames(mir_MuVar)%in%rownames(filteredVals)),]

filtered_check %>%ggplot(aes(x=variance,y=mu))+ geom_point(alpha=0.3)

##### this is the X after filtering out the features with low mu(abundance) and low variance [i.e teh bottom 25% of the distribution]
filt_X= x_1[,!c(rownames(mir_MuVar)%in%rownames(filteredVals))]
write.csv(filt_X,file="X_maxEXP.sub3.6.var25.filtered.miRNA.csv",row.names=TRUE)





##original code
muqie(mir_MuVar, dm=c(1,2), probs=0.5, nsegs=1,
      nprojs=2000, reltol=0.001, plot.it=FALSE, 
      full.return=FALSE, xlab=NULL, ylab=NULL)
meadianMu <- quantile(mir_MuVar$mu,probs = 0.5)
filteredmus  <- mir_MuVar %>% filter(mu<meadianMu)
medianVar= quantile(filteredmus$variance,probs = 0.5)
filteredVals  <- filteredmus %>% filter(variance<medianVar)
filtered_check= mir_MuVar[!c(rownames(mir_MuVar)%in%rownames(filteredVals)),]
