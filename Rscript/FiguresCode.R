# ------------- requires -------------------
library(data.table)
library(ggsci)
library(tidyverse)
library(dplyr)
library(Vennerable)
library(DESeq2)
library(ggpubr)
library(randomForest)
library(pROC)
library(caret)
library(forestmodel)
library(survminer)
library(survival)
source("./Rscript/fun_to_boxplot.R")
source("./Rscript/expr2MITAscore.R")
# -------------Candidate MITA genes, n =1219  ------------
ori_MITAgenes<-fread("./data/0.GOresource.csv",data.table = F)
cd1<-intersect(ori_MITAgenes$Immune.response,ori_MITAgenes$Response.to.external.biotic.stimulus)%>%na.omit()
cd2<-intersect(ori_MITAgenes$Defense.response,ori_MITAgenes$Response.to.external.biotic.stimulus)%>%na.omit()
cd3<-intersect(ori_MITAgenes$Defense.response,intersect(ori_MITAgenes$Immune.response,ori_MITAgenes$Response.to.external.biotic.stimulus))%>%na.omit()
cd_MITAgenes<-data.frame(Gene=unique(c(cd1,cd2,cd3)))
vennlist<-list(Response.to.external.biotic.stimulus=ori_MITAgenes$Response.to.external.biotic.stimulus,
               Immune.response=ori_MITAgenes$Immune.response,
               Defense.response=ori_MITAgenes$Defense.response)

myve<-Venn(vennlist)
plot(myve, doWeights = F)
# ----------- MITA score example：ICI_GC_SYSUCC data ----------------------- 
load("./data/1.ICIgc_sysucc_files.Rdata")
dds <- DESeqDataSetFromMatrix(countData = ICIgc_sysucc_files$expr,
                              colData = ICIgc_sysucc_files$meta,
                              design = ~ Group)
dds_filter <- dds[rowSums(counts(dds))>1, ]
dds_out <- DESeq(dds_filter)
res <- results(dds_out)
res_deseq <- res[order(res$padj),]
res_diff_data <- merge(as.data.frame(res),as.data.frame(counts(dds_out,normalize=TRUE)),by="row.names",sort=FALSE)
colnames(res_diff_data)[1]<-"Gene"
exprdata<-res_diff_data[,-c(2:7)]
DEGs<-res_diff_data[,1:7]
DEGs_filter<-subset(DEGs,abs(log2FoldChange)>=1.5|pvalue<=0.05)
colnames(DEGs_filter)[3]<-"Log2FC"
MITAscore<-expr2MITAscore(expr = exprdata,DEexpr = DEGs_filter,signature = cd_MITAgenes)

# ------------------ Fig.S1 -------------------------

##' Fig.S1A and B: MITAscore for response to microbes
data<-fread("./data/2.anti-microbes_MITAscore.csv",data.table = F)
data$Group<-ifelse(data$HeathyGroup=="Healthy","Healthy",
                   ifelse(data$BacGroup=="Bacterial_infection","Bacterial_infection",
                          ifelse(data$VirGroup=="Viral_infection","Viral_infection","Other_infection")))
data$Group<-factor(data$Group,levels = c("Healthy", "Bacterial_infection","Viral_infection", "Other_infection"))
my_comparisons<-list( c("Healthy", "Bacterial_infection"),
                      c("Healthy", "Viral_infection"),
                      c("Healthy", "Other_infection"),
                     c("Viral_infection","Bacterial_infection"),
                      c("Viral_infection", "Other_infection"),
                      c("Bacterial_infection", "Other_infection"))
p<-fun_to_boxplot(data,group = "Group",variable = "MITAscore",comparisons = my_comparisons)
S1A <-p+scale_fill_manual(values = c("#0E3BF0","#AFFFDF","#377F5B","#D72323"))
## -------- S1A ---------
S1A
###' Random forest
set.seed(1000)
trainIndex<-sample(nrow(data),nrow(data)*0.5)
trainData<-data[trainIndex,]
testData<-data[-trainIndex,]
trainData$HeathyGroup = as.factor(trainData$HeathyGroup)
testData$HeathyGroup = as.factor(testData$HeathyGroup)
Hg_randomforest <- randomForest(HeathyGroup~MITAscore,
                                data = trainData,ntree=300)
pre_ran <- predict(Hg_randomforest,newdata=testData)
HG_ran_roc <- roc(testData$HeathyGroup,as.numeric(pre_ran),
                  ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE)

#####Bacterial infection group
set.seed(1234)
table(data$BacGroup)
trainIndex1<-sample(nrow(data),nrow(data)*0.5)
trainData1<-data[trainIndex1,]
testData1<-data[-trainIndex1,]
trainData1$BacGroup = as.factor(trainData1$BacGroup)
testData1$BacGroup = as.factor(testData1$BacGroup)
Bg_randomforest <- randomForest(BacGroup~MITAscore,
                                data = trainData1)
pre_ran <- predict(Bg_randomforest,newdata=testData1)
BG_ran_roc <- multiclass.roc(testData1$BacGroup,as.numeric(pre_ran))

##Viral infection group
set.seed(123)
table(data$VirGroup)
trainIndex2<-sample(nrow(data),nrow(data)*0.5)
trainData2<-data[trainIndex2,]
testData2<-data[-trainIndex2,]
trainData2$VirGroup = as.factor(trainData2$VirGroup)
testData2$VirGroup = as.factor(testData2$VirGroup)
Vg_randomforest <- randomForest(VirGroup~MITAscore,
                                data = trainData2)
pre_ran <- predict(Vg_randomforest,newdata =testData2)
VG_ran_roc <- multiclass.roc(testData2$VirGroup,as.numeric(pre_ran))

##' ROC plot
roc1<-HG_ran_roc
roc2<-VG_ran_roc$rocs[[1]]
roc3<-BG_ran_roc$rocs[[1]]
auc1<-round(auc(roc1),2)
auc2<-round(auc(roc2),2)
auc3<-round(auc(roc3),2)
auc_res<-data.frame("HC vs Other"=auc1,
                    "Vir vs HC vs Other"=auc2,
                    "Bac vs HC vs Other"=auc3)

## -------- S1B ---------
plot(roc1, lwd = 3, col = "red",print.auc=F, colorize = F,
     auc.polygon=T,
     max.auc.polygon=F,auc.polygon.col=alpha("lightcoral",alpha=0.1),
     print.thres=F,main='RF_ROC')
plot(roc2, lwd = 3, add=TRUE, lty=2, col = "darkorange1",
     print.auc=F, 
     auc.polygon=TRUE, 
     max.auc.polygon=F,auc.polygon.col=alpha("darkorange1",alpha=0.1), 
     print.thres=F)
plot(roc3, lwd = 3, add=TRUE, lty=2, col = "blue",
     print.auc=F, 
     auc.polygon=T,
     max.auc.polygon=F,auc.polygon.col=alpha("lightblue1",alpha=0.4), 
     print.thres=F)
legend(0.7,0.25, 
       legend = c(paste0("HC-Other: AUC = ",auc1), 
                  paste0( "Vir-HC-Other: AUC = ",auc2),
                  paste0("Bac-HC-Other: AUC = ",auc3)),
       lty = c(2,2,2), col = c("red","darkorange","blue"))

##' Fig.S1C: MITAscore for response to tumor
data<-fread("./data/3.TCGA_MITAscore_meta.csv",data.table = F)
data<-subset(data,Stage!="NA")
data$MITAscore<-ifelse(data$MITAscore>=median(data$MITAscore),">=median","<median")
data$Age<-ifelse(data$Age>=65,">=65","<65")
data$Stage<-factor(data$Stage,levels = c("Low","High"))
coxfit1<-coxph(Surv(OStime,OS) ~ MITAscore+Stage+Age+Gender,data = data) 
coxfit2<-coxph(Surv(PFItime,PFI) ~ MITAscore+Stage+Age+Gender,data = data) 
S1C<-forest_model(coxfit1,factor_separate_line=T,
                theme = theme_forest())

## -------- S1C ---------
S1C 

S1D<-forest_model(coxfit2,factor_separate_line=T,
                  theme = theme_forest())
## -------- S1D ---------
S1D
##' Fig.S1D：MITAscore for ICI response
score<-MITAscore$MITAscore
meta<-fread("./data/4.SYSUCC_GC_cli.csv",data.table = F)
data<-merge(score,meta,by="sampleID")
surv_cut <- surv_cutpoint(
  data,
  time = "PFStime",
  event = "PFS",
  variables = c("MRscore")
)
summary(surv_cut)
data$Group<-sample(c("High","Low"),nrow(data),replace = T)
data$Group<-ifelse(data$MRscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
fit<-survfit(Surv(PFStime,PFS) ~ Group,
             data = data)
S1E<-ggsurvplot( fit,
               data=data,
               pval.method = T,combine = F,
               ggtheme = theme(axis.title = element_text(size = 8),
                               axis.text = element_text(size = 8),
                               title = element_text(size = 8),
                               text = element_text(size = 8),
                               axis.line = element_line(size=0.2),
                               panel.background = element_rect(fill = "white")),
               tables.theme = theme(axis.title = element_text(size = 6),
                                    axis.text = element_text(size = 6),
                                    title = element_text(size = 6),
                                    text = element_text(size = 6),
                                    panel.background = element_blank()),
               risk.table = TRUE,
               pval = TRUE,
               title = "PFS",
               palette = c("darkorange","darkblue"),
               legend.title="MITAscore",
               risk.table.col = "strata",
               #surv.median.line = "hv",
               risk.table.y.text.col = T,
               risk.table.y.text = T )
## -------- S1E ---------
S1E

# --------------- get effective MITAgenes, n =12 --------------
load("./data/5.MRgenes_DEGs_allDatasets.Rdata")
vennlist<-list(
  Anti_tumor =unique(c(MRgenes$Response2Tumor$TCGA$Gene,
                      MRgenes$Response2Tumor$LUAD$Gene,
                      MRgenes$Response2Tumor$PAAD$Gene)),
  Anti_microbe=unique(c(MRgenes$Response2Microbes$Bacterial_infection$Gene,
                        MRgenes$Response2Microbes$Viral_infection$Gene,
                        MRgenes$Response2Microbes$Other$Gene)),
  ICI=unique(c(MRgenes$Response2ICI$ICIgc_SYSUCC$Gene,
               MRgenes$Response2ICI$ICIgbm_ZhaoJ$Gene,
               MRgenes$Response2ICI$ICImelanoma_GideTN$Gene,
               MRgenes$Response2ICI$ICImelanoma_Beth$Gene,
               MRgenes$Response2ICI$ICImelanoma_MosheCD45$Gene,
               MRgenes$Response2ICI$ICImelanoma_MosheCD8T$Gene)))

myve<-Venn(vennlist)
plot(myve, doWeights = T)

