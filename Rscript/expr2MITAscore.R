expr2MITAscore<-function(expr,DEexpr,signature){
  MITAgene_DEGs<-subset(DEexpr,Gene%in%signature$Gene)
  MITAgene_DEGs$Regulation<-sample(c("Up","Down"),nrow(MITAgene_DEGs),replace = T)
  MITAgene_DEGs$Regulation<-ifelse(MITAgene_DEGs$Log2FC>=0,"Up","Down")
  DE_Mgene<-MITAgene_DEGs[,c("Gene","Regulation")]
  Total_MITAscore<-sum(subset(MITAgene_DEGs,Log2FC>=0)$Log2FC)+sum(subset(MITAgene_DEGs,Log2FC<=0)$Log2FC)
  message(paste0("Total Score = ",Total_MITAscore))
  upgene<-subset(DE_Mgene,Regulation=="Up")
  downgene<-subset(DE_Mgene,Regulation=="Down")
  Upexpr<-expr[which(expr$Gene%in%upgene$Gene),]
  Downexpr<-expr[which(expr$Gene%in%downgene$Gene),]
  Upexpr_average<-data.frame(sampleID = colnames(Upexpr)[-1],Upscore=apply(Upexpr[,-1], 2, mean))
  Downexpr_average<-data.frame(sampleID = colnames(Downexpr)[-1],Downscore=apply(Downexpr[,-1], 2, mean))
  DE_average<-merge(Upexpr_average,Downexpr_average,by="sampleID")
  MITAscore<-DE_average%>%mutate(MITAscore=Upscore-Downscore)
  return(list(MITAgene_DEGs=MITAgene_DEGs,MITAscore=MITAscore))
}