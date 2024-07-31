MITAscore2Survival<-function(MITAscore,meta){
  survdata<-merge(MITAscore,meta,by="sampleID")
  surv_cut <- surv_cutpoint(
    survdata,
    time = "OStime",
    event = "OS",
    variables = c("MITAscore")
  )
  summary(surv_cut)
  survdata$Group<-sample(c("High","Low"),nrow(MITAscore),replace = T)
  survdata$Group<-ifelse(survdata$MITAscore>=as.numeric(summary(surv_cut)[1]),"High","Low")
  fit<-survfit(Surv(OStime,OS) ~ Group,
               data = survdata)
  p<-ggsurvplot( fit,
                 data=survdata,
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
                                      panel.background = element_blank(),
                                      panel.grid.major = element_line(size=1)),
                 #risk.table = TRUE,
                 pval = TRUE,
                 title = "OS",
                 palette = c("#BA0F20", "#1A128A"),
                 #facet.by = "Efficacy",
                 legend.title="MITAscore",
                 risk.table.col = "strata",
                 surv.median.line = "hv",
                 risk.table.y.text.col = T,
                 risk.table.y.text = FALSE )
  return(list(p,fit))
}
