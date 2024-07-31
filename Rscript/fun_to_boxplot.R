fun_to_boxplot <- function(data, group, variable,comparisons) {
  p <- ggboxplot(data, x=group, y=variable,fill = group, 
                 add = "jitter")+
    stat_compare_means(comparisons = comparisons,label = "...p.signif.."
                       )+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title.x = element_blank())
  return(p)
}
