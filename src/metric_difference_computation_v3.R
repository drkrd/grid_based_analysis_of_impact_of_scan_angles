library(ggplot2)

remove_outliers <- function(x, na.rm = TRUE, ...)
{
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 3 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}



bplot <- function(rst,nm)
{
  lst <- c()
  rst_stack <- stack(rst, output_allp95)
  rst_df <- as.data.frame(rst_stack)
  rst_df_7m <- rst_df[rst_df$zq95>7,]
  diff <- list()
  for(i in 1:4)
  {
    txt <- paste("cl1","cl",toString(i+1),sep="")
    df_ss <- rst_df_7m[,c(1,(i+1))]
    df_ss <- df_ss[complete.cases(df_ss),]
    diffr <- c(df_ss[,1]-df_ss[,2])
    diff[[txt]] <- remove_outliers(diffr)
  }
  
  diff <- lapply(diff, `length<-`, max(lengths(diff)))
  
  dd  <-  as.data.frame(matrix(unlist(diff), nrow=length(unlist(diff[1]))))
  colnames(dd) <- names(diff)
  
  y1 <- boxplot.stats(dd$cl1cl2)$stats[c(1,5)]
  y2 <- boxplot.stats(dd$cl1cl3)$stats[c(1,5)]
  y3 <- boxplot.stats(dd$cl1cl4)$stats[c(1,5)]
  y4 <- boxplot.stats(dd$cl1cl5)$stats[c(1,5)]
  y5 <- boxplot.stats(dd$cl2cl3)$stats[c(1,5)]
  y6 <- boxplot.stats(dd$cl2cl4)$stats[c(1,5)]
  y7 <- boxplot.stats(dd$cl2cl5)$stats[c(1,5)]
  y8 <- boxplot.stats(dd$cl3cl4)$stats[c(1,5)]
  y9 <- boxplot.stats(dd$cl3cl5)$stats[c(1,5)]
  y10<- boxplot.stats(dd$cl4cl5)$stats[c(1,5)]
  
  y_all <- rbind(y1,y2,y3,y4)
  #print(y_all)
  #print(min(y_all[,1]))
  #print(max(y_all[,2]))
  
  rst_bp <- ggplot(stack(dd), aes(x = ind, y = values))+
    geom_boxplot()+ggtitle(toString(nm))+
    coord_cartesian(ylim=c(min(y_all[,1]),max(y_all[,2])))
  
  lst <- list(rst_bp, dd)
  
  return(lst)
  
}


bp_sk <- bplot(output_sk,"sk")

bp_sk[[1]]


ggsave("bp_ent.png", plot = bp_ent, dpi = 600)

for (name in names){
  nm <- paste0(name,".png")
  ggsave(plot = get(name), filename = nm, dpi = 600 )
}