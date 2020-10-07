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
  
  
  rst_bp <- ggplot(stack(dd), aes(x = ind, y = values))+geom_boxplot()+ggtitle(toString(nm))
  
  return(rst_bp)
  
}


bp_p70_6 <- bplot(output_p70_6,"p70_6")

bp_p70_6


ggsave("bp_ent.png", plot = bp_ent, dpi = 600)

for (name in names){
  nm <- paste0(name,".png")
  ggsave(plot = get(name), filename = nm, dpi = 600 )
}