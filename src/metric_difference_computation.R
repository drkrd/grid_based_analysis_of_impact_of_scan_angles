library(ggplot2)
library(raster)
library(lidr)


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
    diff[[txt]] <- c(df_ss[,1]-df_ss[,2])
  }
  
  
  
  diff <- lapply(diff, `length<-`, max(lengths(diff)))
  
  dd  <-  as.data.frame(matrix(unlist(diff), nrow=length(unlist(diff[1]))))
  colnames(dd) <- names(diff)
  
  rst_bp <- ggplot(stack(dd), aes(x = ind, y = values))+geom_boxplot()+ggtitle(toString(nm))
  
  return(rst_bp)
  
}


bp_p701 <- bplot(output_p70,"p70")

bp_p701
