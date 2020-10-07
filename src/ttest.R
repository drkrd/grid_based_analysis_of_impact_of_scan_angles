library(stringr)

nms <- ls(pattern = "output_")
name=nms[2]
t_df <- {}
for (name in nms)
{
  rst_stack <- stack(get(name), allp95)
  rst_df <- as.data.frame(rst_stack)
  rst_df_7m <- rst_df[rst_df$zq95>7,]
  
  nm <- str_remove(name,"output_")
  for(i in 1:4)
  {
    txt <- paste0("cl1","cl",toString(i+1),sep="")
    v2 <- paste0("rst_df_7m$cl",toString(i+1),sep="")
    df_ss <- rst_df_7m[,c(1,(i+1))]
    df_ss <- df_ss[complete.cases(df_ss),]
    ttest <- t.test(df_ss[,1], df_ss[,2], paired = TRUE, alternative = "two.sided", conf.level = 0.95)
    t_df <- rbind(t_df,c(nm,txt,ttest$p.value, ttest$estimate[[1]], ttest$conf.int[1], ttest$conf.int[2]))
  }
  t_df <- as.data.frame(t_df)
  
}

colnames(t_df) <- c("Metric", "Diff", "P-Value", "Mean_diff", "CI1", "CI2")

diff$cl1cl2$p.value





