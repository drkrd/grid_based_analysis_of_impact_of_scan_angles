splot <- function(rst,nm)
{
  rst_stack <- stack(rst, p95)
  rst_df <- as.data.frame(rst_stack, xy=TRUE)
  rst_df_7m <- rst_df[rst_df$zq95>7,]
  rst_df_7m <- rst_df_7m[,c(1:6)]
  rst_vals <- rst_df_7m[complete.cases(rst_df_7m),]
  
  #model1 <- lm(cl2~cl1, data = rst_vals)
  #model2 <- lm(cl3~cl1, data = rst_vals)
  #model3 <- lm(cl4~cl1, data = rst_vals)
  
  rst_vals_2 <- melt(rst_vals, id.vars = c("x","y","cl1"))
  
  plt <- ggplot(data = rst_vals_2) +
    aes(x = cl1,
        y = value,
        colour = variable,
        shape = variable) +
    geom_point(size = 1.6) +
    geom_abline() +
    geom_smooth(method = "lm") +
    xlab(label = "Class 1") +
    ylab(label = "Class i") +
    ggtitle(nm) +
    scale_color_brewer(palette="Dark2")#+
    #theme(legend.title = "Class i")
      
  return(plt)
  
}

plot_mean <- splot(met_mean, "max")
plot_mean






for (i in rasters)
{
  nm = paste0("plot_",i)
  assign(nm, splot(get(i), sub('.*_','',i)))
  
}



str_match("met_max", "met_(.*?)")[2]

sub('.*plot_met_','',"plot_met_max")


plots <- ls(pattern = "plot_met_")

for (plot in plots)
{
  plt <- get(plot)

  nm <- paste0("sp_",sub('.*plot_met_','',plot), ".png")

  ggsave(plot = plt,
         filename = nm,
         dpi = 600)
  
}

