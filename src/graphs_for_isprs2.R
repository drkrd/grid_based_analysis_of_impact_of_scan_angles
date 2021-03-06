rst_stack <- stack(met_max, p95)
rst_df <- as.data.frame(rst_stack, xy=TRUE)
rst_df_7m <- rst_df[rst_df$zq95>7,]
rst_df_7m <- rst_df_7m[,c(1:6)]
rst_df_7m <- rst_df_7m[complete.cases(rst_df_7m),]
var.test(rst_df_7m[,3], rst_df_7m[,4], alternative = "two.sided", conf.level = 0.95)
t.test(rst_df_7m[,3], rst_df_7m[,6], paired = TRUE, alternative = "two.sided", conf.level = 0.95)

stk <- stack(nos,p95)
nos_df <- as.data.frame(stk)
nos_df <- nos_df[nos_df$zq95>7,]
nos_df <- nos_df[complete.cases(nos_df),]
nos1 <- nos_df[nos_df$V1==5,]
####################################################################################################


library(ggplot2)
library(reshape2)
library(ggpubr)
library(lidR)
library(car)
library(tidyr)
library(dplyr)
library(conflicted)



func_scatter <- function(met, nm)
{
  rst_stack <- raster::stack(met, p95)
  rst_df <- extract(rst_stack, shp, method='simple', cellnumbers=TRUE)
  rst_df <- as.data.frame(rbind(rst_df[[1]], rst_df[[2]], rst_df[[3]], rst_df[[4]], rst_df[[5]]))
  rst_vals <- rst_df %>% 
    drop_na(zq95) %>% 
    dplyr::select(-c(cl5)) %>% 
    dplyr::filter(zq95>7) %>%
    dplyr::select(-zq95) %>% 
    dplyr::rename("Class_1" = cl1,
                  "Class 2" = cl2,
                  "Class 3" = cl3,
                  "Class 4" = cl4) %>% 
    drop_na()
  
  
  
  rst_vals_re <- melt(rst_vals, id.vars = 1:3)
  model <- lm(data = rst_vals, rst_vals[,4]~Class_1)
  linhyp3 <- linearHypothesis(model, c("(Intercept)=0","Class_1=1"))
  
  scplot <- ggscatter(
    data = rst_vals_re,
    x = "Class_1",
    y = "value",
    xlab = "Class 1",
    ylab = "Class i",
    title = nm,
    legend.title = "Classes",
    add = "reg.line",
    color = "variable",
    palette = "jco",
    size = 1,
    ggtheme = theme_pubr()) +
    font("title", size = 18)+
    font("xlab", size = 14, face = "bold")+
    font("ylab", size = 14, face = "bold")+
    font("xy.text", size = 13, face = "bold")+
    font("legend.title", size = 14, face = "bold")+
    font("legend.text", size = 12)+
  geom_abline()+
    
  
  return(scplot)
}



scplot <- ggscatter(
  data = all_data,
  x = "ba_pred_ols1",
  y = "sum_ba_hec",
  xlab = "Predicted BA (m²/hec)",
  ylab = "Reference BA (m²/hec)",
  title = ,
  legend.title = "Classes",
  add = "reg.line",
  palette = "jco",
  size = 2,
  cor.coef = FALSE,
  ggtheme = theme_pubr())+
  geom_abline(linetype = "dashed")+
  font("title", size = 18)+
  font("xlab", size = 16)+
  font("ylab", size = 16)+
  font("xy.text", size = 13)+
  font("legend.title", size = 14)+
  font("legend.text", size = 12)+
  border()+
  stat_cor(aes(label = ..rr.label..), size=5)+
  geom_text(aes(label=id_placette),hjust=-0.5, vjust=0.5, size=4)

scplot

 
nplot_gf2_rip <- func_scatter(met_gfrac2,"gf_riparian")

nplot_gf2_rip



plots1 <- ggarrange(nplot_gf,
                   nplot_ru,
                   common.legend = TRUE)

plots2 <- ggarrange(nplot_p10,
                    nplot_p30,
                    common.legend = TRUE)

plots3 <- ggarrange(nplot_p50,
                    nplot_p90,
                    common.legend = TRUE)

plots4 <- ggarrange(nplot_max,
                    nplot_mean,
                    common.legend = TRUE)


plots5 <- ggarrange(nplot_cv,
                    nplot_sd,
                    common.legend = TRUE)



ggexport(plot_gf2_rip, 
         filename = "plot_gf2_rip.png",
         res = NA,
         width = 800,
         height = 800,
         pointsize = 8)

getwd()



ggsave(plot = plots,
       filename = "plots.png",
       width = 12,
       height = 8,
       dpi = 600)



for(name in lst){
  writeRaster(get(name), filename = name, format="GTiff")
}










ggplot(data = rst_vals)+
  aes(x=cl1, y=cl2)+
  geom_point(size=1.5)+
  geom_abline()+
  geom_smooth(method = "lm", se = FALSE)+
  scale_color_jco()

