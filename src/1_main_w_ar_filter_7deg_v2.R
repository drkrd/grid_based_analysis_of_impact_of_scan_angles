library(lidR)
library(dplyr)
library(future)
library(e1071)



lasc <-readLAScatalog("C:/1_Work/2_Ciron/Data/ALS/norm/test_areas/tile3/ALS_FRISBEE_Ciron_EPSG2154_453200_6356600_norm.las")
start_time <- Sys.time()
rs <- 50 #resolution of the grid
area_filt <- 2250

plan(multisession, workers = 4L)
opt_chunk_buffer(lasc) <- 0
opt_chunk_size(lasc)   <- 100
opt_stop_early(lasc) <- FALSE
plot(lasc, chunk_pattern = TRUE)
opt_select(lasc)       <- "xyzScanAngleRankgpstimeReturnNumber"
#opt_filter(lasc) <- "-drop_z_below 1"


opt <- list(raster_alignment = rs, automerge = TRUE)



class_list <- list(0:7, 7:14, 14:21, 21:28, 28:35, 35:42, 42:49)
class_median <- list(0, 10.5, 17.5, 24.5, 31.5, 38.5, 45.5 )
names(class_list) = letters[seq_along(class_list)]
breaks = c(min(class_list[[1]]), sapply(class_list, max))
med_data = data.frame(median = unlist(class_median),
                      class = names(class_list))
med_data <- mutate(med_data, class = as.vector(class))




met_calc = function(cluster, res)
{
  las = readLAS(cluster)
  if (is.empty(las))
    return(NULL)
  las <- lasflightline(las, dt = 30)
  out <- grid_metrics(las, mymets(X, Y, Z, ScanAngleRank, flightlineID, ReturnNumber), rs)
  bbox <- raster::extent(cluster)
  out <- raster::crop(out, bbox)
  return(out)
}



area_calc = function(dfr)
{
  #print(length(dfr$x))
  ch_pts <- chull(dfr$x,dfr$y)
  ch_pts <- c(ch_pts, ch_pts[1])
  dfr <- dfr[ch_pts,]
  dfr <- dfr %>% 
    select(1:2) 
  ch_poly <- Polygon(dfr, hole=F)
  return(ch_poly@area)
}




gfrac <- function(z,rn)
{
  dfr <- as.data.frame(cbind(z,rn))
  dfr <- dfr[dfr$rn==1,]
  dfrn <- as.data.frame(filter(dfr, (z<=1.5 & rn==1)))
  num=nrow(dfrn)
  den=nrow(dfr)
  gf <- num/den
  return(gf)
}




mymets = function(x, y, z, sc, flid, rn)
{
  #x,y,Z co-ordinates
  #sc-scan angle
  #flid-flight line ID
  
  
  
  
  dframe <- as.data.frame(cbind(x, y, z, sc, flid, rn))
  #dframe <- dframe[dframe$z>=1,]
  flist <- unique(dframe$flid)
  #print(flist)
  #added the following because some flight lines had less than 3 points. Area computation was returning an error.
  dframe <- dframe %>%
    group_by(flid) %>%
    filter(n() >= 1000) %>% 
    ungroup()
  #print(dim(dframe))
  if (nrow(dframe) > 0)
  {
    flist <- unique(dframe$flid) #get the unique flight lines
    
    meanlist <- c()
    fl <- data.frame()
    
    
    
    #For all the unique flight lines, the means of the scan angles of the points in a flight line are computed
    for (i in flist)
    {
      dframe1 <- dframe[dframe$flid == i,]
      ch_ar <- area_calc(dframe1)
      mean_val <- (abs(mean(dframe1$sc)))
      fl <- rbind(fl, c(i, mean_val, ch_ar))
    }
    
    
    names(fl) <- c("flist", "meanlist", "area")
    fl <- fl %>% filter(area>area_filt)
    print(fl)
    
    
    if (nrow(fl) > 0)
    {
      #Each class will have multiple flight lines. From among the flight lines (in a class), the flight line whose mean angle
      #closest to the median of the class is picked as a representative
      fin_flist <-   fl %>%
        # assign classes
        mutate(class = cut(meanlist,
                           breaks = breaks,
                           labels = names(class_list)),
               class = as.vector(class)) %>%
        # get medians
        left_join(med_data, by = "class") %>%
        # within each class...
        group_by(class) %>%
        # keep the row with the smallest absolute difference to the median
        slice(which.min(abs(meanlist - median))) %>%
        # sort in original order
        arrange(flist)
      fin_flist <- as.data.frame(fin_flist)
      #Get the classes which are not present in fin_flist
      miss <- setdiff(class_median, fin_flist$median)
      
      
      
      
      #To ensure the final table fin_flist has 5 rows with each representing a class, the missing rows are replaced with flight line
      # no. "99" and class "x"
      count = 91
      for (i in miss)
      {
        fin_flist <- rbind(fin_flist, c(count, 99, 0, as.vector("x"), i))
        count = count + 1
      }
      fin_flist <- fin_flist %>% arrange(median)
      print(fin_flist)
      
      
      
      
      #Computation of the metrics and compiling in a list. Each value in a list is the computed metric for a given class
      met_list <- c()
      for (i in 1:length(fin_flist$flist))
      {
        if (fin_flist$flist[i] < 90)
        {
          dframe2 <- dframe[dframe$flid == fin_flist$flist[i],]
          met_list <- c(met_list, mean(dframe2$z, na.rm=TRUE))
        }
        else
        {
          met_list <- c(met_list, NULL)
        }
      }
      mets = list(
        cl1 = met_list[1],
        cl2 = met_list[2],
        cl3 = met_list[3],
        cl4 = met_list[4],
        cl5 = met_list[5]
      )
      return(mets)
    }
  }
}





test3_mean <- catalog_apply(lasc, met_calc, res = rs, .options = opt)

end_time <- Sys.time()

end_time - start_time
