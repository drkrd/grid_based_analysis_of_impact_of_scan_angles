library(lidR)
library(dplyr)
library(future)
library(e1071)

start_time <- Sys.time()

lasc <- readLAScatalog("D:/1_Work/4_Bauges/LAS_Hauteur_500m/las/")
las_check(lasc)

rs <- 30 #resolution of the grid
area_filt <- 0.9*rs*rs #min area threshold for the area, in a pixel, sampled by a flight line 


##defining processing parameters
plan(multisession, workers = 4L)
opt_chunk_buffer(lasc) <- 10
opt_chunk_size(lasc)   <- 720
opt_stop_early(lasc) <- TRUE
#opt_filter(lasc) <- "-drop_withheld" #note sure about this. Buffer issue is still causing confusion
#plot(lasc, chunk_pattern = TRUE)
opt_select(lasc) <- "xyzClassificationScanAngleRankgpstimeReturnNumber"
opt <- list(raster_alignment = rs, automerge = TRUE)


##creating tables containing the classes, medians etc
class_list <- list(0:10, 10:20, 20:30, 30:40, 40:50)
class_median <- list(0, 15, 25, 35, 45)
names(class_list) = letters[seq_along(class_list)]
breaks = c(min(class_list[[1]]), sapply(class_list, max))
med_data = data.frame(median = unlist(class_median),
                      class = names(class_list))
med_data <- mutate(med_data, class = as.vector(class))



##lidr syntax.
##Only by using catalog_apply template was I able to compute the flightline information on the fly
##The retrieve_flightlines function works on las object which is read each time

##grid_metrics didn't seem feasible as any user_defined function receives only the attributes and 
#NOT a las object which could be processed or analysed using functions built for las objects
##I could not figure out how to use these attributes and compute flight line information using the function

##This seemed logical to me as, at the level of a pixel or cell, I am only interested
##in the unique flight lines from which the cells were sampled. Here, this happens at the level of the cluster
##that is read in readLAS

met_calc = function(cluster, res)
{
  las = readLAS(cluster)
  if (is.empty(las))
    return(NULL)
  las <- retrieve_flightlines(las, dt = 30)
  out <-
    grid_metrics(las, mymets(X, Y, Z, ScanAngleRank, flightlineID, ReturnNumber), rs)
  bbox   <- raster::extent(cluster)
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



##All computations happen in the following function
mymets = function(x, y, z, sc, flid, rn)
{
  #x,y,Z co-ordinates
  #sc-scan angle
  #flid-flight line ID
  #rn-return number
  
  
  
  dframe <- as.data.frame(cbind(x, y, z, sc, flid, rn))
  #added the following because some flight lines had less than 3 points. Area computation was returning an error.
  dframe <- dframe %>%
    group_by(flid) %>%
    filter(n() >= 1000) %>% #1000 is an arbitrary number. Generally, a flight line that covers an entire area of 30m x 30m has no. of points far greater than 1000
    ungroup()

  
  
  ##if there is not data left in an area after the previous filtering, 
  #we jump to else condition to return a stack of five NA values
  if (nrow(dframe) > 0)
  {
    flist <- unique(dframe$flid) #get the unique flight lines
    meanlist <- c()
    fl <- data.frame()
    ##for all the unique flight lines the following are computed:
    #the means of the scan angles of the points in a flight line are computed
    #the area of the pixel covered by the flight line
    for (i in flist)
    {
      dframe1 <- dframe[which(dframe$flid == i),]
      ch_ar <- area_calc(dframe1)
      mean_val <- (abs(mean(dframe1$sc)))
      fl <- rbind(fl, c(i, mean_val, ch_ar))
    }
    fl <- as.data.frame(fl)
    
    names(fl) <- c("flist", "meanlist", "area")
    ##drop all flight lines that do not satisfy the threshold criteria i.e. 90% of the area of a cell
    fl <- fl %>% filter(area>area_filt)
    
    print(fl)
  
    ##if there are no flight lines in a cell from which the cell was sampled insufficiently, we skip to else condition to return
    #a stack of NA values
    
    ##from among the flight lines that were retained, we now pick one flight each for the classes. 
    
    ##if there are five classes, we will pick five flight lines, one per class.
    ##if there are multiple flight lines per class, the flight lines whose mean scan angle, for the cell, is closest to class median 
    ##if there are no flight lines in a class, we will return NA value for that class

    if (nrow(fl) > 0)
    {
      ##The light lines are classified based on the mean scan angle
      fin_flist <-   fl %>%
        # assign classes
        mutate(class = cut(meanlist,
                           breaks = breaks,
                           labels = names(class_list)),
               class = as.vector(class)) %>%
        # get medians from initial data
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
      
      
      fin_flist$flist <- as.numeric(fin_flist$flist)
      fin_flist$meanlist <- as.numeric(fin_flist$meanlist)
      fin_flist$area <- as.numeric(fin_flist$area)
      fin_flist$median <- as.numeric(fin_flist$median)
      
      print(fin_flist)
      

      

      
      #Computation of the metrics and compiling in a list. Each value in a list is the computed metric for a given class
      met_list <- vector(mode="numeric", length=5)
      for (i in 1:length(fin_flist$flist))
      {
        if (fin_flist$flist[i] < 90)
        {
          dframe2 <- dframe[which(dframe$flid == fin_flist$flist[i]),]
          met_list[i] <- func_varch(dframe2$z, dframe2$rn)
        }
        else
        {
          met_list[i] <- NA_real_
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
    else
    {
      mets = list(
        cl1 = NA_real_,
        cl2 = NA_real_,
        cl3 = NA_real_,
        cl4 = NA_real_,
        cl5 = NA_real_
      )
      return(mets)
    }
  }
  else
  {
    mets = list(
      cl1 = NA_real_,
      cl2 = NA_real_,
      cl3 = NA_real_,
      cl4 = NA_real_,
      cl5 = NA_real_
    )
    return(mets)
  }
}



test_mean <- catalog_apply(lasc, met_calc, res = rs, .options = opt)

end_time <- Sys.time()

end_time - start_time
?grid_canopy
