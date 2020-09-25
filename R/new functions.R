



do.calculate.area <- function(dat,
                              region){
  
  ### Setup
  
      # dat <- spatial.dat
      # region <- 'regions'
      
  ### Preparation
  
      roi.names <- names(spatial.dat)
  
  ### Processing
    
      area.list <- list()
      
      for(i in roi.names){
        # i <- 8
        
        poly.names <- spatial.dat[[i]]$MASKS[[region]]$polygons@data
        areas <- area(spatial.dat[[i]]$MASKS[[region]]$polygons)
        
        areas <- as.data.table(t(areas))
        names(areas) <- as.character(poly.names[[1]])
        area.list[[i]] <- areas
      }
      
      area.list
      
      area.res <- rbindlist(area.list, fill = TRUE)
      area.res <- cbind(roi.names, area.res)
      names(area.res)[1] <- "ROI"

  ### Return
      
      return(area.res)
  
}







########## cells per region per sample ##########
do.rmv.na = function(dat) {
  # either of the following for loops
  
  # by name :
  for (j in names(dat))
    set(dat,which(is.na(dat[[j]])),j,0)
  
  # or by number (slightly faster than by name) :
  # for (j in seq_len(ncol(dat)))
  #   set(dat,which(is.na(dat[[j]])),j,0)
}


do.region.analysis <- function(dat,
                               sample.col,
                               pop.col,
                               region.col,
                               measure.col,
                               
                               area.table,
                               
                               func = 'mean'){
  
  ### Test data
  
      # dat <- cell.dat
      # sample.col <- "Patient"
      # pop.col <- "cell_type_annot"
      # region.col <- "regions_annot"
      # measure.col <- cellular.cols
      # 
      # area.table <- area.res
      # 
      # func <- 'mean'
      
  ### Setup
  
      pops <- unique(dat[[pop.col]])
      samples <- unique(dat[[sample.col]])
      regions <- unique(dat[[region.col]])
  
  ### Counts of pops per region in each sample
      
      all.counts <- data.table()
      
      for(i in samples){
        # i <- samples[[1]]
        
        ## Subset sample data
            samp.dat <- dat[dat[[sample.col]] == i,]
        
        ## Preparation
            samp.counts <- as.data.table(regions)
            names(samp.counts)[1] <- "REGION"
            
        ## Loop for each population across regions
        
            reg.res <- as.data.table(pops)
            names(reg.res)[1] <- "POPULATIONS"
            
            for(a in regions){
              # a <- regions[[1]]
              
              reg.dat <- samp.dat[samp.dat[[region.col]] == a,]
              counts <- reg.dat[, .(count = .N), by = pop.col]
              names(counts)[1] <- "POPULATIONS"
              names(counts)[2] <- a
              
              reg.res <- do.add.cols(reg.res, "POPULATIONS", counts, "POPULATIONS")
              do.rmv.na(reg.res)
            }
        
        ## Additional versions -- distribution and composition
            
            # Type A -- for each cell type, where is it located (row proportions) -- DISTRIBUTION
                
                a.res <- data.table()
                
                for(o in c(1:nrow(reg.res))){
                  # o <- 1
                  nme <- reg.res[o,1]
                  
                  rw.ttl <- sum(reg.res[o,-1])
                  res <- reg.res[o,-1]/rw.ttl
                  res <- res*100
                  
                  a.res <- rbind(a.res, cbind(nme, res))
                  
                  rm(o)
                  rm(nme)
                  rm(rw.ttl)
                  rm(res)
                }
                
                do.rmv.na(a.res)
                a.res
            
            # Type B -- for each region, what cells are in it (column proportions) -- COMPOSITION
                
                b.res <- as.data.table(reg.res[,1])
                
                for(o in c(2:length(names(reg.res)))){
                  # o <- 2
                  nme <- names(reg.res)[o]
                  
                  col.ttl <- sum(reg.res[,..o])
                  res <- reg.res[,..o]/col.ttl
                  res <- res*100
                  
                  b.res <- cbind(b.res, res)
                  
                  rm(o)
                  rm(nme)
                  rm(col.ttl)
                  rm(res)
                }
                
                do.rmv.na(b.res)
                b.res
            
        ## Wrap up
            
            # COUNTS
                
                reg.res
                reg.res.long <- melt(setDT(reg.res), id.vars = c("POPULATIONS"), variable.name = "REGION")
                reg.res.long
                
                reg.res.long.new <- data.table()
                reg.res.long.new$measure <- paste0(reg.res.long$REGION, " -- Cells per region - ", reg.res.long$POPULATIONS, " | Cells per region")
                reg.res.long.new$counts <- reg.res.long$value
                reg.res.long.new
                
                reg.res.long.new <- dcast(melt(reg.res.long.new, id.vars = "measure"), variable ~ measure)
                reg.res.long.new$variable <- i
                names(reg.res.long.new)[1] <- sample.col
                
                
            # COUNTS / AREA
            
                reg.res.by.area <- reg.res
                
                for(u in regions){
                  # u <- regions[[1]]
                  
                  ar <- area.table[area.table[[sample.col]] == i,..u]
                  ar <- ar[[1]]
                  
                  reg.res.by.area[[u]] <- reg.res.by.area[[u]] / ar * 10000 # per 100 um^2
                }

                reg.res.area.long <- melt(setDT(reg.res.by.area), id.vars = c("POPULATIONS"), variable.name = "REGION")
                reg.res.area.long
                
                reg.res.area.long.new <- data.table()
                reg.res.area.long.new$measure <- paste0(reg.res.area.long$REGION, " -- Cells per 100 um^2 of region - ", reg.res.area.long$POPULATIONS, " | Cells per 100 um^2 of region")
                reg.res.area.long.new$counts <- reg.res.area.long$value
                reg.res.area.long.new
                
                reg.res.area.long.new <- dcast(melt(reg.res.area.long.new, id.vars = "measure"), variable ~ measure)
                reg.res.area.long.new$variable <- NULL

            # DISTRIBUTION
                
                a.res
                a.res.long <- melt(setDT(a.res), id.vars = c("POPULATIONS"), variable.name = "REGION")
                a.res.long
                
                a.res.long.new <- data.table()
                a.res.long.new$measure <- paste0(a.res.long$REGION, " -- Percentage of ", a.res.long$POPULATIONS, " in region | Percent of cell type")
                a.res.long.new$counts <- a.res.long$value
                a.res.long.new
                
                a.res.long.new <- dcast(melt(a.res.long.new, id.vars = "measure"), variable ~ measure)
                
                a.res.long.new$variable <- NULL
                
                # a.res.long.new$variable <- i
                # names(a.res.long.new)[1] <- sample.col
                
            # COMPOSITION
                
                b.res
                b.res.long <- melt(setDT(b.res), id.vars = c("POPULATIONS"), variable.name = "REGION")
                b.res.long
                
                b.res.long.new <- data.table()
                b.res.long.new$measure <- paste0(b.res.long$REGION, " -- Composition of region - ", b.res.long$POPULATIONS, " | Percent of cells in region")
                b.res.long.new$counts <- b.res.long$value
                b.res.long.new
                
                b.res.long.new <- dcast(melt(b.res.long.new, id.vars = "measure"), variable ~ measure)
                
                b.res.long.new$variable <- NULL
                
                #b.res.long.new$variable <- i
                #names(b.res.long.new)[1] <- sample.col
        
        ## Adjustments to add it to complete dataset
            
            # long.temp <- melt(setDT(reg.res), id.vars = c("POPULATIONS"), variable.name = "REGION")
            # long.temp
            #
            # long <- data.table()
            # long$measure <- paste0(long.temp$POPULATIONS, " in ", long.temp$REGION)
            # long$counts <- long.temp$value
            #
            # long <- dcast(melt(long, id.vars = "measure"), variable ~ measure)
            # long$variable <- i
            # names(long)[1] <- sample.col
            
        ## Add to 'add.counts'
            
            all.counts <- rbind(all.counts, cbind(reg.res.long.new, reg.res.area.long.new, a.res.long.new, b.res.long.new))
            
            # rm(reg.res.long.new)
            # rm(a.res.long.new)
            # rm(b.res.long.new)
            
            message("... Sample ", i, " complete")
      }

  return(all.counts)
  
}



