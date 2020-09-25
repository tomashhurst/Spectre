#' Read TIFF files into R and and create spatial data object
#'
#' @usage read.spatial.files()
#'
#' @param rois vector of ROI names (directory names)
#' @param roi.loc working directory for ROIs
#' @param multi.tiff Default FALSE
#' @param correct.extent Default TRUE
#' @param flip.y Default TRUE
#' @param value.modifier Default 65535
#' @param ext Default = ".tiff"
#'
#' @return Returns a spatial data object.
#'
#' @examples
#'
#' @import data.table
#'
#' @export

read.spatial.files <- function(rois,
                               roi.loc = getwd(),
                               multi.tiff = FALSE,
                               correct.extent = TRUE,
                               flip.y = TRUE,
                               value.modifier = 65535,
                               ext = ".tiff"){

  ### Checks

      if(multi.tiff == TRUE){
        if(length(grep(".tif", rois)) != length(rois)){
          stop("It appears that your list of ROIs are not TIFF stack files, and might be directories full of single TIFFs (i.e. one TIFF per channel. If this is correct, please use 'multi.tiff = FALSE'")
        }
      }

  ### Setup
      message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")
      ROI.list <- list()
      spatial.dat <- list()
      setwd(roi.loc)


  ### Loop for ROIs -- one TIFF per ROI (i.e. tiff stack)

      if(multi.tiff == TRUE){
        message("Multi.tiff is not currently supported")
        #
        #   setwd(roi.loc)
        #   ROI.list <- list()
        #
        #   for(i in rois){
        #     # i <- rois[[1]]
        #     message(paste0("Reading TIFF file for ", i))
        #
        #     temp <- readTIFF(i, all = TRUE)
        #     str(temp)
        #
        #     temp <- raster(temp[[1]])
        #
        #     raster.stack <- stack(active.roi)
        #
        #
        #
        #
        #
        #     raster.stack <- stack(active.roi)
        #     ROI.list[[i]]$rasters <- raster.stack
        #   }
      }

  ### Loop for ROIs -- one FOLDER per ROI

      if(multi.tiff == FALSE){
        for(i in rois){
          # i <- rois[[1]]
          message(paste0("Reading TIFF files for ", i))

          setwd(roi.loc)
          setwd(i)
          tiffs <- list.files(pattern = ext)

          ## TIFF loop
          active.roi <- list()

          for(a in tiffs){
            # a <- tiffs[[1]]
            active.roi[[a]] <- readTIFF(a)
            active.roi[[a]] <- raster(active.roi[[a]])

            if(correct.extent == TRUE){
              extent(active.roi[[a]]) <- c(0, dim(active.roi[[a]])[2], 0,dim(active.roi[[a]])[1]) # Y axis - X axis
            }

            if(flip.y == TRUE){
              active.roi[[a]] <- flip(active.roi[[a]], 'y')
            }

            values(active.roi[[a]]) <- values(active.roi[[a]])*value.modifier

            for(n in c(1:length(names(active.roi)))){
              names(active.roi)[n] <- gsub(ext, "", names(active.roi)[n])
            }
          }

          raster.stack <- stack(active.roi)
          ROI.list[[i]]$RASTERS <- raster.stack
        }
      }


  ### Return
      message("Constructing of spatial data object complete")
      return(ROI.list)

}
