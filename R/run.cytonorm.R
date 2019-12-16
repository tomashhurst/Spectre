#' run.cytonorm - ...
#'
#' @usage run.cytonorm(x, ...)
#'
#' @param x data.frame. Input sample. No default.
#'
#' This function...
#'
#' @export

run.cytonorm <- function(x,
                         )
{
  ## Demo data

}



# library(devtools)
# install_github('saeyslab/CytoNorm')


library(devtools)
library(CytoNorm)


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# library(BiocManager)
#
# BiocManager::install("FlowSOM") ### FlowSOM

# devtools::install_github(("saeyslab/FlowSOM"))

library(FlowSOM)
library(Spectre)
library(flowCore)
library(umap)

# Identifying the data

dir <- system.file("extdata", package = "CytoNorm")
files <- list.files(dir, pattern = "fcs$")
data <- data.frame(File = files,
                   Path = file.path(dir, files),
                   Type = stringr::str_match(files, "_([12]).fcs")[, 2],
                   Batch = stringr::str_match(files, "PTLG[0-9]*")[, 1],
                   stringsAsFactors = FALSE)
data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]

train_data <- dplyr::filter(data, Type == "Train")
validation_data <- dplyr::filter(data, Type == "Validation")

ff <- flowCore::read.FCS(data$Path[1])
channels <- flowCore::colnames(ff)[c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47,
                                     40, 44, 33, 17, 11, 18, 51, 14, 23, 32, 10,
                                     49, 27, 24, 31, 42, 37, 39, 34, 41, 26, 30,
                                     28, 29, 25, 35)]
transformList <- flowCore::transformList(channels,
                                         cytofTransform)
transformList.reverse <- flowCore::transformList(channels,
                                                 cytofTransform.reverse)



## TA - UMAP of all reference batches

  # run umap (subsampled?)

  ff

  df <- exprs(ff)
  df <- as.data.frame(df)
  df
  dim(df)

  Spectre::run.umap(x = df,
                    use.cols = channels,
                    umap.seed = 42)

  df <- cbind(df, umap.res)

  names(df)
  colour.plot(d = df,
              x.axis = "UMAP_42_X",
              y.axis = "UMAP_42_Y",
              col.axis = "Nd145Di",
              title = "Nd145Di - CD4")




  # assess
  #


# Testing whether clustering is appropriate

fsom <- prepareFlowSOM(train_data$Path,
                       channels,
                       nCells = 6000,
                       FlowSOM.params = list(xdim = 5,
                                             ydim = 5,
                                             nClus = 10,
                                             scale = FALSE),
                       transformList = transformList,
                       seed = 1)


    fsom$metaclustering

### TA - plot clusters on UMP

    # run umap (subsampled ?)
        # visual assessment


    # CV
        # numerical assessment
        cvs <- testCV(fsom,
                      cluster_values = c(5,10))

        cvs$pctgs$`10`

         ?testCV

            #
            # testCV <- function(fsom,
            #                    cluster_values = 3:50,
            #                    plot = TRUE,
            #                    verbose = FALSE) {
            #
            #   ###
            #   cluster_values <- c(5, 10, 15)
            #   plot <- TRUE
            #   verbose <- FALSE
            #
            #
            #
            #   nClus <- fsom$FlowSOM$map$nNodes
            #   cluster_labels <- FlowSOM::GetClusters(fsom)
            #
            #   # Determine metacluster labels
            #   meta_cl <- list()
            #   for(mc in cluster_values){
            #     if(verbose) message("Computing ", mc, " metaclusters")
            #     meta_cl[[as.character(mc)]] <-
            #       FlowSOM::metaClustering_consensus(fsom$FlowSOM$map$codes,
            #                                         mc,
            #                                         seed = 1)
            #   }
            #   meta_cl[[as.character(nClus)]] <- seq_len(nClus)
            #
            #   # Percentages assigned to each of the clusters per file
            #   pctgs <- list()
            #   for(mc in as.character(c(cluster_values, nClus))){
            #     counts <- matrix(0,
            #                      nrow = length(unique(fsom$FlowSOM$data[,"File"])),
            #                      ncol = as.numeric(mc),
            #                      dimnames = list(unique(fsom$FlowSOM$data[,"File"]),
            #                                      as.character(seq_len(as.numeric(mc)))))
            #     tmp <- table(fsom$FlowSOM$data[,"File"],
            #                  meta_cl[[mc]][cluster_labels])
            #     counts[rownames(tmp), colnames(tmp)] <- tmp
            #     pctgs[[mc]] <- t(apply(counts, 1,
            #                            function(x){ 100 * x/sum(x) }))
            #   }
            #
            #   # Coefficient of variation for each of the percentages
            #   cvs <- list()
            #   for(mc in as.character(c(cluster_values, nClus))){
            #     cvs[[mc]] <- apply(pctgs[[mc]],
            #                        2,
            #                        function(x){ stats::sd(x) / mean(x)})
            #   }
            #   res <- named.list(pctgs, cvs, meta_cl)
            #   if(plot){
            #     PlotOverviewCV(fsom,
            #                    res)
            #   }
            #
            #   return(res)
            # }

# Trainin the model

model <- CytoNorm.train(files = train_data$Path,
                        labels = train_data$Batch,
                        channels = channels,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = 6000,
                                              xdim = 5,
                                              ydim = 5,
                                              nClus = 10,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)

# Normalising the data

CytoNorm.normalize(model = model,
                   files = validation_data$Path,
                   labels = validation_data$Batch,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "Normalized",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE)


