#' @title convert_tibble_to_RSE
#'
#' @description Converts NOMe-seq data in tibble format to Ranged Summarized Experiment format.
#'
#' @details Converts NOMe-seq data in tibble format to Ranged Summarized Experiment format.
#'
#' @param NomeMatrix A tibble where each row corresponds to a sample ROI combination,
#' which should contain the following columns:
#' nFragsAnalyzed = the number fragments (read pairs) that was available
#' names = the names of the ROIs that were analyzed
#' SampleName = the names of the samples that were analzyed
#' GCH_DataMatrix = for every combination of ROI and sample, a matrix where the columns
#' correspond to genomic positions from start to end of the ROI and the rows to fragments.
#' 1 = protected from GpC methylation. 0 = not protected from GpC methylation. For example, the ouput of the fetch-NOMe R package.
#' @param ROIs_gr A GRanges object describing the regions of interest (ROIs).
#' @param annots A data.fame with a column for sample names and a column for group names.
#'
#' @return A Ranged Summarized Experiment with an entry for each ROI. The rowData contains information about each ROI,
#' including a ROIgroup.The assays contain:nFragsAnalyzed, describing the number of fragments that were analyzed for each sample/ROI combination.
#' reads, containg a Gpos object for each sample/ROI combination, with a position for each base in the ROI and two metadata columns: protection,
#' a sparse matrix where TRUE stands for Cs protected from methylation, and methylation, where TRUE stands for methylated Cs.
#'
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays assay<- colData rowRanges
#' @importFrom Matrix Matrix
#' @importFrom tidyr pivot_wider
#' @importFrom GenomicRanges GPos seqnames start end strand pos mcols mcols<-
#' @importFrom BiocGenerics cbind
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
convert_tibble_to_RSE <- function(NomeMatrix,ROIs_gr,annots){

  #fragment count assays
  #assay_types <- c("nFragsFetched", "nFragsNonUnique", "nFragsBisFailed", "nFragsAnalyzed")
  assay_types <- c( "nFragsAnalyzed")

  assay_list <- list()
  for (i in seq_along(assay_types)){
    NomeMatrix_wide <- NomeMatrix[NomeMatrix$names %in% names(ROIs_gr),] %>% pivot_wider(id_cols = "names",names_from = .data$SampleName,values_from = assay_types[i])
    NomeMatrix_wide2 <- as.matrix(NomeMatrix_wide[,-1])
    row.names(NomeMatrix_wide2) <- NomeMatrix_wide$names
    assay_list[[i]] <- NomeMatrix_wide2
  }
  names(assay_list) <- assay_types

  #sort the ROIs_gr to be the same order as the assay rownames
  ROIs_gr <- ROIs_gr[row.names(assay_list[[1]]),]

  #combine to RSE
  se <- SummarizedExperiment(colData = annots,
                             rowRanges = ROIs_gr,
                             assays = assay_list)

  #make the GPos objects
  ROIs <- names(ROIs_gr)
  samples <- annots$samples
  gr_list1 <- list()
  for (r in seq_along(ROIs)){

    #use the ROI GRanges to generate a GPos object with an entry for each base
    gpos1 <- GPos(seqnames=seqnames(ROIs_gr)[r], pos=start(ROIs_gr)[r]:end(ROIs_gr)[r], strand=strand(ROIs_gr)[r], stitch=NA)


    #reverse the positions if the motif is on the - strand
    if(as.character(strand(gpos1))[1]=="-"){
      gpos1 <- rev(gpos1)
    }

    #extract the row of the NomeMatrix for the current ROI
    NomeMatrix1 <- NomeMatrix[NomeMatrix$names==ROIs[r],]

    #loop through the samples and make a list of GCH methylation matrices (rows=positions, columns=fragments)
    GCH_matrix_list <- rep(list(matrix(nrow=length(gpos1))),length(samples))
    for (s in seq_along(samples)){
      if(length(NomeMatrix1$GCH_DataMatrix[NomeMatrix1$SampleName==samples[s]])==0){
        GCH_matrix_list[[s]] <- matrix(nrow=length(gpos1),ncol=1,NA)
      }
      else if (is.null(dim(NomeMatrix1$GCH_DataMatrix[NomeMatrix1$SampleName==samples[s]][[1]])[1])){
        GCH_matrix_list[[s]] <- matrix(nrow=length(gpos1),ncol=1,NomeMatrix1$GCH_DataMatrix[NomeMatrix1$SampleName==samples[s]][[1]])
      }

      else if (dim(NomeMatrix1$GCH_DataMatrix[NomeMatrix1$SampleName==samples[s]][[1]])[1] == 0){
        GCH_matrix_list[[s]] <- matrix(nrow=length(gpos1),ncol=1,NA)
      } else {
        GCH_matrix_list[[s]] <- t(NomeMatrix1$GCH_DataMatrix[NomeMatrix1$SampleName==samples[s]][[1]])
      }
      names(GCH_matrix_list) <- samples

    }
    #add this list as mcols to the GPos object
    mcols(gpos1)  <- GCH_matrix_list[which(sapply(GCH_matrix_list,is.null)==FALSE)]

    #add the Gpos object to the list of GPos objects
    gr_list1[[r]] <- gpos1
  }
  names(gr_list1) <- ROIs


  #add GPos objects
  reads <- Reduce(cbind, lapply(annots$samples, (\(sample)
                                                 lapply(seq_along(gr_list1), function(idx) {
                                                   # Hack...GRanges doesn't support list ops
                                                   roi <- gr_list1[idx]
                                                   roi_name <- names(roi)
                                                   x <- gr_list1[roi_name][[1]][, sample]
                                                   # GRangesList(x)
                                                 }))))


  assay(se, "reads", withDimnames = FALSE) <- reads


  ########### Turn all matrices into sparse matrices
  SE_List <- list()
  for(s in seq_along(samples)){
    SE_List[[s]] <- se[,samples[s]]
    assay(SE_List[[s]], "reads", withDimnames = FALSE)  <- cbind(lapply(assays(SE_List[[s]])[["reads"]],makeSparse))
  }
  NomeSE <- do.call(BiocGenerics::cbind,SE_List)
  return(NomeSE)
}

#function to make matrices sparse:
makeSparse <- function(sereads){

  # x1 <- assays(NomeSE)[["reads"]][[1]]
  x1 <- sereads
  #protection
  x <- mcols(x1)[,1]
  x[is.na(x)] <- 0
  xt <- x==1
  protection <- Matrix(xt,sparse=TRUE)
  #methylation
  x <- mcols(x1)[,1]
  x[is.na(x)] <- 1
  xt <- x==0
  methylation <- Matrix(xt,sparse=TRUE)

  #add as mcols
  mcols(x1) <- NULL
  x1$protection <- protection
  x1$methylation <- methylation
  sereads <- x1
}

