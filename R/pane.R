#' pane
#'
#' pane() function takes a PCA as input and analyses it via NNLS to describe admixed individuals as a mixture of sources groups.
#' @param pca_input R data.frame of PCA with N PCs.
#' @param as_file R data.frame with two columns: POP and A/S, POP column lists all populations to be considered, A/S indicatea whether the population should be considered as Admixed ('A') or as Source ('S').
#' @param sources R vector indicating the groups that should be considered as Sources.
#' @param admixed R vector indicating the groups that should be considered as Admixed.
#' @param pc_met String indicating the summary function to be used to summarize the coordinates within groups (can be "mean","median","sum","var","length") or any summary function that can be read directly by R function aggregate(). Default is 'mean'.
#' @param pc_weights R data.frame with one column and no header, listing values to be used as weights. Number of values should be identical to the number of the Principal Components considered by PANE.
#' @examples
#' pca = read_eigen('data/TOY.pca.evec') #OR
#' pca = read_flash('data/TOY_flash.pca')
#' example_AS = read.table('data/Example_AS', header = T)
#' pane(pca_input = pca, as_file = example_AS) #OR
#' pane(pca_input = pca, sources = c('Source1','Source2','Source3'), admixed = c('Admixed1','Admixed2')
#' @return Returns a list containing the ancestries proportions per each Admixed group
#' @export

pane <- function(pca_input, as_file, sources = NULL, admixed = NULL, pc_met = "mean", pc_weights = NULL) {

  if(missing(pca_input)) {
    stop('ERROR: pane() needs a PCA matrix as input')
  }

  else if(is.null(sources)&is.null(admixed)&missing(as_file)) {

    stop('ERROR: pane() needs a reference file to assign sources and admixed samples. Assign a file to as_file, OR a vector to source_list AND vector to admixed_list arguments')

  } else if(is.null(sources)&is.null(admixed)&!missing(as_file)) {

    # Extracting only columns with POP and PCs information
    subset_pca = select(pca_input, contains("P"))
    if(ncol(subset_pca)<2) {stop("Can't find any columns with PC1..PCX labels")}

    # Grouping samples by POP name and apply summary statistic given by pc_met to their PCs values
    # Method approach

    if (!is.null(pc_weights)) {
      print('Principal Components will now be weighted based on the values provided')
      if (nrow(pc_weights) != ncol(subset_pca)-1) {stop('ERROR: mismatch between number of variance values and principal components') }

      pca_weighted = subset_pca[,2:ncol(subset_pca)]*t(pc_weights)
      subset_pca = cbind(subset_pca[,1],pca_weighted)
      names(subset_pca)[1] = "POP"
    }

    pca_aggregated = aggregate(subset_pca[,-c(1)], by= list(subset_pca$POP), FUN = pc_met)

    names(pca_aggregated)[names(pca_aggregated) == 'Group.1'] <- 'POP'


    ## Selecting Sources and Admixed based on as_file file

    Admixed_pop = as.character(as_file$POP[as_file$A.S == 'A'])
    Source_pop = as.character(as_file$POP[as_file$A.S == 'S'])

    ## Check on AS_file labels

    as_file$mismatch <- as_file$A.S %in% c("A", "S")
    false_val = as_file[as_file$mismatch == FALSE,]

    if (nrow(false_val) >= 1) {
      pop_names = false_val[,1]
      print(paste0("No A/S assignation, this group will be skipped: ", pop_names))
    }

    ### Source copying vector

    Source_cv = pca_aggregated[pca_aggregated$POP %in% Source_pop,]
    if (nrow(Source_cv) == 0) {stop("Your source vector is empty, check your input files")}
    if(any(is.na(Source_cv))) {
      print(Source_cv)
      stop("The source vector contains NA, check your input files")}

    ### Admixed copying vector

    Admixed_cv = pca_aggregated[pca_aggregated$POP %in% Admixed_pop,]
    if (nrow(Admixed_cv) == 0) {stop("Your admixed vector is empty, check your input files")}
    if(any(is.na(Admixed_cv))) {
      print(Admixed_cv)
      stop("The admixed vector contains NA, check your input files")}

    ## NNLS analyses
    pane_RES = nnls.mat2(donors = as.matrix(Source_cv[,-c(1)]),recipients = as.matrix(Admixed_cv[,-c(1)]))

    ## Rename cold and rows for nnls.mat
    #colnames(pane_RES) = Source_pop
    #rownames(pane_RES) = Admixed_pop

    ## Rename cols and rows for nnls.mat2
    colnames(pane_RES[[1]]) = Source_cv$POP
    rownames(pane_RES[[1]]) = Admixed_cv$POP

    ## printing/returning output matrix
    return(pane_RES)
  }  else if(missing(as_file)&!is.null(sources)&!is.null(admixed)) {

    # Extracting only columns with POP and PCs information
    subset_pca = select(pca_input, contains("P"))
    if(ncol(subset_pca)<2) {stop("Can't find any columns with PC1..PCX labels")}

    # Grouping samples by POP name and apply summary statistic given by pc_met to their PCs values
    # Method approach

    if (!is.null(pc_weights)) {
      print('Principal Components will now be weighted based on the values provided')
      if (nrow(pc_weights) != ncol(subset_pca)-1) {stop('ERROR: mismatch between number of variance values and principal components') }

      pca_weighted = subset_pca[,2:ncol(subset_pca)]*t(pc_weights)
      subset_pca = cbind(subset_pca[,1],pca_weighted)
      names(subset_pca)[1] = "POP"
    }

    ## Aggregate
    pca_aggregated = aggregate(subset_pca[,-c(1)], by= list(subset_pca$POP), FUN = mean)
    names(pca_aggregated)[names(pca_aggregated) == 'Group.1'] <- 'POP'

    ## Group_by
    #pca_grouped_avg <- subset_pca %>% group_by(POP) %>% summarise_all(funs(mean=mean))
    #pca_grouped_int = as.matrix(sapply(pca_grouped_avg, as.numeric))

    ## Selecting Sources and Admixed based on input lists

    ### Source copying vector

    Source_cv = pca_aggregated[pca_aggregated$POP %in% sources,]
    if (nrow(Source_cv) == 0) {stop("Your source vector is empty, check your input files")}
    if(any(is.na(Source_cv))) {
      print(Source_cv)
      stop("The source vector contains NA, check your input files")}

    ### Admixed copying vector
    Admixed_cv = pca_aggregated[pca_aggregated$POP %in% admixed,]
    if (nrow(Admixed_cv) == 0) {stop("Your admixed vector is empty, check your input files")}
    if(any(is.na(Admixed_cv))) {
      print(Admixed_cv)
      stop("The admixed vector contains NA, check your input files")}

    ## NNLS analyses
    pane_RES = nnls.mat2(donors = as.matrix(Source_cv[,-c(1)]),recipients = as.matrix(Admixed_cv[,-c(1)]))

    ## Rename cols and rows
    colnames(pane_RES[[1]]) = Source_cv$POP
    rownames(pane_RES[[1]]) = Admixed_cv$POP

    ## printing/returning output matrix
    return(pane_RES)
  } else {
    stop('ERROR: optional arguments were not found.')
  }
}


#' write_pane
#'
#' write_pane allows to save pane results in a table-like format.
#' @param pane_input R list returned by pane() function
#' @param output_name string containing the file output name
#' @examples
#' \dontrun{
#' pca = read_eigen(pca_input = 'data/TOY.pca.evec')
#' example_as = read.table('data/Example_AS', header=TRUE)
#' pane_results <- pane(pca_input = pca, as_file = example_as)
#' write_pane(pane_results, output_name = 'my_dir/my_pane_results.txt')
#' }
#' @export

write_pane <- function(pane_input,output_name) {

  df_pane = data.frame(pane_input[[1]])
  df_residuals = data.frame(pane_input[[2]])
  names(df_residuals) = 'Residuals'

  final_df = cbind(df_pane,df_residuals)

  write.table(final_df,file=output_name)
}
