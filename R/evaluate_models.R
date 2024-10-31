#' evaluate_models
#'
#' evaluate_models function takes a PCA as input and analyses it via NNLS to describe admixed individuals as a mixture of sources groups.
#' @param pca_input R data.frame of PCA with N PCs
#' @param as_file R data.frame with two columns: POP and A/S, POP column lists all populations to be considered, A/S indicatea whether the population should be considered as Admixed ('A') or as Source ('S')
#' @param sources R vector indicating the groups that should be considered as Sources
#' @param admixed R vector indicating the groups that should be considered as Admixed
#' @examples
#' pca = read_eigen('data/TOY.pca.evec') #OR
#' pca = read_flash('data/TOY_flash.pca')
#' example_AS = read.table('data/Example_AS', header = T)
#' pane(pca_input = pca, as_file = example_AS) #OR
#' pane(pca_input = pca, sources = c('Source1','Source2','Source3'), admixed = c('Admixed1','Admixed2')
#' @return Returns a list containing the ancestries proportions per each Admixed group
#' @export

evaluate_models <- function(pca_input, as_file, sources = NULL, admixed = NULL) {

  if(missing(pca_input)) {
    stop('ERROR: negative_assign() needs a PCA matrix as input')
  }

  else if(is.null(sources)&is.null(admixed)&missing(as_file)) {

    stop('ERROR: negative_assign() needs a reference file to assign sources and admixed samples. Assign a file to as_file, OR a vector to source_list AND vector to admixed_list arguments')

  } else if(is.null(sources)&is.null(admixed)&!missing(as_file)) {

    # Extracting only columns with POP and PCs information
    subset_pca = select(pca_input, contains("P"))
    if(ncol(subset_pca)<2) {stop("Can't find any columns with PC1..PCX labels")}

    # Grouping samples by POP name and average their PCs values

    ## Aggregate
    pca_aggregated = aggregate(subset_pca[,-c(1)], by= list(subset_pca$POP), FUN = mean)
    names(pca_aggregated)[names(pca_aggregated) == 'Group.1'] <- 'POP'

    ## Selecting Sources and Admixed based on as_file file

    Admixed_pop = as.character(as_file$POP[as_file$A.S == 'A'])
    Source_pop = as.character(as_file$POP[as_file$A.S == 'S'])

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
    pane_RES = evaluate_nnls.mat2(donors = as.matrix(Source_cv[,-c(1)]),recipients = as.matrix(Admixed_cv[,-c(1)]))

    ## Rename cols and rows

    # Cols
    colnames(pane_RES) = append(Source_cv$POP, "Residuals")

    # Rows
    initial_rownames = rownames(pane_RES)
    new_rownames = c()

    for (pane_iteration in 1:length(initial_rownames)) {

      # This will keep only the str after "Target_X"
      keep_str_iteration = sub("Target_[0-9]+","\\1", initial_rownames[pane_iteration])
      get_target_N = gsub(".*_([0-9]+)_.*", "\\1", initial_rownames[pane_iteration])

      new_string = paste0(Admixed_cv$POP[as.numeric(get_target_N)],keep_str_iteration)
      new_rownames = rbind(new_rownames,new_string) }

    rownames(pane_RES) = new_rownames[,1]

    ## printing/returning output matrix
    return(pane_RES)
  }  else if(missing(as_file)&!is.null(sources)&!is.null(admixed)) {

    # Extracting only columns with POP and PCs information
    subset_pca = select(pca_input, contains("P"))
    if(ncol(subset_pca)<2) {stop("Can't find any columns with PC1..PCX labels")}

    # Grouping samples by POP name and average their PCs values

    ## Aggregate
    pca_aggregated = aggregate(subset_pca[,-c(1)], by= list(subset_pca$POP), FUN = mean)
    names(pca_aggregated)[names(pca_aggregated) == 'Group.1'] <- 'POP'

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
    pane_RES = evaluate_nnls.mat2(donors = as.matrix(Source_cv[,-c(1)]),recipients = as.matrix(Admixed_cv[,-c(1)]))

    ## Rename cols and rows

    # Cols
    colnames(pane_RES) <- append(Source_cv$POP, "Residuals")

    # Rows
    initial_rownames = rownames(pane_RES)
    new_rownames = c()

    for (pane_iteration in 1:length(initial_rownames)) {

      # This will keep only the str after "Target_X"
      keep_str_iteration = sub("Target_[0-9]+","\\1", initial_rownames[pane_iteration])
      get_target_N = gsub(".*_([0-9]+)_.*", "\\1", initial_rownames[pane_iteration])

      new_string = paste0(Admixed_cv$POP[as.numeric(get_target_N)],keep_str_iteration)
      new_rownames = rbind(new_rownames,new_string) }

    rownames(pane_RES) = new_rownames[,1]

    ## printing/returning output matrix
    return(pane_RES)
  } else {
    stop('ERROR: optional arguments were not found.')
  }
}
