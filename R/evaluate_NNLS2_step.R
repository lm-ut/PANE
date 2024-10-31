evaluate_getoverall_tentativefits=function(predmat,fitdata){
  restrict=1
  rep=1
  i=1
  t = 1
  all_fits = list()
  while(rep==1){
    q=evaluate_getfit(predmat,fitdata,restrict=i)
    fits=evaluate_get_alltentativefits(predmat,fitdata,restrict=i)
    all_fits[[t]] <- fits
    t = t+1
    
    if(q$x[i]>0) rep=0
    i=i+1
  }
  
  return(all_fits)

}

evaluate_get_alltentativefits=function(predmat,fitdata,restrict=1){
  empty_list=list()
  temp=matrix(predmat[-restrict,],ncol=dim(predmat)[2])
  for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]
  fitdata2=fitdata-predmat[restrict,]
  v=nnls(t(temp),fitdata2)
  x=v$x
  newx=1:nrow(predmat)
  newx[!((1:nrow(predmat))==restrict)]=x
  newx[restrict]=1-sum(x)
  v$x=newx
  
  nnls_fits_residuals = append(v$x, v$deviance)
  empty_list <- nnls_fits_residuals
  return(empty_list)
}

evaluate_getfit=function(predmat,fitdata,restrict=1){
  temp=matrix(predmat[-restrict,],ncol=dim(predmat)[2])
  for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]
  fitdata2=fitdata-predmat[restrict,]
  v=nnls(t(temp),fitdata2)
  x=v$x
  
  newx=1:nrow(predmat)
  newx[!((1:nrow(predmat))==restrict)]=x
  newx[restrict]=1-sum(x)
  
  v$x=newx
  names(v$x)=rownames(predmat)
  return(v)
}


#' evaluate_nnls.mat2
#'
#' evaluate_nnls.mat2 function solves nonnegative least squares problems. It requires two matrices, one ('donors') refers to the source groups, the second ('recipients') refers to the admixed groups
#' @param donors Matrix with reference groups
#' @param recipients Matrix with target groups
#' @examples
#' \dontrun{
#' nnls.mat2(donors = my_source_individuals_matrix,recipients = my_admixed_individuals_matrix)
#' }
#' @return Returns matrix describing the admixed groups as a mixture of the source groups, along with the residuals
#' @export

evaluate_nnls.mat2 <- function(donors,recipients){
  
  q = list()
  r=c()
  for (rec in 1:nrow(recipients)){
    q[[rec]]=evaluate_getoverall_tentativefits(predmat=donors,fitdata=recipients[rec,]) }
  
  for (target in 1:length(q)) {
    for (tentative_fits_per_target in 1:length(q[[target]])) {
      iter_text = paste0("Target_", target, "_iteration_", tentative_fits_per_target)
      iter_rows = append(iter_text, q[[target]][tentative_fits_per_target][[1]])
      r = rbind(r,iter_rows)
    }
  }
  
  # Reassigning rownames, now matrix is r_rownames
  r_rownames <- r[,-1]

  if (is.null(dim(r_rownames))) {
    r_mat = matrix(r[,-c(1)],nrow=1)
    rownames(r_mat) = "Target_1_iteration_1"
    r_rownames_int <- array(as.numeric(sub(",",".",r_mat)), dim(r_mat), dimnames(r_mat))
    return (r_rownames_int)
    #stop("There are no multiple models to evaluate. Please simply run pane() or add source or target groups.")

  } else {
    rownames(r_rownames) <- r[,1]
    r_rownames_int <- array(as.numeric(sub(",", ".", r_rownames)), dim(r_rownames), dimnames(r_rownames))
    return (r_rownames_int)

  }
}
