getfit=function(predmat,fitdata,restrict=1){

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

getoverallfit=function(predmat,fitdata){
  restrict=1
  rep=1
  i=1
  while(rep==1){
    q=getfit(predmat,fitdata,restrict=i)

    if(q$x[i]>0) rep=0
    i=i+1
  }

  return(q)
}

newgetfit <- function(predmat,fitdata,restrict=1){
  temp=predmat[-restrict,]
  for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]
  fitdata2=fitdata-predmat[restrict,]
  v=nnls(t(temp),fitdata2)
  x=v$x
  newx=1:nrow(predmat)
  newx[-restrict]=x
  newx[restrict]=1-sum(x)
  v$x=newx

  return(v)}


getoverallfit2=function(predmat,fitdata){
  i=1
  q=newgetfit(predmat,fitdata,restrict=i)
  while(q$x[i]<0){
    print(q)
    q=newgetfit(predmat,fitdata,restrict=i)
    i=i+1
    print(q)
  }
  names(q$x)=rownames(predmat)
  return(q)
}

#' nnls.mat2
#'
#' nnls.mat2 function solves nonnegative least squares problems. It requires two matrices, one ('donors') refers to the source groups, the second ('recipients') refers to the admixed groups
#' @param donors Matrix with reference groups
#' @param recipients Matrix with target groups
#' @examples
#' \dontrun{
#' nnls.mat2(donors = my_source_individuals_matrix,recipients = my_admixed_individuals_matrix)
#' }
#' @return Returns matrix describing the admixed groups as a mixture of the source groups, along with the residuals
#' @export

nnls.mat2 <- function(donors,recipients){
  mat <- matrix(nrow=nrow(recipients),ncol=nrow(donors),dimnames =list(rownames(recipients),rownames(donors)))
  r=c()
  for (i in 1:nrow(mat)){
    # changed getoverallfit2 with getoverallfit
    q=getoverallfit(predmat=donors,fitdata=recipients[i,])
    mat[i,]=q$x
    r[i]=q$deviance}
  return (list(mat,r))
}


# nnls.mat  function solves nonnegative least squares problems. It requires two matrices, one ('donors') refers to the reference groups, the second ('recipients') refers to the target groups
# donors Matrix with reference groups
# recipients Matrix with target groups
# nnls.mat(donors = my_source_individuals_matrix,recipients = my_admixed_individuals_matrix)
# Returns matrix describing the target groups as a mixture of the reference groups, WITHOUT the residuals

nnls.mat=function(donors,recipients){
  tmp=t(sapply(1:nrow(recipients),function(X) getoverallfit2(predmat=donors,fitdata=recipients[X,])$x))
  rownames(tmp)=rownames(recipients)
  return(tmp)
}


