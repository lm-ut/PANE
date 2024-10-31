#' read_eigen
#'
#' read_eigen reads smartpca output from EIGENSTRAT, where the first column contains Family ID, the second sample ID, and the following the PCs
#' @param pca_input points to the directory and pca.evec file
#' @examples
#' \dontrun{
#' read_eigen(pca_input = 'data/TOY.pca.evec')
#' }
#' @return Returns a PCA matrix with new header: POP ID PC1 PCN CC
#' @export

read_eigen <- function(pca_input) {
pca = read.table(text = gsub(":", " ", readLines(pca_input)))
pcx <- c()
for(i in 1:ncol(pca[,-c(1,2,ncol(pca))])) {
    pcx[[i]] <- paste0("PC",i)}
colnames(pca)=c("POP","IND",pcx,"CC")
return(pca) }

#' read_flash
#'
#' read_flash reads PCA output from flashpca, where the first column contains Family ID, the second sample ID, and the following the PCs
#' @param pca_input point to the directory and flash pca file
#' @examples
#' \dontrun{
#' read_flash(pca_input = 'data/TOY_flash.pca')
#' }
#' @return Returns a PCA matrix with new header: POP ID PC1 PCN
#' @export

read_flash <- function(pca_input) {
  pca = read.table(pca_input, header=TRUE, sep="")
  pcx <- c()
  for(i in 1:ncol(pca[,-c(1,2)])) {
    pcx[[i]] <- paste0("PC",i)}
  colnames(pca)=c("POP","IND",pcx)
  return(pca)
}

