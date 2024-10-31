#' pcs_distances
#'
#' pcs_distances allows to estimate the cumulative euclidean distances of the PCs between the Sources, and plot the results (return_plot = 'YES').
#' @param pca_input Dataframe or table with PCA results
#' @param output_name String indicating output name
#' @param sources_file Two columns table: S1 and S2. Contains the pairs of sources that will be compared.
#' @param return_plot NULL by default or 'YES' to plot with ggplot2
#' @examples
#' \dontrun{
#' Source_D = read.table('data/Sources_Distances', header =T)
#' pca = read_flash(pca_input = 'data/TOY_flash.pca')
#' pcs_distances(pca, 'data/output_distances', sources_file = Source_D, return_plot = 'YES')
#' }
#' @return Returns an output_name.csv file with the cumulative distances
#' @return Returns a output_name.pdf file with the cumulative distances
#' @export

pcs_distances <- function(pca_input, output_name, sources_file, return_plot = NULL) {

### GROUP BY POP
  subset_pca = select(pca_input, contains("P"))
  if(ncol(subset_pca)<2) {stop("Can't find any column with PC1..PCX labels, check your header.")}
  pca_aggregated = aggregate(subset_pca[,-c(1)], by= list(subset_pca$POP), FUN = mean)
  names(pca_aggregated)[names(pca_aggregated) == 'Group.1'] <- 'POP'
  grouped = pca_aggregated
  grouped_reassigned <- grouped[,-1]
  rownames(grouped_reassigned) <- grouped[,1]

### Estimate N of PCS
  x = select(pca, contains("PC"))
  n_pcs = ncol(x)

### Loop over pcs

distances = data.frame(nrows=n_pcs-1)
for (pop in 1:nrow(sources_file)) {
  print(paste0("Calculating euclidean distances between: ",sources_file$S1[pop]," ",sources_file$S2[pop]))
  sources = grouped_reassigned[c(as.character(sources_file$S1[pop]),as.character(sources_file$S2[pop])),]
  if(any(is.na(sources))) {warning("Cannot find PCs for: ",sources_file$S1[pop]," ",sources_file$S2[pop])}

  df = data.frame()
  for (pcs in 2:n_pcs){
    # D is estimated from PC1 to PCS
    d = dist(sources[,1:pcs], method = "euclidean")
    #if (is.na(d)) {print(paste0("Your distances are NA, check ", sources_file$S1[pop]," ",sources_file$S2[pop]," in your source/PCA file"))}
    #print(paste0("PC1-PC",pcs," ",d))
    df[pcs-1,1] <- c(paste0("PC1-PC",pcs))
    df[pcs-1,2] <- d
  }

  colnames(df) <- c("PCs",paste0(as.character(sources_file$S1[pop]),"-",as.character(sources_file$S2[pop])))
  distances = cbind(distances,df)
}

distances_sub = distances[,-1]

dupl = duplicated(as.list(distances_sub))
distances_nodupl = distances_sub[!dupl]
final_d <- distances_nodupl %>% select_if(~ !any(is.na(.)))

write.table(final_d, file = paste0(output_name,".csv"), row.names = FALSE)
print(paste0("Eulidean distances csv file in: ", output_name, ".csv"))

if (!is.null(return_plot)) {

### PLOT

final_d$PCs.ordered <- factor(final_d$PCs, levels=final_d$PCs)

# Open the PDF device at some file
plot_list = list()

### Thanks Gregor Thomas on Stackoverflow
### https://stackoverflow.com/questions/52919899/ggplot2-display-every-nth-value-on-discrete-axis
### For this function to display every nth value on discrete axis

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

for (i in (2:(ncol(final_d)-1))) {
p <- ggplot(data=final_d)+
  geom_point(aes(x=PCs.ordered, y=final_d[,i])) +
  ggtitle(as.character(colnames(final_d)[i]))+
  xlab("PCs") + ylab("Cumulative distance between PCs") +
  theme(axis.text.x = element_text(size = 4, angle = 45, vjust = 1, hjust=1)) +

#  coord_fixed(ratio=400) +

### Gregor Thomas on Stackoverflow
  scale_x_discrete(breaks = every_nth(n = 2)) +
  NULL

plot_list[[i]] = p
#print(as.character(colnames(final_d)[i]))
}

print(paste0("Eulidean distances plot in: ", output_name, ".pdf"))

pdf(paste0(output_name,".pdf"))
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off() }
}

