#' plot_pane
#'
#' plot_pane allows to plot pane assignations as barplot, using either R base (default) or ggplot2 (type_ggplot = 'YES')
#' @param pane_result The matrix returned by pane()
#' @param output_name output name for pane pdf plot
#' @param type_ggplot NULL by default or 'YES' to plot with ggplot2
#' @return Returns an pane_plot.pdf file in the working directory
#' @export

plot_pane <- function(pane_result, output_name, type_ggplot = NULL) {

  #### R BASE PLOT ####

  if(is.null(type_ggplot)) {

    # Creating an empty pdf
    pdf(file = paste0(output_name,'.pdf'))
    par(mar=c(6, 4, 4.5, 5.2), xpd=TRUE)

    # Creating a barplot
    x = barplot(t(pane_result[[1]]),
                main="pane ancestry assignations",
                #las = 2,
                #cex.names = .6,
                xaxt = 'n',
                #xlab = 'Target groups',
                ylab = 'Proportion of ancestry',
                col = terrain.colors(nrow(pane_result[[1]])))

    # Creating the label list for the x axis
    labs <- paste(names(table(rownames(pane_result[[1]]))))

    # Rotating the labels in the x axis 45°
    text(cex=0.7, x=x-1.5, y=-0.12, labs, xpd=TRUE, srt=45)

    # Creating the legend
    legend('topright', inset = c(-0.16,0),
           legend=colnames(pane_result[[1]]),
           fill=terrain.colors(nrow(pane_result[[1]])))

    dev.off() } else if(type_ggplot == 'YES') {

  #### GGPLOT2 PLOT ####

      # Reading pane() output as a dataframe, while transposing it
      z = as.data.frame(t(pane_result[[1]]))

      # Converting the dataframe in a 3 column dataframe, for easier plotting
      converted_dataframe = NULL

      for (i in 1:ncol(z)) {
        tmp = z[,i, drop=FALSE]
        tmp$Target = colnames(tmp)
        tmp$Donors <- rownames(tmp)
        colnames(tmp) = c('Proportion','Target','Donors')
        rownames(tmp) <- NULL
        converted_dataframe <- rbind(converted_dataframe,tmp)
      }

      # Plotting the converted dataframe with ggplot2 as barplot

      ggplot_pane_plot <- ggplot(data = converted_dataframe, aes(x=Target, y=Proportion, fill=Donors)) +
        geom_bar(stat='identity') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        NULL

      #Saving the plot in a pdf format

      ggsave(paste0(output_name,'.pdf'))
    } else {
      print("Please indicate either type_ggplot = NULL or type_ggplot = 'YES' as optional argument ")
    }
}


