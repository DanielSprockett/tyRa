#' @title plots the fit of an OTU table to \code{sncm}
#'
#' @description This funtion plots the output from \code{\link{fit_sncm}}

#' @param spp_out the output from \code{\link{fit_sncm}}
#'
#' @param fill (optional) can either be set to fit_class to color by prediction (default), or by a taxonomic level
#' available in \code{\link[phyloseq]{rank_names}}
#'
#' @param title (optional) the title of the plot.
#'
#' @return This function returns a plot of the fit of the OTU table to sncm
#'
#' @seealso
#' \code{\link{fit_sncm}}
plot_sncm_fit <- function(spp_out = spp_out, fill = NULL, title = NULL){

  tax_levels <- colnames(spp_out$predictions)[7:length(colnames(spp_out$predictions))-1]

  if(is.null(fill)){
    fill <- "fit_class"
  }

  r2_val <- paste("r^2 ==", round(spp_out$fitstats$Rsqr,4))
  m_val <- paste("m ==", round(spp_out$fitstats$m,4))
  df <- data.frame(t(table(spp_out$predictions$fit_class)))
  df <- df[,c(2,3)]
  colnames(df) <- c("Prediction", "AVS Abundance")

  p <- ggplot(data=spp_out$predictions)

  if(fill == "fit_class"){
    p <- p + geom_point(aes(x = log(p), y = freq, fill=eval(parse(text=fill))), shape =21, color="black", size =2, alpha=0.75)
    p <- p + scale_fill_manual(
      name = "Prediction",
      values = c("Above prediction" = "cyan4", "As predicted" = "blue", "Below prediction" = "goldenrod", "NA" = "white"),
      breaks = c("Above prediction", "As predicted", "Below prediction", "NA"),
      labels = c(paste0("Above prediction (",round((df[1,2]/spp_out$fitstats$Richness)*100, 1),"%)"),
                 paste0("As predicted (",round((df[2,2]/spp_out$fitstats$Richness)*100, 1),"%)"),
                 paste0("Below Prediction (",round((df[3,2]/spp_out$fitstats$Richness)*100, 1),"%)"),
                 paste0("NA (",df[4,2],")")))

  }else if (fill %in% tax_levels){
    p <- p + geom_point(aes(x = log(p), y = freq, fill=eval(parse(text=fill))), shape =21, color="black", size =2, alpha=0.75)
    p <- p + scale_fill_discrete(name = "Taxon")

  } else{
    print(paste0("fill variable: ", fill, " is not a valid taxonomic level or fit_class"))
  }

  p <- p + geom_line(aes(x = log(p), y = freq.pred), color = "blue")
  p <- p + geom_line(aes(x = log(p), y = pred.lwr), color = "blue", linetype="dashed")
  p <- p + geom_line(aes(x = log(p), y = pred.upr), color = "blue", linetype="dashed")
  p <- p + xlab("log(Mean Relative Abundance)")
  p <- p + ylab("Frequency")
  p <- p + ggtitle(title)
  p <- p + annotate("text", x=mean(log(spp_out$predictions$p), na.rm = TRUE), y=0.95, size=5, label = r2_val, parse=TRUE)
  p <- p + annotate("text", x=mean(log(spp_out$predictions$p), na.rm = TRUE), y=0.9, size=5, label = m_val, parse=TRUE)

  return(p)
}



