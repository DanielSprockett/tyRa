#' @title Fits \code{sncm} to an OTU table
#'
#' @description Fits the neutral model from
#' \href{http://onlinelibrary.wiley.com/doi/10.1111/j.1462-2920.2005.00956.x/abstract}{Sloan \emph{et al.} (2006)}
#' to an OTU table and returns several fitting statistics as well as predicted occurrence frequencies for each OTU
#' based on their abundance in the metacommunity. The author of this function is Adam Burns (\email{aburns2@@uoregon.edu}),
#' and was originally published in \href{https://www.nature.com/articles/ismej2015142}{Burns \emph{et al.} (2016)}.
#'
#' @param spp A community table for communities of interest with local communities/samples as rows and taxa as columns.
#' All samples must be rarefied to the same depth.
#'
#' @param pool (optional) A community table for defining source community.
#'
#' @param taxon (optional) A table listing the taxonomic calls for each otu, with OTU ids as row names
#' and taxonomic classifications as columns.
#'
#' @return This function returns list of two elements.
#' The first element, spp.out$fitstats, contains various fitting stats.
#' The second element contains the predicted occurrence frequencies for each OTU/ASV, as well as their \code{fit_class}
#'
#' @seealso
#' \code{\link{plot_sncm_fit}}
#'
#' @examples
#' spp <- otu_table(ps)@.Data
#' spp2 <- otu_table(ps2)@.Data
#' spp.out <- fit_sncm(spp, pool=NULL, taxon=NULL)
#' spp.out <- fit_sncm(spp, pool=spp2, taxon=data.frame(tax_table(ps)))
#'
fit_sncm <- function(spp, pool=NULL, taxon=NULL){

  options(warn=-1)

  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))

  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }

  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]

  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  #Removes rows with any zero (absent in either source pool or local communities)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]

  #Calculate the limit of detection
  d = 1/N

  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.001))
  m.ci <- confint(m.fit, 'm', level=0.95)

  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))

  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))

  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)

  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))

  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

  ##Results
  fitstats <- data.frame(
    m=as.numeric(coef(m.fit)),
    m.ci=as.numeric(coef(m.fit)-m.ci[1]),
    poisLL=as.numeric(pois.mle@details$value),
    Rsqr=as.numeric(Rsqr), # measuring fit, #comparing fit differing datasets to the same model
    Rsqr.pois=as.numeric(Rsqr.pois),
    RMSE=as.numeric(RMSE), # measuring fit #comparing fit differing datasets to the same model
    RMSE.pois=as.numeric(RMSE.pois),
    AIC.pois=as.numeric(aic.pois),  #comparing differing models to the dataset
    BIC.pois=as.numeric(bic.pois), #comparing differing models to the dataset
    N=as.numeric(N),
    Samples=as.numeric(nrow(spp)),
    Richness=as.numeric(length(p)),
    Detect=as.numeric(d))

  A <- cbind(p, freq, freq.pred, pred.ci[,2:3])
  A <- as.data.frame(A)
  colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr')
  if(is.null(taxon)){
    B <- A[order(A[,1]),]
  } else {
    B <- merge(A, taxon, by=0, all=TRUE)
    row.names(B) <- B[,1]
    B <- B[,-1]
    B <- B[order(B[,1]),]
  }
  B <- B[!is.na(B$freq),]
  # fit_class for graphing
  B$fit_class <-"As predicted"
  B[which(B$freq < B$pred.lwr),"fit_class"]<- "Below prediction"
  B[which(B$freq > B$pred.upr),"fit_class"]<- "Above prediction"
  B[which(is.na(B$freq)),"fit_class"]<- "NA"

  # combine fit stats and predicitons into list
  i <- list(fitstats, B)
  names(i) <- c("fitstats", "predictions")
  return(i)
}


