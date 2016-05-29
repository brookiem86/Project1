#' Mantel Haenszel DIF Function for epmr
#' 
#' mhdif
#' 
#' This is a function that identifies differential item functioning using the
#' Mantel Haenszel method
#' 
#' This is a function that identifies differential item functioning (DIF) 
#' for items using the Mantel Haenszel (MH) method. The MH statistic evaluates
#' if the item responses are independent of group membership after taking 
#' into account the performance on the test.
#' 
#' The equation for the MH statistic was obtained from Holland and Thayer
#' (1988).  The value of 0.5 reflects Yates’s correction for continuity and
#' is included in the MH formula.
#' 
#' The common odds ratio is an indication of the degree of association 
#' (Mantel & Haenszel, 1959). Values greater than 1 indicate that on 
#' average the reference group performed better than the focal group on
#' the item. The common odds ratio is commonly transformed into different
#' scales that are symmetric around 0. ETS’s delta scale is a common 
#' transformation used in practice in order to better understand the 
#' magnitude of DIF (Paek & Holland, 2015).
#' 
#' @param idata matrix or data.frame of binary scored item responses
#' @param group numeric vector of values identifying group membership
#' missing values are not allowed for group
#' 
#' @return A data frame containing the MH statistic, common odds ratio, 
#' and delta scale values for each item 
#' 
#' @examples
#' mydata <- ltm::LSAT	
#' g <- runif(1000, min = 0, max = 1)
#' g <- round(g, digits = 0)
#' mydata <- cbind(mydata, g)
#' 
#' mhdif(mydata[ , 1:5], mydata[ , 6]) 
#' 
mhdif <- function(idata, group) {
  vals <- MH.stat <- c.o.ratio <- delta.scale <- out <- NULL
  
  for(item in 1:ncol(idata)) {
    exam.tot <- rowSums(idata, na.rm = TRUE)
    scores <- sort(unique(exam.tot))
    exs <- 1:nrow(idata)
    
    for(a in 1:length(scores)) {
      At <- length(exs[exam.tot == scores[a] & group == 0 & 
                         idata[, item] == 1])
      Bt <- length(exs[exam.tot == scores[a] & group == 0 & 
                         idata[, item] == 0])
      Ct <- length(exs[exam.tot == scores[a] & group == 1 & 
                         idata[, item] == 1])
      Dt <- length(exs[exam.tot == scores[a] & group == 1 & 
                         idata[, item] == 0])
      nRt <- length(exs[exam.tot == scores[a] & group == 0])
      nFt <- length(exs[exam.tot == scores[a] & group == 1])
      nCt <- length(exs[exam.tot == scores[a] & idata[, item] == 1])
      nWt <- length(exs[exam.tot == scores[a] & idata[, item] == 0])
      nT <- length(exs[exam.tot == scores[a]])
      
      vals <- rbind(vals, c(At, Bt, Ct, Dt, nRt, nFt, nCt, nWt, nT))
      colnames(vals) <- c("At", "Bt", "Ct", "Dt", "nRt", "nFt", "nCt", 
                          "nWt", "nT")
    }  
    
    CE.At <- vals[, "nRt"] * vals[, "nCt"] / vals[, "nT"]
    
    var.At.n <- vals[, "nRt"] * vals[, "nFt"] * vals[, "nCt"] * vals[, "nWt"]
    var.At.d <- ((vals[, "nT"])^2) * ((vals[, "nT"]) - 1)
    
    
    MH.stat[[item]] <- ((abs(sum(vals[, "At"] - CE.At)) -.5)^2) /
      sum(var.At.n / var.At.d)
    
    c.o.ratio.n <- sum((vals[, "At"] * vals[, "Dt"]) / vals[, "nT"])
    c.o.ratio.d <- sum((vals[, "Bt"] * vals[, "Ct"]) / vals[, "nT"])        
    c.o.ratio[[item]] <- c.o.ratio.n / c.o.ratio.d
    
    delta.scale[[item]] <- (-2.35) * log(c.o.ratio.n / c.o.ratio.d)
    
  }
  
  out <- data.frame(MH.stat, c.o.ratio, delta.scale)
  out <- round(out, digits = 4)
  
  return(out)
}
