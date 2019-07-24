#' Summarize results from evaluating methods through simulation with function fitModels
#'
#' @param df results data frame as output from function fitModels
#'
#' @return data frame of simulation results averaged over number of data sets
#' 
#' @export
#'
summarizeSimulation <- function(df){
  
  results.tab <- df
  
  colnames(results.tab) <- c("method", "intercept", "slope", "R2", "RMSE", "Cvg", "tsr", "fsr", "brier",
                             "tsr.int", "fsr.int", "brier.int")
  
  lin.res <- results.tab[grep(pattern = ".lin", x = results.tab$method, value = FALSE, fixed = TRUE),]
  nonlin.res <- results.tab[grep(pattern = ".nonlin", x = results.tab$method, value = FALSE, fixed = TRUE),]
  prof.res <- results.tab[grep(pattern = ".prof", x = results.tab$method, value = FALSE, fixed = TRUE),]
  
  # linear exposure-response function
  jack.se <- function(y){
    if(anyNA(y)) y <- y[-which(is.na(y))] # get rid of NA's
    n <- length(y)
    t <- rep(NA, n)
    for(i in 1:n){
      t[i] <- mean(y[-i])
    }
    t.bar <- mean(t)
    jack <- sqrt( ((n-1)/n)*sum((t-t.bar)^2) )
    return(jack)
  }
  
  AME <- apply(lin.res[grep(pattern = "AME", x = lin.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  EMM <- apply(lin.res[grep(pattern = "EMM", x = lin.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  UDR <- apply(lin.res[grep(pattern = "UDR", x = lin.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  PReM <- apply(lin.res[grep(pattern = "PReM", x = lin.res$method, value = FALSE, fixed = TRUE),-1],
                2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  BKMR <- apply(lin.res[grep(pattern = "BKMR", x = lin.res$method, value = FALSE, fixed = TRUE),-1], 
                2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  LM.main <- apply(lin.res[grep(pattern = "main.lm", x = lin.res$method, value = FALSE, fixed = TRUE),-1], 
                   2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  LM.int <- apply(lin.res[grep(pattern = "int.lm", x = lin.res$method, value = FALSE, fixed = TRUE),-1], 
                  2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  linear <- rbind(AME[1,], EMM[1,], UDR[1,], PReM[1,], BKMR[1,], LM.main[1,], LM.int[1,])
  linear.se <- rbind(AME[2,], EMM[2,], UDR[2,], PReM[2,], BKMR[2,], LM.main[2,], LM.int[2,])
  
  
  # non-linear exposure-response function
  AME <- apply(nonlin.res[grep(pattern = "AME", x = nonlin.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  EMM <- apply(nonlin.res[grep(pattern = "EMM", x = nonlin.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  UDR <- apply(nonlin.res[grep(pattern = "UDR", x = nonlin.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  PReM <- apply(nonlin.res[grep(pattern = "PReM", x = nonlin.res$method, value = FALSE, fixed = TRUE),-1],
                2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  BKMR <- apply(nonlin.res[grep(pattern = "BKMR", x = nonlin.res$method, value = FALSE, fixed = TRUE),-1], 
                2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  LM.main <- apply(nonlin.res[grep(pattern = "main.lm", x = nonlin.res$method, value = FALSE, fixed = TRUE),-1], 
                   2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  LM.int <- apply(nonlin.res[grep(pattern = "int.lm", x = nonlin.res$method, value = FALSE, fixed = TRUE),-1], 
                  2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  nonlin <- rbind(AME[1,], EMM[1,], UDR[1,], PReM[1,], BKMR[1,], LM.main[1,], LM.int[1,])
  nonlin.se <- rbind(AME[2,], EMM[2,], UDR[2,], PReM[2,], BKMR[2,], LM.main[2,], LM.int[2,])
  
  
  
  # profile exposure-response function
  AME <- apply(prof.res[grep(pattern = "AME", x = prof.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  EMM <- apply(prof.res[grep(pattern = "EMM", x = prof.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  UDR <- apply(prof.res[grep(pattern = "UDR", x = prof.res$method, value = FALSE, fixed = TRUE),-1], 
               2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  PReM <- apply(prof.res[grep(pattern = "PReM", x = prof.res$method, value = FALSE, fixed = TRUE),-1],
                2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  BKMR <- apply(prof.res[grep(pattern = "BKMR", x = prof.res$method, value = FALSE, fixed = TRUE),-1], 
                2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  LM.main <- apply(prof.res[grep(pattern = "main.lm", x = prof.res$method, value = FALSE, fixed = TRUE),-1], 
                   2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  LM.int <- apply(prof.res[grep(pattern = "int.lm", x = prof.res$method, value = FALSE, fixed = TRUE),-1], 
                  2, FUN = function(x) c(mean(x, na.rm = TRUE), jack.se(x)))
  prof <- rbind(AME[1,], EMM[1,], UDR[1,], PReM[1,], BKMR[1,], LM.main[1,], LM.int[1,])
  prof.se <- rbind(AME[2,], EMM[2,], UDR[2,], PReM[2,], BKMR[2,], LM.main[2,], LM.int[2,])
  
  linear <- round(linear, 2)
  nonlin <- round(nonlin, 2)
  prof <- round(prof, 2)
  
  linear.se <- round(linear.se, 2)
  nonlin.se <- round(nonlin.se, 2)
  prof.se <- round(prof.se, 2)
  
  df.all <- data.frame(rbind(NA, linear, NA, nonlin, NA, prof)) 
  rownames(df.all) <- NULL
  df <- df.all[,(4:10)]
  df <- cbind( rep(c("", "NPBr", 
                     "NPB", 
                     "UPR", 
                     "SPR",
                     "BKMR",
                     "LM",
                     "LM-int"), 3), df)
  colnames(df)[1] <- "method"
  df2 <- cbind(NA, df)
  # df2 <- add_column(df, NA, .before = 1)
  
  df2[,1][which(df2[,1] == "NA")] <- ""
  
  df.print <- df2[,-7] # remove brier scores
  df.print[1,1] <- "linear scenario"
  df.print[9,1] <- "nonlinear scenario"
  df.print[17,1] <- "profiles scenario"
  
  return(df.print)
}

