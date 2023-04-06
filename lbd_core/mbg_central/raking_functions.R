#' Standard Raking Function



#' @title Makes summary rasters from cell pred
#'
#' @description wrapper function for make_cell_pred_summary
#' 
#' @param cell_pred cell pred object 
#' @param simple_raster admin raster matching the cell pred object
#' @param summary_measures charater vector - functions to summarize by. see default below
#' @param raked boolean default `T` - adds "raked" or "unraked" to name of output
#'
#' @return a named list with a raster brick for each summary measure
#' 
make_summary_rasters <- function(cell_pred,
                                 simple_raster,
                                 summary_measures = c('mean','cirange','lower','upper'),
                                 raked = T) {
  
  message('\nMaking summary rasters from cell pred')
  outputlist <-
    lapply(summary_measures, function(x) {
      message("   ", x, ":")
      make_cell_pred_summary(draw_level_cell_pred = cell_pred,
                             mask                 = simple_raster,
                             return_as_raster     = TRUE,
                             summary_stat         = x)
    })
  names(outputlist) <- paste0(summary_measures, ifelse(raked, "_raked", "_unraked"), "_raster")
  
  message("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  message("           Summary Rasters Complete" )
  message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  
  return(outputlist)
}




LogitFindK <- function(gbdval, pixelval, weightval, MaxIter = 40, MaxJump = 10, FunTol = 1e-5, approx_0_1){
  
  # Logit raking won't work with any values of 0 or 1 in cell_pred
  # Adjust values slightly to avoid -Inf or Inf in NewPixelVal
  if (approx_0_1) {
    pixelval[pixelval == 0] <- 1e-10
    pixelval[pixelval == 1] <- 1-(1e-10)
  }
  
  NumIter <- ceiling(-log2(FunTol / MaxJump))
  
  if(NumIter > MaxIter){
    stop(paste("Maximum number of iterations is less than projected iterations required:", NumIter / MaxIter))
  }
  
  CurrentError <- EvalDiff(gbdval, pixelval, weightval)
  if (CurrentError > 0){
    Range <- c(0, MaxJump)
  } else {
    Range <- c(-MaxJump, 0)
  }
  
  a <- Range[1]
  b <- Range[2]
  F_a <- EvalDiff(gbdval, NewPixelVal(a, pixelval), weightval)
  F_b <- EvalDiff(gbdval, NewPixelVal(b, pixelval), weightval)
  
  if (F_a * F_b > 0){
    stop("Your estimates are WAY too far away from GBD")
  } else {
    i <- 1
    c <- (a + b) / 2
    F_c <- EvalDiff(gbdval, NewPixelVal(c, pixelval), weightval)
    Success <- (abs(F_c) <= FunTol) 
    while (!Success & i < NumIter){
      if (sign(F_c) == sign(F_a)){
        a <- c
        F_a <- F_c
      } else {
        b <- c
        F_b <- F_c
      }
      c <- (a + b) / 2
      F_c <- EvalDiff(gbdval, NewPixelVal(c, pixelval), weightval)
      Success <- (abs(F_c) <= FunTol)
      i <- i + 1
    }
    if (Success){
      return(c)
    } else {
      return(sprintf("Need more iterations, current output: K = %g, F(K) = %g",c, F_c))
    }
  }
}


NewPixelVal <- function(K, vals){
  return(ilogit(logit(vals)+K))
}


NewPixelValLog <- function(K, vals){
  return(exp(log(vals)+K))
}


#vals: matrix
#weightval: vector
NewEst <- function(vals, weightval){
  
  vals = vals * weightval
  vals = apply(vals, 2, sum)
  vals = vals / sum(weightval)
  
  return(mean(vals))
}


EvalDiff <- function(gbdval, vals, weightval){
  return(gbdval - NewEst(vals, weightval))
}


logit <- function(x) {
  log(x/(1-x))
}


ilogit <- function(x) {
  exp(x)/(1+exp(x))
}


#overwrite cellIdx function to reduce outward dependancies
cellIdx =function(x) which(!is.na(getValues(x[[1]])))
