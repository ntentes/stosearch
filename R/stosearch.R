################################################################################################
################################################################################################
######## Stochastic Search Functions - An R Package by Konstantinos Ntentes ####################
################################################################################################
################################################################################################



################################################################################################
# Simulate standard brownian motion with mu = 0 and sigma^2 = t
################################################################################################


## input: a positive integer n, a positive integer t
## output: a vector of length n+1 modelling the movement of
##         standard brownian motion w distr. N(0, t), B(0)=0

SBM_gen <- function(n, start_y = 0, t = 1){
  sbm <- c(start_y, start_y+cumsum(rnorm(n, sd = sqrt(t/n))))
  return(sbm)
}



################################################################################################
## Search for an integer target starting at 0 (one searcher)
################################################################################################


SBM_search_pro <- function(x0 = 20, t = 1000, n = 1e5){
  vec <- SBM_gen(n, start_y = x0, t)
  steps <- -1
  if (x0 <= 0) {
    targ_hit <- vec >= 0
    hit <- which(targ_hit)
    if (length(hit) == 0){
      result <- list("TargetFound"=FALSE,
                     "TimeUnits"=t,
                     "SearchPath"=vec[1:n])
    } else {
      steps <- hit[1]
      cut_vec <- vec[1:steps]
      result <- list("TargetFound"=TRUE,
                     "TimeUnits"=(length(cut_vec)*(t/n)),
                     "SearchPath"=c(cut_vec, rep(NaN, (n-steps))))
    }
  } else {
    targ_hit <- vec <= 0
    hit <- which(targ_hit)
    if (length(hit) == 0){
      result <- list("TargetFound"=FALSE,
                     "TimeUnits"=t,
                     "SearchPath"=vec[1:n])
    } else {
      steps <- hit[1]
      cut_vec <- vec[1:steps]
      result <- list("TargetFound"=TRUE,
                     "TimeUnits"=(length(cut_vec)*(t/n)),
                     "SearchPath"=c(cut_vec, rep(NA, (n-steps))))
    }
  }
  return(result)
}


################################################################################################
## stochastic search with poisson resetting (one searcher) at times according to an
## exponential distribution at optimal rate 1/(2*x0^2)
################################################################################################


SBM_search_SR_pro <- function(x0 = 20, t = 1000, n = 1e5, rate = 'optimal'){
  if (rate == 'optimal'){
    reset_dist_rate <- 1/(2*(x0^2))
  } else {
    reset_dist_rate <- rate
  }
  reset_times <- rexp(n = 1, rate = reset_dist_rate)
  while (sum(reset_times)<t){
    reset_times <- c(reset_times, rexp(n = 1, rate = reset_dist_rate))
  }
  reset_steps <- round((reset_times*n)/t)
  reset_steps_copy <- reset_steps
  vec <- SBM_gen(n, start_y = x0, t)[1:reset_steps[1]]
  while (length(vec)<n){
    if (is.na(reset_steps[2])==FALSE){
      reset_steps <- reset_steps[-1]
    }
    cutoff <- reset_steps[1]
    new_bm <- SBM_gen(n, start_y = x0, t)
    vec <- c(vec, new_bm[1:cutoff])
  }
  vec <- vec[1:(n+1)]
  steps <- -1
  if (x0 <= 0) {
    targ_hit <- vec >= 0
    hit <- which(targ_hit)
    if (length(hit) == 0){
      num_resets <- length(reset_times)-1
      result <- list("TargetFound" = FALSE,
                     "TimeUnits" = t,
                     "NumResets" = num_resets,
                     "SearchPath" = vec[1:n],
                     "ResetTimes" = if (num_resets==0){
                       0
                     } else {
                       cumsum(reset_times[1:num_resets])
                     })
    } else {
      steps <- hit[1]
      cut_vec <- vec[1:steps]
      loc_resets <- which((steps<=cumsum(reset_steps_copy)))
      num_resets <- loc_resets[1]-1
      result <- list("TargetFound" = TRUE,
                     "TimeUnits" = (length(cut_vec)/n)*t,
                     "NumResets" = num_resets,
                     "SearchPath" = c(cut_vec, rep(NA, (n-steps))),
                     "ResetTimes" = if (num_resets==0){
                       0
                     } else {
                       cumsum(reset_times[1:num_resets])
                     })
    }
  } else {
    targ_hit <- vec <= 0
    hit <- which(targ_hit)
    if (length(hit) == 0){
      num_resets <- length(reset_times)-1
      result <- list("TargetFound" = FALSE,
                     "TimeUnits" = t,
                     "NumResets" = num_resets,
                     "SearchPath" = vec[1:n],
                     "ResetTimes" = if (num_resets==0){
                       0
                     } else {
                       cumsum(reset_times[1:num_resets])
                     })
    } else {
      steps <- hit[1]
      cut_vec <- vec[1:steps]
      loc_resets <- which((steps<=cumsum(reset_steps_copy)))
      num_resets <- loc_resets[1]-1
      result <- list("TargetFound" = TRUE,
                     "TimeUnits" = (length(cut_vec)/n)*t,
                     "NumResets" = num_resets,
                     "SearchPath" = c(cut_vec, rep(NA, (n-steps))),
                     "ResetTimes" = if (num_resets==0){
                       0
                     } else {
                       cumsum(reset_times[1:num_resets])
                     })
    }
  }
  return(result)
}



################################################################################################
## Search for an integer target starting at 0 (one searcher), with deterministic resetting
################################################################################################


SBM_search_DR_pro <- function(x0 = 20, t = 1000, n = 1e5, det_reset = 200){
  terms_before_reset <- (n/t)*det_reset
  vec <- SBM_gen(n, start_y = x0, t)[1:(terms_before_reset+1)]
  while (length(vec)<n){
    vec <- c(vec, SBM_gen(n, start_y = x0, t)[1:terms_before_reset])
  }
  steps <- -1
  if (x0 <= 0) {
    vec <- vec[1:n]
    targ_hit <- vec >= 0
    hit <- which(targ_hit)
    if (length(hit) == 0){
      result <- list("TargetFound" = FALSE,
                     "TimeUnits" = t,
                     "NumResets" = ifelse(t%%det_reset==0, t/det_reset-1, floor(t/det_reset)),
                     "SearchPath" = vec[1:n])
    } else {
      steps <- hit[1]
      cut_vec <- vec[1:steps]
      path_vec <- rep(NA, n)
      result <- list("TargetFound" = TRUE,
                     "TimeUnits" = (t/n)*steps,
                     "NumResets" = floor((steps/n)/det_reset*t),
                     "SearchPath" = c(cut_vec, path_vec)[1:n])
    }
  } else {
    vec <- vec[1:n]
    targ_hit <- vec <= 0
    hit <- which(targ_hit)
    if (length(hit) == 0){
      result <- list("TargetFound" = FALSE,
                     "TimeUnits" = t,
                     "NumResets" = ifelse(t%%det_reset==0, t/det_reset-1, floor(t/det_reset)),
                     "SearchPath" = vec[1:n])
    } else {
      steps <- hit[1]
      cut_vec <- vec[1:steps]
      path_vec <- rep(NA, n)
      result <- list("TargetFound" = TRUE,
                     "TimeUnits" = (t/n)*steps,
                     "NumResets" = floor((steps/n)/det_reset*t),
                     "SearchPath" = c(cut_vec, path_vec)[1:n])
    }
  }
  return(result)
}



################################################################################################
## Search for an integer target starting at 0 (N searchers), no graphics
################################################################################################



multisearch_sbm_pro <- function(x0 = 20, t = 1000, n = 1e5, searchers = 2){
  search_array <- numeric()
  steps <- 0
  for (i in 1:searchers){
    iter_search <- SBM_search_pro(x0, t, n)
    time_to_n <- round(iter_search$TimeUnits * (n/t))
    search_array <- rbind(search_array, iter_search$SearchPath)
    if ((iter_search$TargetFound == TRUE)&&(steps==0)){
      steps <- time_to_n
    } else if ((iter_search$TargetFound == TRUE)&&(steps>time_to_n)){
      steps <- time_to_n
    }
  }
  if (steps==0){
    result <- list("TargetFound" = FALSE,
                   "TimeUnits" = t,
                   "SearchCost" = t*searchers,
                   "SearchPaths" = search_array,
                   "NumResets" = NA,
                   "ResetTimes" = NA)
  } else {
    search_array[ , (steps+1):n] <- 0L
    result <- list("TargetFound" = TRUE,
                   "TimeUnits" = steps*(t/n),
                   "SearchCost" = steps*(t/n)*searchers,
                   "SearchPaths" = search_array,
                   "NumResets" = NA,
                   "ResetTimes" = NA)
  }
  return(result)
}




################################################################################################
## Search for an integer target starting at 0 (N searchers)
## with stochastic resetting, no graphics
################################################################################################




multisearch_SR_pro <- function(x0 = 20, t = 1000, n = 1e5, rate = 'optimal', searchers = 2){
  if (rate == 'optimal'){
    rate <- 1/(2*(x0^2))
  }
  search_array <- numeric()
  steps <- 0
  reset_times <- numeric()
  num_resets <- 0
  for (i in 1:searchers){
    iter_search <- SBM_search_SR_pro(x0, t, n, rate)
    time_to_n <- round(iter_search$TimeUnits * (n/t))
    search_array <- rbind(search_array, iter_search$SearchPath)
    if (iter_search$NumResets!=0){
      reset_times <- c(reset_times, iter_search$ResetTimes)
      num_resets <- num_resets + iter_search$NumResets
    }
    if ((iter_search$TargetFound == TRUE)&&(steps==0)){
      steps <- time_to_n
    } else if ((iter_search$TargetFound == TRUE)&&(steps>time_to_n)){
      steps <- time_to_n
    }
  }
  if (steps==0){
    result <- list("TargetFound" = FALSE,
                   "TimeUnits" = t,
                   "SearchCost" = t*searchers,
                   "SearchPaths" = search_array,
                   "NumResets" = num_resets,
                   "ResetTimes" = reset_times)
  } else {
    reset_times <- reset_times[reset_times<=(steps*(t/n))]
    search_array[ , (steps+1):n] <- 0L
    result <- list("TargetFound" = TRUE,
                   "TimeUnits" = steps*(t/n),
                   "SearchCost" = steps*(t/n)*searchers,
                   "SearchPaths" = search_array,
                   "NumResets" = length(reset_times),
                   "ResetTimes" = reset_times)
  }
  return(result)
}


################################################################################################
## Search for an integer target starting at 0 (N searchers)
## with deterministic resetting, no graphics
################################################################################################



multisearch_DR_pro <- function(x0 = 20, t = 1000, n = 1e5, det_reset = 200, searchers = 2){
  search_array <- numeric()
  steps <- 0
  for (i in 1:searchers){
    iter_search <- SBM_search_DR_pro(x0, t, n, det_reset)
    time_to_n <- round(iter_search$TimeUnits * (n/t))
    search_array <- rbind(search_array, iter_search$SearchPath)
    if ((iter_search$TargetFound == TRUE)&&(steps==0)){
      steps <- time_to_n
    } else if ((iter_search$TargetFound == TRUE)&&(steps>time_to_n)){
      steps <- time_to_n
    }
  }
  if (steps==0){
    result <- list("TargetFound" = FALSE,
                   "TimeUnits" = t,
                   "SearchCost" = t*searchers,
                   "SearchPaths" = search_array,
                   "NumResets" = floor(t/det_reset)*searchers,
                   "ResetTimes" = if ((floor(t/det_reset))>=1) {
                     (seq(from=1, to=(floor(t/det_reset)), by=1))*det_reset
                   } else {
                     NA
                   })
  } else {
    search_array[ , (steps+1):n] <- 0L
    result <- list("TargetFound" = TRUE,
                   "TimeUnits" = steps*(t/n),
                   "SearchCost" = steps*(t/n)*searchers,
                   "SearchPaths" = search_array,
                   "NumResets" = floor(((steps/n)*t)/det_reset)*searchers,
                   "ResetTimes" = if (floor(((steps/n)*t)/det_reset)==0){
                     NA
                   } else {
                     seq(from=1, to=floor(((steps/n)*t)/det_reset), by=1)*det_reset
                   })
  }
  return(result)
}


################################################################################################
## Search for an integer target starting at 0 (N searchers), with graphics
################################################################################################


multisearch_sbm <- function(x0 = 20, t = 1000, n = 1e5, searchers = 2, output_data=FALSE){
  msearch <- multisearch_sbm_pro(x0, t, n, searchers)
  steps <- msearch$TimeUnits * (n/t)
  if (searchers==1){
    plot_matrix <- msearch$SearchPaths[ 1:searchers, 1:steps]
  } else{
    plot_matrix <- t(msearch$SearchPaths[ 1:searchers, 1:steps])
  }
  t_axis <- seq(from = 0, to = (steps-1)*t/n, by = t/n)
  clim <- if(x0<=0) {
    c(floor(min(plot_matrix)),  as.numeric(0+1))
  } else {
    c(as.numeric(0-1), ceiling(max(plot_matrix)))
  }
  matplot(x=t_axis, y=plot_matrix, type = 'l', col = 1:searchers,
          xlim = c(0, t), ylim = clim, lty = 1,
          xlab = "Time Units", ylab = "Distance",
          main = "Stochastic Search: Standard")
  abline(h=0, col = 'red', lty = 'dashed')
  if (msearch$TargetFound == TRUE){
    abline(v=msearch$TimeUnits, col = 'darkolivegreen', lty = 'dashed')
  }
  if (output_data==TRUE){
    return(msearch)
  }
}


################################################################################################
## Search for an integer target starting at 0 (N searchers)
## with stochastic resetting, with graphics
################################################################################################


multisearch_SR <- function(x0 = 20, t = 1000, n = 1e5, rate = 'optimal', searchers = 2, output_data = FALSE){
  if (rate == 'optimal'){
    rate <- 1/(2*(x0^2))
  }
  msearch <- multisearch_SR_pro(x0, t, n, rate, searchers)
  steps <- msearch$TimeUnits * (n/t)
  if (searchers==1){
    plot_matrix <- msearch$SearchPaths[ 1:searchers, 1:steps]
  } else{
    plot_matrix <- t(msearch$SearchPaths[ 1:searchers, 1:steps])
  }
  t_axis <- seq(from = 0, to = (steps-1)*t/n, by = t/n)
  clim <- if(x0<=0) {
    c(floor(min(plot_matrix)),  as.numeric(0+1))
  } else {
    c(as.numeric(0-1), ceiling(max(plot_matrix)))
  }
  matplot(x=t_axis, y=plot_matrix, type = 'l', col = 1:searchers,
          xlim = c(0, t), ylim = clim, lty = 1,
          xlab = "Time Units", ylab = "Distance",
          main = "Stochastic Search: Stochastic Resetting")
  abline(h=0, col = 'red', lty = 'dashed')
  if (msearch$TargetFound == TRUE){
    abline(v=msearch$TimeUnits, col = 'darkolivegreen', lty = 'dashed')
  }
  if (msearch$NumResets!=0){
    abline(v=msearch$ResetTimes, col = 388, lty = 3)
  }
  if (output_data==TRUE){
    return(msearch)
  }
}


################################################################################################
## Search for an integer target starting at 0 (N searchers)
## with deterministic resetting, with graphics
################################################################################################


multisearch_DR <- function(x0 = 20, t = 1000, n = 1e5, det_reset = 'optimal', searchers = 2, output_data = FALSE){
  if (det_reset=='optimal'){
    det_reset <- round(0.5*(x0^2)/0.458)
  }
  msearch <- multisearch_DR_pro(x0, t, n, det_reset, searchers)
  steps <- msearch$TimeUnits * (n/t)
  if (searchers==1){
    plot_matrix <- msearch$SearchPaths[ 1:searchers, 1:steps]
  } else {
    plot_matrix <- t(msearch$SearchPaths[ 1:searchers, 1:steps])
  }
  t_axis <- seq(from = 0, to = (steps-1)*t/n, by = t/n)
  clim <- if(x0<=0) {
    c(floor(min(plot_matrix)),  as.numeric(0+1))
  } else {
    c(as.numeric(0-1), ceiling(max(plot_matrix)))
  }
  matplot(x=t_axis, y=plot_matrix, type = 'l', col = 1:searchers,
          xlim = c(0, t), ylim = clim, lty=1,
          xlab = "Time Units", ylab = "Distance",
          main = "Stochastic Search: Deterministic Resetting")
  abline(h=0, col = 'red', lty = 'dashed')
  if (msearch$TargetFound == TRUE){
    abline(v=msearch$TimeUnits, col = 'darkolivegreen', lty = 'dashed')
  }
  reset_lines <- if (ifelse(t%%det_reset==0, t/det_reset-1, floor(t/det_reset))==0){
    0
  } else {1:(ifelse((t%%det_reset)==0, t/det_reset-1, floor(t/det_reset)))
  }
  if (reset_lines[1]!=0){
    abline(v=(reset_lines*det_reset), col = 388, lty = 3)
  }
  if (output_data==TRUE){
    return(msearch)
  }
}




################################################################################################
################################################################################################
############ Culminating Function - Main Search Function #######################################
################################################################################################
################################################################################################



#' Stochastic Search with resetting
#'
#' Allows a user to specify a positive starting point and simulate (data or visually or both) 1 to N searchers moving according to standard brownian motion ($\mu = 0$ and $\sigma^2 = t$) with or without deterministic or stochastic resetting according to rates set by the user. Search will stop when any searcher finds the target (the origin).
#'
#' @param x0 positive numeric initial value for all searchers (vector of length 1) (searchers will also reset to this value if resetting is chosen)
#' @param t time units alotted to searchers to complete the search
#' @param searchers number of searchers (integer from 1 to N) (warning: simulation will take a long time with a large number of searchers)
#' @param reset_type type of resetting, argument accepts "none", "stochastic" (each searcher is independently and randomly reset according to a poisson process with specified rate) or "deterministic" (all searchers are reset deterministically according to the specified rate).
#' @param reset_rate reciprocal of mean reset time. Accepts "optimal" or a positive numeric vector of length 1. "optimal" is rate = 1/(2*(x0)^2) which is the optimal resetting rate according to Bhat, De Bacco, Redner (2016).
#' @param n precision of brownian motion, as $n \rightarrow \inf$, brownian motion is calculated in continuous time.
#' @param output_data accepts TRUE or FALSE, if false, does not output list of collected data re: simulation.
#' @param output_visual accepts TRUE or FALSE, if false, does not output visual plot of the search
#'
#' @return if output_data==TRUE, then function returns a list of named vectors: "TargetFound" = T or F, "TimeUnits" = time elapsed until a searcher finds the target, "SearchCost" = TimeUnits * # of searchers, "SearchPaths" = a matrix of numerics where the rows (of length n) are the search path of each searcher, "NumResets" = number of times a searcher was reset (provides NA is no resetting is selected), "ResetTimes" = a vector of reset times (NA is provided if no resetting is selected). If output_visual==TRUE, a plot of the search is printed, with each searcher a different colour, a reset marked by a vertical dotted blue line, and a completed search marked by a vertical dashed green line.
#'
#' @export



stochastic_search <- function(x0 = 40,
                              t = 1000,
                              searchers = 1,
                              reset_type = 'none',
                              reset_rate = 'optimal',  ## reciprocal of mean reset time
                              n = 1e5,
                              output_data = TRUE,
                              output_visual = TRUE){
  if (reset_rate=='optimal'){
    reset_rate <- 1/(2*(x0^2))
  }
  if (reset_type=='none'){
    if (output_visual==TRUE){
      multisearch_sbm(x0, t, n, searchers, output_data)
    } else {
      if (output_data==TRUE){
        multisearch_sbm_pro(x0, t, n, searchers)
      } else {
        return(NA)
      }
    }
  } else if (reset_type=='stochastic'){
    if (output_visual==TRUE){
      multisearch_SR(x0, t, n, reset_rate, searchers, output_data)
    } else {
      if (output_data==TRUE){
        multisearch_SR_pro(x0, t, n, reset_rate, searchers)
      } else {
        return(NA)
      }
    }
  } else if (reset_type=='deterministic'){
    if (output_visual==TRUE){
      multisearch_DR(x0, t, n, det_reset=(1/reset_rate), searchers, output_data)
    } else {
      if (output_data==TRUE){
        multisearch_DR_pro(x0, t, n, det_reset=(1/reset_rate), searchers)
      } else {
        return(NA)
      }
    }
  } else {
    return(NA)
  }
}


