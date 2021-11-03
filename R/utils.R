# Utils functions for the physical system of ODEs model

# ICO-Merge fix the grammatical errors and the incomplete descriptions

# Some functions are prepared to be used with the 'deSolve' package or with the
# 'ode_var_params' wrapper. This functions will be identifiable by having
# always the parameters 'time', 'y' and 'params' and begining with an
# 'with(as.list(c(y, params)), {...}'. This functions define the equations of
# the model and should respect the 'deSolve' syntax.

# Simple sigmoid function.
sigmoid <- function(t){
  1 / (1 + exp(1)^-t)
}


coupled_Sc_Ca_t1_t2_t3_vol <- function(time, y, params){
  with(as.list(c(y, params)), {
    T_3 <- t_min + sigmoid(m_c) * (t_max - t_min)
    f_2 <-  (A_1 / y[1]) * (y[3] - y[4]) + eps_2 * (T_3^4 - y[4]^4)
    dT_2 <- f_2 / (rho_2 * C_p2)
    #Q_in <- (pi * delta_p * r^4) / (8 * mu_1 * l)
    Q_in <- (pi * vol * (2*r)^2) / (4)
    delta_T1 <- A_6 * Q_in * (Tin - y[3])
    f_1 <-  (A_1 / y[1]) * (y[4] - y[3])
    dT_1 <- (f_1 + rho_1 * C_p1 * C_in * Tin - 2 * rho_1 * C_p1 * C_out * y[3]) / (vol * rho_1 * C_p1) #+ delta_T1
    k_1 <- A_3 * exp(-A_5 / (R * y[3]))
    k_2 <- A_3_p * exp(-A_5_p / (R * y[3]))
    dS_c <- A_2 * k_1 * y[2]
    dC_a <- -A_4 * k_2 * y[2] + (C_in/vol) * C_ain - (C_out/vol) * y[2]
    list(c(dS_c, dC_a, dT_1, dT_2))
  })
}

# Security check to see if the parameters are of either length 1 or length
# equal to times.
initial_params_check <- function(params, times){
  invisible(lapply(params, function(x){
    if(length(x) != 1 && length(x) != length(times))
      stop("The length of the parameters has to be either 1 or the length of 'times'")
  }))
}

# Function that gets the parameters of the i iteration.
get_params <- function(params_m, i){
  if(nrow(params_m) == 1)
    res <- params_m[1,]
  else
    res <- params_m[i,]

  return(as.list(res))
}

# This function allows for a list with vectors of parameters, so that
# interventions or varying parameters over time can be used
#
# y = initial values of the modelled variables
# times = sequence of the instants to be modelled
# foo = the function that defines the equations of the model
# params = list with the values of the parameters, one value if constant all
#          the time, a vector of equal length to times ioc
# mu_1 = mean of the Gaussian noise added to T_1. Both mean and sd set to 0 equal no noise.
# mu_2 = mean of the Gaussian noise added to T_2
# sig_1 = standard deviation of the Gaussian noise added to T_1
# sig_2 = standard deviation of the Gaussian noise added to T_2
# returns a matrix with the value of each variable in each instant of time
ode_var_params <- function(y, times, foo, params, mu_1 = 0, mu_2 = 0, sigma_1 = 0.2, sigma_2 = 5){
  # Security checks
  initial_params_check(params, times)

  res <- matrix(nrow = length(times), ncol = 1 + length(y), data = 0.0)
  colnames(res) <- c("time", names(y))
  params_m <- do.call(cbind, params)
  times_it <- c(1,2)
  y_it <- y
  for(i in 1:length(times)){
    res_it <- deSolve::ode(y = y_it, times = times_it, func = foo, parms = get_params(params_m, i))
    res_it[,4] <- res_it[,4] + rnorm(1, mu_1, sigma_1) # White noise for T_1
    res_it[,5] <- res_it[,5] + rnorm(1, mu_2, sigma_2) # White noise for T_2
    res[i,] <- res_it[2,]
    y_it <- as.vector(res_it[2,-1])
    names(y_it) <- colnames(res_it[-1])
  }
  res[,1] <- times
  class(res) <- c("deSolve", "matrix")

  return(res)
}

plot_model_change <- function(orig, pred, model_ch){
  mask <- c(F, sapply(2:length(model_ch), function(x){model_ch[x-1] != model_ch[x]}))
  mask[!mask] <- mask[!mask] * NA
  mask <- pred * mask
  plot(ts(orig))
  lines(pred, col = "blue")
  lines(mask, col = "red", type = "dot")
}

plot_model_change_plotly <- function(orig, pred, model_ch){
  x_line <- 1:length(orig)
  mask <- c(F, sapply(2:length(model_ch), function(x){model_ch[x-1] != model_ch[x]}))
  mask[!mask] <- mask[!mask] * NA
  mask <- pred * mask
  fig <- plot_ly(x = x_line, y = orig, type = 'scatter', mode = 'lines+markers', name = "Fluid \ntemperature", 
                 line = list(color = 'rgba(0, 0, 0, .9)'), marker = list(color = 'rgba(0, 0, 255, 0)'))
  fig <- fig %>% add_trace(y = pred,  mode = 'lines+markers', name = 'Prediction', 
                           line = list(color = 'rgba(255, 0, 0, .9)'), marker = list(size = 0)) 
  fig <- fig %>% add_trace(y = mask, type = 'scatter', mode = 'lines+markers', name = 'Model \nchange',
                           marker = list(size = 11, color = 'rgba(0, 0, 255, .9)'), line = list(color = 'rgba(0, 0, 255, 0)'))
  fig <- fig %>% layout(
                   xaxis = list(
                     zerolinewidth = 2,
                     linecolor = 'black',
                     gridcolor = '5555',
                     title =  list(text ="Time", font = list(size = 25))),
                   yaxis = list(
                     zerolinewidth = 0,
                     gridcolor = '5555',
                     title = list(text ="ºC", font = list(size = 25))),
                   legend = list(font= list(size = 25)))
  
  fig
}

# Create a combination of m disjoint sets of k elements each from the original set n
# It will only create 1 combination, not all the existing ones like a leave-3-out cross-validation
cross_sets <- function(n, k){
  res <- vector(mode = "list", length = floor(length(n) / k))
  
  for(i in 1:length(res)){
    ids <- sample.int(length(n), size = k, replace = F)
    res[[i]] <- n[ids]
    n <- n[-ids]
  }
  
  return(res)
}

# Create the results file by sinking an initial line and the values of the parameters in the elipsis argument
initialize_results_file <- function(res_file, ...){
  params <- as.list(substitute(list(...)))
  sep <- "\n---------------------------------------\n\n"
  sink(res_file)
  cat("Beginning the experiment. The parameters are: \n")
  for(i in 2:length(params)){
    cat(paste0(names(params)[i], ": ", params[i], "\n"))
  }
  cat(sep)
  sink(NULL)
}

mae <- function(orig, pred){
  return(sum(abs(orig-pred)/length(orig)))
}

forecast_cycle_intervals <- function(f_dt_test, model, id_var, ...){
  cycles <- f_dt_test[, unique(get(id_var))]
  reps <- f_dt_test[, ceiling(dim(.SD)[1] / pred_len), by=id_var]$V1
  res_mae <- vector(mode = "numeric", length=sum(reps))
  for(i in cycles){
    ini <- 1
    # TODO
  }
}

launch_single_model <- function(f_dt_train, f_dt_test, obj_vars, method, min_ind, max_depth,
                                n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs){
  cat("Single model training:\n")
  tmp <- Sys.time()
  model_net <- dbnR::learn_dbn_struc(dt_train, size, method = method, f_dt = f_dt_train, n_it = n_it,
                                     n_ind = n_ind, gb_cte = gb_cte, lb_cte = lb_cte, r_probs = r_probs,
                                     v_probs = v_probs, cte = cte)
  model_fit <- dbnR::fit_dbn_params(model_net, f_dt_train)
  
  cat(paste0("Elapsed time: ", tmp - Sys.time()))
  
  print("Forecasting time for single net: ")
  res_net <- suppressWarnings(dbnR::forecast_ts(f_dt_test, model_fit, size = size, 
                                                obj_vars = obj_vars, ini = ini, len = len, prov_ev = ev_vars))
  
  cat("\n---------------------------------------\n\n")
}

launch_hybrid_model <- function(){
  
}

train_test_iteration <- function(dt, id_var, test_id, obj_vars, 
                                 method = "psoho", min_ind = 300, max_depth = 8, n_it = 100,
                                 n_ind = 100, gb_cte = 0.3, lb_cte = 0.7, cte = F, 
                                 r_probs = c(-0.5, 1.5), v_probs = c(10,65,25), prune_val = 0.015){
  
  res_mae <- vector(mode = "numeric", length = 5)
  pred_vars <- names(dt_train)[!(names(dt_train) %in% obj_vars)]
  dt_train <- dt[!(profile_id %in% test_id)]
  dt_test <- dt[profile_id %in% test_id]
  
  f_dt_train <- dbnR::filtered_fold_dt(dt_train, size, id_var)
  f_dt_test <- dbnR::filtered_fold_dt(dt_test, size, id_var, clear_id_var = F)
  dt_train[, (id_var) := NULL]
  #dt_test[, (id_var) := NULL]
  
  res_mae[1] <- launch_single_model()
  res_mae[2] <- launch_hybrid_model()
  res_mae[3] <- launch_hybrid_model()
  res_mae[4] <- launch_hybrid_model()
  res_mae[5] <- launch_hybrid_model()
  
  return(res_mae)
}

full_exp_motor_run <- function(dt, id_var, cross_sets, obj_vars, seed = NULL,
                               method = "psoho", min_ind = 300, max_depth = 8, n_it = 100,
                               n_ind = 100, gb_cte = 0.3, lb_cte = 0.7, cte = F, 
                               r_probs = c(-0.5, 1.5), v_probs = c(10,65,25), prune_val = 0.015,
                               res_file = "full_run_results.txt"){
  
  
  res_mae <- matrix(nrow = length(cross_sets), ncol = 5)
  set.seed(seed)
  
  initialize_results_file(res_file = "full_run_results.txt", cross_sets, obj_vars, seed,
                          method, min_ind, max_depth, n_it, n_ind, gb_cte, lb_ct, cte, 
                          r_probs, v_probs, prune_val)
  
  for(i in 1:length(cross_sets)){
    res_mae[i,] <- train_test_iteration(dt, id_var, cross_sets[[i]], obj_vars, seed, method, min_ind, max_depth, 
                                        n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
    
  }
  
}

