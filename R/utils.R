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
                     title =  list(text ="Time", font = list(size = 30)),
                     tickfont = list(size = 25)),
                   yaxis = list(
                     zerolinewidth = 0,
                     gridcolor = '5555',
                     title = list(text ="ºC", font = list(size = 30)),
                     tickfont = list(size = 25)),
                   legend = list(font= list(size = 30)))

  fig
}

plot_two_cycles_plotly <- function(dt, id_var="profile_id", obj_var="pm", y1_id=12, y2_id=10, freq=0.5){
  y1 <- dt[get(id_var) == y1_id, get(obj_var)]
  x1 <- (1:length(y1))*freq
  y2 <- dt[get(id_var) == y2_id, get(obj_var)]
  x2 <- (1:length(y2))*freq

  y1_fig <- plot_ly(x = x1, y = y1, type = 'scatter', mode = 'lines', line=list(color='black'), name = paste0("Session ", y1_id))%>%
    layout(
      xaxis = list(
        zerolinecolor = '#ffff',
        zerolinewidth = 2,
        linecolor = 'black',
        gridcolor = '5555',
        title =  list(text ="Time", font = list(size = 30)),
        tickfont = list(size = 25)),
      yaxis = list(
        linecolor = 'black',
        zerolinecolor = '#ffff',
        zerolinewidth = 2,
        gridcolor = '5555',
        title = list(text ="ºC", font = list(size = 30)),
        tickfont = list(size = 25)),
      legend = list(font= list(size = 30)))

  y2_fig <- plot_ly(x = x2, y = y2, type = 'scatter', mode = 'lines', line=list(color='black'), name = paste0("Session ", y2_id))%>%
    layout(
      xaxis = list(
        zerolinecolor = '#ffff',
        zerolinewidth = 2,
        linecolor = 'black',
        gridcolor = '5555',
        title =  list(text ="Time", font = list(size = 30)),
        tickfont = list(size = 25)),
      yaxis = list(
        linecolor = 'black',
        zerolinecolor = '#ffff',
        zerolinewidth = 2,
        gridcolor = '5555',
        title = list(text ="ºC", font = list(size = 30)),
        tickfont = list(size = 25)),
      legend = list(font= list(size = 30)))


  fig <- subplot(y1_fig, y2_fig, nrows = 1)
  y1_fig
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
  params <- list(...)
  param_names <- as.list(substitute(list(...)))[-1L]
  sep <- "\n---------------------------------------\n\n"
  sink(res_file)
  cat("Beginning the experiment. The parameters are: \n")
  for(i in 2:length(params)){
    cat(paste0(param_names[i], ": ", params[[i]], "\n"))
  }
  cat(sep)
  sink()
}

print_current_results <- function(res_file, res_matrix, it){
  sep <- "\n---------------------------------------\n\n"
  sink(res_file, append = T)
  if(it == -1){ # Final results
    cat("The final results of the experiment are: \n")
    print(res_matrix)
  }
  else{
    cat(paste0("The results of fold number ", it, " are: \n"))
    print(res_matrix)
    cat(sep)
  }

  sink()
}

mae <- function(orig, pred){
  return(sum(abs(orig-pred)/length(orig)))
}

# Experiment pipeline functions

forecast_cycle_intervals_single <- function(f_dt_test, model_fit, id_var, size, obj_vars, prov_ev, pred_len){
  cycles <- f_dt_test[, unique(get(id_var))]
  reps <- f_dt_test[, ceiling(dim(.SD)[1] / pred_len), by=id_var]$V1
  res_matrix <- matrix(nrow = sum(reps), ncol = 2)
  global_rep <- 1

  for(i in 1:length(cycles)){
    ini <- 1
    din_pred_len <- pred_len
    for(j in 1:reps[i]){
      if(ini + din_pred_len > dim(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F])[1]) # Last iteration of each cycle is probably smaller
        din_pred_len <- dim(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F])[1] - ini + 1
      span <- Sys.time()
      res_cycle <- suppressWarnings(dbnR::forecast_ts(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F],
                                                      model_fit, size = size, obj_vars = obj_vars,
                                                      ini = ini, len = din_pred_len, prov_ev = prov_ev,
                                                      print_res = F, plot_res = F))
      res_matrix[global_rep, 2] <- span - Sys.time()
      res_matrix[global_rep, 1] <- mae(res_cycle$orig[,get(obj_vars)], res_cycle$pred[,get(obj_vars)])
      global_rep <- global_rep + 1
      ini <- ini + din_pred_len
    }
  }

  return(list(mean_res = apply(res_matrix, 2, mean), mae = res_matrix[,1]))
}

launch_single_model <- function(dt_train, f_dt_train, f_dt_test, id_var, obj_vars,
                                prov_ev, pred_len, size, method, min_ind, max_depth,
                                n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs){
  res_matrix <- matrix(nrow = 1, ncol = 3)
  span <- Sys.time() # Size were? TODO
  model_net <- dbnR::learn_dbn_struc(dt_train, size, method = method, f_dt = f_dt_train, n_it = n_it,
                                     n_ind = n_ind, gb_cte = gb_cte, lb_cte = lb_cte, r_probs = r_probs,
                                     v_probs = v_probs, cte = cte)
  model_fit <- dbnR::fit_dbn_params(model_net, f_dt_train)
  res_matrix[1,3] <- span - Sys.time()
  fore_results <- forecast_cycle_intervals_single(f_dt_test, model_fit, id_var, size, obj_vars, prov_ev, pred_len)
  res_matrix[1,1:2] <- fore_results$mean_res

  return(list(mean_res = res_matrix, mae = fore_results$mae))
}

forecast_cycle_intervals_hybrid <- function(f_dt_test, model, id_var, obj_vars, prov_ev, pred_len){
  cycles <- f_dt_test[, unique(get(id_var))]
  reps <- f_dt_test[, ceiling(dim(.SD)[1] / pred_len), by=id_var]$V1
  res_matrix <- matrix(nrow = sum(reps), ncol = 2)
  global_rep <- 1

  for(i in 1:length(cycles)){
    ini <- 1
    din_pred_len <- pred_len
    for(j in 1:reps[i]){
      if(ini + din_pred_len > dim(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F])[1]) # Last iteration of each cycle is probably smaller
        din_pred_len <- dim(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F])[1] - ini + 1
      span <- Sys.time()
      orig <- f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F]
      res_cycle <- suppressWarnings(model$forecast_ts(orig, obj_vars = obj_vars, ini = ini, len = din_pred_len,
                                                      prov_ev = prov_ev, print_res = F, plot_res = F, debug_m = F))
      res_matrix[global_rep, 2] <- span - Sys.time()
      res_matrix[global_rep, 1] <- mae(orig[ini:(ini+din_pred_len-1),get(obj_vars)], res_cycle[,get(obj_vars)])
      global_rep <- global_rep + 1
      ini <- ini + din_pred_len
    }
  }

  return(list(mean_res = apply(res_matrix, 2, mean), mae = res_matrix[,1]))
}

launch_hybrid_model <- function(dt_train, f_dt_train, f_dt_test, id_var,
                                obj_vars, obj_vars_tree, mv, homogen, prov_ev, pred_len, size, method,
                                min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte,
                                cte, r_probs, v_probs, prune_val){
  res_matrix <- matrix(nrow = 1, ncol = 3)
  train_obj_vars <- sapply(obj_vars_tree, function(x){strsplit(x, "_t_0")[[1]]}, USE.NAMES = F)
  span <- Sys.time()
  model <- mtDBN::mtDBN$new()
  model$fit_model(dt_train, size, method = method, obj_var = train_obj_vars, mv = mv, homogen = homogen,
                  min_ind = min_ind, max_depth = max_depth, f_dt = f_dt_train, n_it = n_it, n_ind = n_ind, gb_cte = gb_cte,
                  lb_cte = lb_cte, cte = cte, r_probs = r_probs, v_probs = v_probs, prune_val = prune_val)
  res_matrix[1,3] <- span - Sys.time()
  fore_results <- forecast_cycle_intervals_hybrid(f_dt_test, model, id_var, obj_vars, prov_ev, pred_len)
  res_matrix[1,1:2] <- fore_results$mean_res

  return(list(mean_res = res_matrix, mae = fore_results$mae))
}

train_test_iteration <- function(dt, id_var, test_id, obj_vars, obj_var_univ, obj_var_multiv, prov_ev, size = 3,
                                 method = "psoho", min_ind = 300, max_depth = 8, n_it = 100,
                                 n_ind = 100, gb_cte = 0.3, lb_cte = 0.7, cte = F,
                                 r_probs = c(-0.5, 1.5), v_probs = c(10,65,25), prune_val = 0.015, pred_len = 20){

  res_matrix <- matrix(nrow = 5, ncol = 3)
  colnames(res_matrix) <- c("MAE", "exec_time", "train_time")
  dt_train <- dt[!(get(id_var) %in% test_id)]
  dt_test <- dt[get(id_var) %in% test_id]

  f_dt_train <- dbnR::filtered_fold_dt(dt_train, size, id_var)
  dt_train[, eval(id_var) := NULL]
  f_dt_test <- dbnR::filtered_fold_dt(dt_test, size, id_var, clear_id_var = F)
  res_mae <- matrix(nrow = sum(f_dt_test[, ceiling(dim(.SD)[1] / pred_len), by=id_var]$V1), ncol = 5) # We have to calculate the total number of reps to know the rows

  res_tmp <- launch_single_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, prov_ev, pred_len, size, method, min_ind, max_depth,
                                 n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs)
  res_matrix[1,] <- res_tmp$mean_res
  res_mae[,1] <- res_tmp$mae

  res_tmp <- launch_hybrid_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, obj_var_univ, F, T, prov_ev, pred_len, size, method,
                                 min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
  res_matrix[2,] <- res_tmp$mean_res
  res_mae[,2] <- res_tmp$mae

  res_tmp <- launch_hybrid_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, obj_var_multiv, T, T, prov_ev, pred_len, size, method,
                                 min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
  res_matrix[3,] <- res_tmp$mean_res
  res_mae[,3] <- res_tmp$mae

  res_tmp <- launch_hybrid_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, obj_var_univ, F, F, prov_ev, pred_len, size, method,
                                 min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
  res_matrix[4,] <- res_tmp$mean_res
  res_mae[,4] <- res_tmp$mae

  res_tmp <- launch_hybrid_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, obj_var_multiv, T, F, prov_ev, pred_len, size, method,
                                 min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
  res_matrix[5,] <- res_tmp$mean_res
  res_mae[,5] <- res_tmp$mae

  return(list(mean_res = res_matrix, mae = res_mae))
}

full_exp_run <- function(dt, id_var, obj_vars, obj_var_univ, obj_var_multiv,
                         prov_ev, res_file, mae_file, pred_len, fold_len,
                         seed = NULL, size = 2, method = "psoho",
                         min_ind = 300, max_depth = 8, n_it = 100,
                         n_ind = 100, gb_cte = 0.3, lb_cte = 0.7, cte = F,
                         r_probs = c(-0.5, 1.5), v_probs = c(10,65,25),
                         prune_val = 0.015){

  res_matrix <- matrix(nrow = 5, ncol = 3, 0) # 5 different models, 3 columns: MAE, exec_time and train_time
  set.seed(seed)
  cv_sets <- cross_sets(dt[, unique(get(id_var))], fold_len)
  res_mae <- matrix(nrow = 0, ncol = 5) # I cannot know how many rows do I need without folding the test dataset and counting the number of needed repetitions, so I'll rbind
  colnames(res_mae) <- c("baseline", "m1", "m2", "m3", "m4")

  initialize_results_file(res_file, cv_sets, obj_vars, obj_var_univ, obj_var_multiv, prov_ev,
                          seed, size, method, min_ind, max_depth, n_it, n_ind, gb_cte,
                          lb_cte, cte, r_probs, v_probs, prune_val)

  for(i in 1:length(cv_sets)){
    message(paste0("Currently on the fold number ", i, " out of ", length(cv_sets)))
    res_tmp <- train_test_iteration(dt, id_var, cv_sets[[i]], obj_vars, obj_var_univ,
                                    obj_var_multiv, prov_ev, size, method, min_ind, max_depth,
                                    n_it, n_ind, gb_cte, lb_cte, cte, r_probs,
                                    v_probs, prune_val, pred_len)
    print_current_results(res_file, res_tmp$mean_res, i)

    res_matrix <- res_matrix + res_tmp$mean_res
    res_mae <- rbind(res_mae, res_tmp$mae)
  }

  fwrite(as.data.table(res_mae), file = mae_file)
  res_matrix <- res_matrix / length(cv_sets)
  print_current_results(res_file, res_matrix, -1)
}

# In case the total results are for some reason not saved in the results file, this function recovers
# the intermediate tables and aggregates them to obtain the final results
recover_results <- function(){
  res_matrix <- matrix(nrow = 5, ncol = 3, 0)
  n_hits <- 0
  regex <- "^[[:punct:][:digit:][:punct:][:punct:]]+ +[[:digit:]|[:punct:]]+ +[[:digit:]|[:punct:]]+ +[[:digit:]|[:punct:]]+$"

  for(line in readLines("full_run_results.txt")){
    if(length(grep(regex, line))){ # If it's a row of a intermediate matrix
      n_hits <- n_hits + 1
      if(n_hits %% 5 == 0)
        res_matrix[5,] <- res_matrix[5,] + as.numeric(strsplit(line," +")[[1]][2:4])
      else
        res_matrix[n_hits%%5,] <- res_matrix[n_hits%%5,] + as.numeric(strsplit(line," +")[[1]][2:4])
    }
  }

  res_matrix <- res_matrix / (n_hits / 5)
  print(res_matrix)

  return(res_matrix)
}
