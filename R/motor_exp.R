# forecast_cycle_intervals_single <- function(f_dt_test, model, id_var, obj_vars, pred_len){
#   cycles <- f_dt_test[, unique(get(id_var))]
#   reps <- f_dt_test[, ceiling(dim(.SD)[1] / pred_len), by=id_var]$V1
#   res_matrix <- matrix(nrow = sum(reps), ncol = 2)
#   global_rep <- 1
#
#   for(i in 1:length(cycles)){
#     ini <- 1
#     din_pred_len <- pred_len
#     for(j in 1:reps[i]){
#       if(ini + din_pred_len > dim(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F])[1]) # Last iteration of each cycle is probably smaller
#         din_pred_len <- dim(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F])[1] - ini + 1
#       span <- Sys.time()
#       res_cycle <- suppressWarnings(dbnR::forecast_ts(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F],
#                                                       model_fit, size = size, obj_vars = obj_vars,
#                                                       ini = ini, len = din_pred_len, prov_ev = ev_vars,
#                                                       print_res = F, plot_res = F))
#       res_matrix[global_rep, 2] <- span - Sys.time()
#       res_matrix[global_rep, 1] <- mae(res_cycle$orig[,get(obj_vars)], res_cycle$pred[,get(obj_vars)])
#       global_rep <- global_rep + 1
#       ini <- ini + din_pred_len
#     }
#   }
#
#   return(apply(res_matrix, 2, mean))
# }
#
# launch_single_model <- function(f_dt_train, f_dt_test, id_var, obj_vars,
#                                 pred_len, size, method, min_ind, max_depth,
#                                 n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs){
#   res_matrix <- matrix(nrow = 1, ncol = 3)
#   span <- Sys.time() # Size were? TODO
#   model_net <- dbnR::learn_dbn_struc(dt_train, size, method = method, f_dt = f_dt_train, n_it = n_it,
#                                      n_ind = n_ind, gb_cte = gb_cte, lb_cte = lb_cte, r_probs = r_probs,
#                                      v_probs = v_probs, cte = cte)
#   model_fit <- dbnR::fit_dbn_params(model_net, f_dt_train)
#   res_matrix[1,3] <- span - Sys.time()
#   res_matrix[1,1:2] <- forecast_cycle_intervals_single(f_dt_test, model_fit, id_var, obj_vars, pred_len)
#
#   return(res_matrix)
# }
#
# forecast_cycle_intervals_hybrid <- function(f_dt_test, model, id_var, obj_vars, pred_len){
#   cycles <- f_dt_test[, unique(get(id_var))]
#   reps <- f_dt_test[, ceiling(dim(.SD)[1] / pred_len), by=id_var]$V1
#   res_matrix <- matrix(nrow = sum(reps), ncol = 2)
#   global_rep <- 1
#
#   for(i in 1:length(cycles)){
#     ini <- 1
#     din_pred_len <- pred_len
#     for(j in 1:reps[i]){
#       if(ini + din_pred_len > dim(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F])[1]) # Last iteration of each cycle is probably smaller
#         din_pred_len <- dim(f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F])[1] - ini + 1
#       span <- Sys.time()
#       orig <- f_dt_test[get(id_var) == cycles[i], !eval(id_var), with=F]
#       res_cycle <- suppressWarnings(model$forecast_ts(orig, obj_vars = obj_vars, ini = ini, len = din_pred_len,
#                                                       prov_ev = ev_vars, print_res = F, plot_res = F, debug_m = F))
#       res_matrix[global_rep, 2] <- span - Sys.time()
#       res_matrix[global_rep, 1] <- mae(orig[ini:(ini+din_pred_len-1),get(obj_vars)], res_cycle[,get(obj_vars)])
#       global_rep <- global_rep + 1
#       ini <- ini + din_pred_len
#     }
#   }
#
#   return(apply(res_matrix, 2, mean))
# }
#
# launch_hybrid_model <- function(dt_train, f_dt_train, f_dt_test, id_var,
#                                 obj_vars, mv, homogen, pred_len, size, method,
#                                 min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte,
#                                 cte, r_probs, v_probs, prune_val){
#   res_matrix <- matrix(nrow = 1, ncol = 3)
#   train_obj_vars <- sapply(obj_vars, function(x){strsplit(x, "_t_0")[[1]]}, USE.NAMES = F)
#   span <- Sys.time()
#   model <- mtDBN::mtDBN$new()
#   model$fit_model(dt_train, size, method = method, obj_var = train_obj_vars, mv = mv, homogen = homogen,
#                   min_ind = min_ind, max_depth = max_depth, f_dt = f_dt_train, n_it = n_it, n_ind = n_ind, gb_cte = gb_cte,
#                   lb_cte = lb_cte, cte = cte, r_probs = r_probs, v_probs = v_probs, prune_val = prune_val)
#   res_matrix[1,3] <- span - Sys.time()
#   res_matrix[1,1:2] <- forecast_cycle_intervals_hybrid(f_dt_test, model, id_var, obj_vars, pred_len)
#
#   return(res_matrix)
# }
#
# train_test_iteration <- function(dt, id_var, test_id, obj_vars, size = 3,
#                                  method = "psoho", min_ind = 300, max_depth = 8, n_it = 100,
#                                  n_ind = 100, gb_cte = 0.3, lb_cte = 0.7, cte = F,
#                                  r_probs = c(-0.5, 1.5), v_probs = c(10,65,25), prune_val = 0.015, pred_len = 20){
#
#   res_matrix <- matrix(nrow = 5, ncol = 3)
#   colnames(res_matrix) <- c("MAE", "exec_time", "train_time")
#   pred_vars <- names(dt_train)[!(names(dt_train) %in% obj_vars)]
#   dt_train <- dt[!(profile_id %in% test_id)]
#   dt_test <- dt[profile_id %in% test_id]
#
#   f_dt_train <- dbnR::filtered_fold_dt(dt_train, size, id_var)
#   dt_train[, eval(id_var) := NULL]
#   f_dt_test <- dbnR::filtered_fold_dt(dt_test, size, id_var, clear_id_var = F)
#
#   res_matrix[1,] <- launch_single_model(f_dt_train, f_dt_test, id_var, obj_vars, pred_len, size, method, min_ind, max_depth,
#                                         n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs)
#   res_matrix[2,] <- launch_hybrid_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, F, T, pred_len, size, method,
#                                         min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
#   res_matrix[3,] <- launch_hybrid_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, T, T, pred_len, size, method,
#                                         min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
#   res_matrix[4,] <- launch_hybrid_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, T, F, pred_len, size, method,
#                                         min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
#   res_matrix[5,] <- launch_hybrid_model(dt_train, f_dt_train, f_dt_test, id_var, obj_vars, F, F, pred_len, size, method,
#                                         min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val)
#
#   return(res_matrix)
# }
#
#
#
# full_exp_motor_run <- function(dt, id_var, obj_vars, seed = NULL, size = 3,
#                                method = "psoho", min_ind = 300, max_depth = 8, n_it = 100,
#                                n_ind = 100, gb_cte = 0.3, lb_cte = 0.7, cte = F,
#                                r_probs = c(-0.5, 1.5), v_probs = c(10,65,25), prune_val = 0.015, pred_len = 20,
#                                fold_len = 3, res_file = "full_run_motor_results.txt"){
#
#   res_matrix <- matrix(nrow = 5, ncol = 3, 0) # 5 different models, 3 columns: MAE, exec_time and train_time
#   set.seed(seed)
#   cv_sets <- cross_sets(dt[, unique(get(id_var))], fold_len)
#
#   initialize_results_file(res_file, cv_sets, obj_vars, seed,
#                           method, min_ind, max_depth, n_it, n_ind, gb_cte, lb_cte, cte,
#                           r_probs, v_probs, prune_val)
#
#   for(i in 1:length(cv_sets)){
#     message(paste0("Currently on the fold number ", i, " out of ", length(cv_sets)))
#     res_tmp <- train_test_iteration(dt, id_var, cv_sets[[i]], obj_vars, size, method, min_ind, max_depth,
#                                     n_it, n_ind, gb_cte, lb_cte, cte, r_probs, v_probs, prune_val, pred_len)
#     print_current_results(res_file, res_tmp, i)
#     res_matrix <- res_matrix + res_tmp
#   }
#
#   res_matrix <- res_matrix / length(cv_sets)
#   print_current_results(res_file, res_matrix, -1)
# }
#
# main_prep_and_run_motor <- function(){
#   dt <- data.table::fread("./dataset/motor_measures_v2.csv") # 0.5 secs between rows
#   id_var <- "profile_id"
#
#   dt[, torque := NULL]
#   dt <- dbnR::reduce_freq(dt, 30, 0.5, id_var) # 30 secs between rows
#   obj_vars <- "pm_t_0"
#   full_exp_motor_run(dt, id_var, obj_vars, seed = 42, size = 3, n_it = 100)
# }
#
# # In case the total results are for some reason not saved in the results file, this function recovers
# # the intermediate tables and aggregates them to obtain the final results
# recover_results <- function(){
#   res_matrix <- matrix(nrow = 5, ncol = 3, 0)
#   n_hits <- 0
#   regex <- "^[[:punct:][:digit:][:punct:][:punct:]]+ +[[:digit:]|[:punct:]]+ +[[:digit:]|[:punct:]]+ +[[:digit:]|[:punct:]]+$"
#
#   for(line in readLines("full_run_results.txt")){
#     if(length(grep(regex, line))){ # If it's a row of a intermediate matrix
#       n_hits <- n_hits + 1
#       if(n_hits %% 5 == 0)
#         res_matrix[5,] <- res_matrix[5,] + as.numeric(strsplit(line," +")[[1]][2:4])
#       else
#         res_matrix[n_hits%%5,] <- res_matrix[n_hits%%5,] + as.numeric(strsplit(line," +")[[1]][2:4])
#     }
#   }
#
#   res_matrix <- res_matrix / (n_hits / 5)
#   print(res_matrix)
#
#   return(res_matrix)
# }
