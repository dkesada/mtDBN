#' @export
main_prep_and_run_motor <- function(){
  dt <- data.table::fread("./dataset/motor_measures_v2.csv") # 0.5 secs between rows
  id_var <- "profile_id"

  dt[, torque := NULL]
  dt <- dbnR::reduce_freq(dt, 30, 0.5, id_var) # 30 secs between rows
  obj_vars <- "pm_t_0"
  obj_var_univ <- "pm"
  obj_var_multiv <- c("pm","stator_tooth", "stator_winding", "stator_yoke")

  full_exp_run(dt = dt, id_var = id_var, obj_vars = obj_vars,
                     obj_var_univ = obj_var_univ, obj_var_multiv = obj_var_multiv,
                     res_file = "full_run_motor_results.txt",
                     mae_file = "full_run_synth_mae.csv", pred_len = 20,
                     fold_len = 3, seed = 42, size = 3, n_it = 1)
}


