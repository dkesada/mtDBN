#' @export
main_prep_and_run_synth <- function(){
  dt <- data.table::fread("./dataset/dt_cycles.csv") # 0.5 secs between rows
  vars_noise <- c("rho_1", "C_p1", "C_in", "vol", "C_ain") # Add some noise to avoid singular matrix errors later
  set.seed(42)
  dt[, eval(vars_noise) := .SD + rnorm(100000, 0, 0.5), .SDcols = vars_noise]
  set.seed(NULL)
  id_var <- "cyc"

  obj_vars <- "T_1_t_0"
  obj_var_univ <- "T_1"
  obj_var_multiv <- c("T_1","T_2", "S_c", "C_a", "C_p1", "rho_1")

  full_exp_run(dt = dt, id_var = id_var, obj_vars = obj_vars,
                     obj_var_univ = obj_var_univ, obj_var_multiv = obj_var_multiv,
                     res_file = "full_run_synth_results.txt",
                     mae_file = "full_run_synth_mae.csv", pred_len = 99,
                     fold_len = 33, seed = 42, size = 2, n_it = 1)
}
