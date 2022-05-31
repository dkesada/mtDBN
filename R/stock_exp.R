#' @export
main_prep_and_run_stock <- function(){
  dt <- data.table::fread("./dataset/TWII_1y.csv") # 1 day between rows
  id_var <- "Date" # We will crossvalidate based on the months
  
  dt[, Date := unlist(strsplit(as.character(dt$Date), "-[0-9]*$"))] # We leave only the year and month of the dates for making CV groups
  dt[, Volume := as.numeric(Volume)]
  setnames(dt, "Adj Close", "Adj_Close")
  obj_vars <- "Open_t_0"
  obj_var_univ <- "Open"
  obj_var_multiv <- c("Open", "Volume")
  prov_ev <- NULL
  # max_min_vals <- max_min_norm(dt, names(dt)[-1])
  dt[, Adj_Close := NULL] # Exactly the same as Close. Correlations of 1 break the inference process
    
  full_exp_run(dt = dt, id_var = id_var, obj_vars = obj_vars,
               obj_var_univ = obj_var_univ, obj_var_multiv = obj_var_multiv,
               prov_ev = prov_ev, res_file = "full_run_stock_results.txt",
               mae_file = "full_run_stock_mae.csv", mape_file = "full_run_stock_mape.csv", pred_len = 1,
               fold_len = 2, seed = 42, size = 3, method = "psoho", n_it = 100, min_ind = 50)
}