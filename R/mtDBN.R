#' R6 class that defines the tree + DBN model with the rpart package
#'
#' The causal lists will be the base of the positions and the velocities
#' in the pso part of the algorithm.
#' @export
mtDBN <- R6::R6Class("mtDBN",
  public = list(
    #' @description
    #' Fit the model tree and each leaf's DBN
    #' @param dt_train a data.table with the training dataset
    #' @param size the size of the networks learned
    #' @param method the structure learning method used
    #' @param obj_var the objective variable for the model tree construction
    #' @param f_dt a previously folded dataset, in case some rows had to be deleted beforehand
    #' @param prune_val complexity parameter for the rpart prune function
    #' @param min_ind the minimum number of instances per leaf node
    #' @param inc the increment added to prune_val each time an invalid tree is generated
    #' @param max_depth maximum depth of the tree
    #' @param ... additional parameters for the structure learning
    #' @return A new 'causlist' object
    fit_model = function(dt_train, size, method, obj_var, f_dt = NULL, prune_val = 0.030, min_ind = 160, inc = 0.005, max_depth = 3, ...){
      # Security checks --ICO-Merge

      private$adjust_tree(dt_train, obj_var, prune_val, min_ind, inc, max_depth)
      private$fit_leaves(dt_train, size, method, f_dt, ...)

    },

    forecast_ts = function(f_dt, obj_vars, ini, len, prov_ev){
      # Security checks --ICO-Merge
      preds_test <- private$forecast_val_data_tree(f_dt, obj_vars, ini, len, prov_ev)
    },

    #' @description
    #' Export the tree model into a .ps file. To convert it into pdf, something
    #' like 'ps2pdf -dDEVICEWIDTHPOINTS=1479 -dDEVICEHEIGHTPOINTS=1598 rtree.ps rtree.pdf'
    #' in linux could be used.
    #' @param exp_dir the path in which to save the file. No path will put the file in the current directory
    #' @param width width of the exported image
    #' @param height height of the exported image
    #' @param paper size of the image. Set to special so that width and height can be used
    #' @param horizontal the orientation of the image. Horizontal if true, vertical otherwise.
    export_tree = function(exp_dir = NULL, width = 20, height = 20, paper = special, horizontal = FALSE){
      if(!is.null(exp_dir)){
        old_path <- getwd()
        setwd(exp_dir)
      }

      rpart::post(private$rtree, width, height, paper, horizontal)

      if(!is.null(exp_dir))
        setwd(old_path)
    }

  ),

  private = list(
    #' @field tree the tree structure
    rtree = NULL,
    #' @field models list of DBN models
    models = NULL,
    #' @field cl Total number of DBN models
    n_models = NULL,
    #' @field vars_t_0 names of the t_0 variables
    vars_t_0 = NULL,
    #' @field f_vars names of the t_0 variables
    f_vars = NULL,
    #' @field size the size of the networks learned
    size = NULL,

    formulate = function(obj, vars){
      res <- paste0(obj, " ~ ")
      res <- Reduce(function(acu, x){paste0(acu, " + ", x)}, vars[-1], ini = paste0(res, vars[1]))
      return(as.formula(res))
    },

    #' @description
    #' Counts the number of instances from the training dataset in each leaf node
    #' of the tree after classifying them.
    #' @param rt the rpart tree object used
    #' @param dt_train the training datset
    #' @return a vector with the number of instances in each class
    sum_unique = function(rt, dt_t_0){
      classif <- treeClust::rpart.predict.leaves(rt, dt_t_0, type = "where")
      u_classif <- unique(classif)
      s_classif <- rep(0, length(u_classif))

      for(i in classif)
        s_classif[u_classif == i] <- s_classif[u_classif == i] + 1

      return(s_classif)
    },

    #' @description
    #' Fits an rpart regression tree to the training dataset. It is then pruned
    #' and has to have a minimum number of instances per leaf node.
    #' @param dt_train the training dataset
    #' @param obj_var the response variable for the splits in the regression tree
    #' @param prune_val complexity parameter for the rpart prune function
    #' @param min_ind the minimum number of instances per leaf node
    #' @param inc the increment added to prune_val each time an invalid tree is generated
    adjust_tree = function(dt_train, obj_var, prune_val, min_ind, inc, max_depth){
      dt_t_0 <- dbnR::time_rename(dt_train)
      private$vars_t_0 <- copy(names(dt_t_0))
      pred_vars <- names(dt_t_0)
      obj_var <- paste0(obj_var, "_t_0")
      pred_vars <- pred_vars[!(pred_vars %in% obj_var)]

      valid <- F
      while (!valid) {
        res <- rpart::rpart(formula = private$formulate(obj_var, pred_vars),
                       data = dt_t_0, method = "anova", control = list(maxdepth = max_depth))
        res <- rpart::prune(res, cp = prune_val)

        valid <- min_ind <= min(private$sum_unique(res, dt_t_0))
        prune_val <- prune_val + inc
      }

      private$rtree <- res
    },

    #' @description
    #' Fits a DBN model to each leaf node of the rpart tree. Each model is fitted
    #' with the instances that correspond to the specific leaf node.
    #' @param dt_train the training dataset
    #' @param size the size of the networks learned
    #' @param method the structure learning method used
    #' @param f_dt a previously folded dataset, in case some rows had to be deleted beforehand
    #' @param ... additional parameters for the structure learning. An already folded dataset can be passed
    #' @import data.table
    fit_leaves = function(dt_train, size, method, f_dt = NULL, ...){
      network <- dbnR::learn_dbn_struc(dt_train, size, method, ...)
      if(is.null(f_dt))
        f_dt <- dbnR::fold_dt(dt_train, size)
      private$f_vars <- copy(names(f_dt))
      private$size <- size
      f_dt[, classif := treeClust::rpart.predict.leaves(private$rtree, f_dt[, .SD, .SDcols = private$vars_t_0], type = "where")] # Code local replacement

      private$n_models <- length(unique(f_dt$classif))
      private$models <- vector(mode = "list", private$n_models)
      names(private$models) <- unique(f_dt$classif) # The names of the list will be the internal number of the leaf nodes. Kind of a dictionary style access, but O(n)

      for(i in unique(f_dt$classif))
        private$models[[paste(i)]] <- dbnR::fit_dbn_params(network, f_dt[classif == i, .SD, .SDcols = private$f_vars])

      f_dt[, classif := NULL]

    },

    #' @description
    #' Classify a data.table with one row and select the appropriate model
    #' @param instance a data.table with one row
    #' @return the DBN model that corresponds to the instance
    get_model = function(instance){
      classif <- treeClust::rpart.predict.leaves(private$rtree, instance[, .SD, .SDcols = private$vars_t_0], type = "where")
      #print(classif)
      return(private$models[[paste(classif)]])
    },

    forecast_val_data_tree = function(f_dt, obj_vars, ini, len, prov_ev){
      res <- matrix(nrow = len, ncol = length(obj_vars), data = 0.0)
      colnames(res) <- obj_vars
      res <- data.table(res)
      instance <- f_dt[ini:(ini+1), .SD]

      # From the dbnR package, unify in a function --ICO-Merge
      var_names <- names(instance)
      vars_pred_idx <- grep("t_0", var_names)
      vars_subs_idx <- grep("t_1", var_names)
      vars_last_idx <- grep(paste0("t_", size-1), var_names)
      vars_pred <- var_names[vars_pred_idx]
      vars_prev <- var_names[-c(vars_pred_idx, vars_subs_idx)]
      vars_post <- var_names[-c(vars_pred_idx, vars_last_idx)]
      vars_ev <- var_names[-vars_pred_idx]
      vars_pred_crop <- vars_pred[!(vars_pred %in% prov_ev)]
      vars_subs_crop <- sub("t_0","t_1", vars_pred_crop)
      prov_ev_subs <- sub("t_0","t_1", prov_ev)

      for(i in 1:len){
        instance[2] <- f_dt[ini + i, .SD]

        # Classify with the tree and choose the model
        fit <- private$get_model(instance[1])

        # Forecast with len 1 with the correct model
        preds <- dbnR::forecast_ts(instance, fit, private$size, vars_pred_crop,
                                   print_res = F, plot_res = F, prov_ev = prov_ev,)

        # Move predictions into evidence
        if(length(vars_post) > 0)
          instance[1, (vars_prev) := .SD, .SDcols = vars_post]
        instance[1, (vars_subs_crop) := preds$pred[, .SD, .SDcols = vars_pred_crop]]
        if(!is.null(prov_ev)){
          instance[1, (prov_ev_subs) := .SD, .SDcols = prov_ev]
          instance[1, (prov_ev) := instance[2, .SD, .SDcols = prov_ev]]
        }

        res[i, (obj_vars) := preds$pred[, .SD, .SDcols = obj_vars]]
      }

      return(res)
    }
  ))
