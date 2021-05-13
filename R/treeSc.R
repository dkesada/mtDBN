#' R6 class that defines a tree scheme
#'
#' The causal lists will be the base of the positions and the velocities
#' in the pso part of the algorithm.
#' @export
treeSc <- R6::R6Class("treeSc",
    public = list(
      #' @description
      #' Create the tree scheme with a rpart tree
      #' @param tree the rpart tree that we want to get the structure from
      #' @return A new 'treeSc' object
      initialize = function(tree){
        cuts <- labels(tree) # The bloody cuts are in the labels of the object. It took me a while to figure this one out
        cuts_e <- parse(text = cuts)
        nodes <- as.numeric(row.names(tree$frame)) # The nodes are in preorder
        private$n_nodes <- length(nodes)
        depth <- private$tree_depth(nodes)
        e <- new.env()
        e$i <- 1
        private$tree_sc <- private$add_node_rec(e, nodes, depth, cuts_e, "(root)")
      },

      get_tree_sc = function(){
        return(private$tree_sc)
      },

      # Given a dataset, return a new data.table with the appropriate model for
      # each row. Used for dividing each instance for different dbn training.
      classify_dt = function(dt){
        idx <- dt[, .I]
        c_names <- sapply(1:private$n_nodes, function(x){paste0("node_",x)}) # Each node will have its own column
        dt[, eval(c_names) := 0]

        dt <- private$classify_dt_rec(dt, idx, private$tree_sc)
        res <- dt[, .SD, .SDcols = c_names]
        dt[, eval(c_names) := NULL]

        return(res)
      },

      #' @description
      #' Given an instance, return the appropriate node
      #' @param inst a row from a data.table
      #' @return
      classify_inst = function(inst){
        node_i <- private$tree_sc

        while(!is.null(node_i$r_node) && !is.null(node_i$l_node)){ # While not in a leaf
          if(nrow(inst[eval(node_i$l_node$cut_e)]))
            node_i <- node_i$l_node
          else
            node_i <- node_i$r_node
        }

        return(node_i)
      }
    ),

    private = list(
      #' @field tree_sc the tree structure scheme
      tree_sc = NULL,
      #' @field n_nodes the total number of nodes in the tree
      n_nodes = NULL,

      #' @description
      #' Returns the depth of each node in a tree.
      #' Unexported function from rpart. I duplicated it because I need to use it
      #' @param nodes the nodes in a tree named with natural depending on their order
      #' @return a vector with the depth of each node in the tree. 0 depth is the root
      tree_depth = function (nodes){
        depth <- floor(log(nodes, base = 2) + 1e-07)
        depth - min(depth)
      },

      #' @description
      #' Fit the model tree and each leaf's DBN
      #' @param e environment with the current index. It is an environment to allow passing by reference
      #' @param nodes names of the nodes in preorder
      #' @param depth depth of each node
      #' @param cuts_e expression of cut in each node
      #' @param parent environment of the parent node
      #' @return the defined subtree
      add_node_rec = function(e, nodes, depth, cuts_e, parent){
        curr_i <- e$i
        curr_node_e <- private$create_node(nodes[curr_i], parent, cuts_e[curr_i])
        e$i <- e$i + 1

        if(curr_i == length(nodes)) # No more nodes
          return(curr_node_e)

        else if(depth[curr_i] < depth[e$i]) # Next node is a child. Left branch
          curr_node_e$l_node <- private$add_node_rec(e, nodes, depth, cuts_e, curr_node_e)

        if(e$i > length(nodes)) # No more nodes after left branch
          return(curr_node_e)

        else if(depth[curr_i] < depth[e$i]) # Next node is a child. Right branch
          curr_node_e$r_node <- private$add_node_rec(e, nodes, depth, cuts_e, curr_node_e)

        return(curr_node_e)
      },

      #' @description
      #' Defines a node in the tree scheme as an environment with a node name,
      #' a left node, a right node, a parent node and a cut rule.
      #' @param name the integer defining the name of the node
      #' @param parent environment of the parent node
      #' @param cut_e expression of the cut in the node
      #' @return the new node environment
      create_node = function(name, parent, cut_e){
        res <- new.env() # I could create an R6 class instead, but I'm feeling lazy
        res$name <- name
        res$l_node <- NULL
        res$r_node <- NULL
        res$p_node <- parent
        res$cut_e <- cut_e

        return(res)
      },

      # Given a dataset and the corresponding node, recursively classify it
      classify_dt_rec = function(dt, idx, node){

        dt[idx, eval(paste0("node_", node$name)) := 1] # Always note the instances of a node

        if(!is.null(node$l_node)){ # If not in a root node
          idx_1 <- dt[idx, which(eval(node$l_node$cut_e))] # Apply the cuts. More readable than in 1 line
          dt <- private$classify_dt_rec(dt, idx_1, node$l_node)

          idx_2 <- dt[idx, which(eval(node$r_node$cut_e))]
          dt <- private$classify_dt_rec(dt, idx_2, node$r_node)
        }

        return(dt)
      }
      )

    )
