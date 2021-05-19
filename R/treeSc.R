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
        private$rtree <- tree
        cuts <- labels(tree) # The bloody cuts are in the labels of the object. It took me a while to figure this one out
        cuts_e <- parse(text = cuts)
        private$node_names <- row.names(tree$frame)
        nodes <- as.numeric(private$node_names) # The nodes are in preorder
        private$n_nodes <- length(nodes)
        depth <- private$tree_depth(nodes)
        e <- new.env()
        e$i <- 1
        private$tree_sc <- private$add_node_rec(e, nodes, depth, cuts_e, "(root)")
      },

      get_tree_sc = function(){
        return(private$tree_sc)
      },

      #' @description
      #' Given a dataset, return a new data.table with the appropriate model for
      #' each row. Used for dividing each instance for different dbn training.
      #' @param dt the dataset to be classified
      #' @return a dataset with the nodes that each instance falls into
      classify_dt = function(dt){
        idx <- dt[, .I]
        dt[, eval(private$node_names) := 0] # Each node will have its own column

        dt <- private$classify_dt_rec(dt, idx, private$tree_sc)
        res <- dt[, .SD, .SDcols = private$node_names]
        dt[, eval(private$node_names) := NULL]

        return(res)
      },

      #' @description
      #' Given an instance, return the appropriate node
      #' @param inst a row from a data.table
      #' @return the environment of the corresponding node
      classify_inst = function(inst){
        node_i <- private$tree_sc

        while(!is.null(node_i$r_node) && !is.null(node_i$l_node)){ # While not in a leaf
          if(nrow(inst[eval(node_i$l_node$cut_e)]))
            node_i <- node_i$l_node
          else
            node_i <- node_i$r_node
        }

        return(node_i)
      },

      #' @description
      #' Getter of the original rpart tree
      #' --ICO-Merge: delete if it remains unused
      #' @return the original rpart object
      get_rtree = function(){
        return(private$rtree)
      }
    ),

    private = list(
      #' @field tree_rpart the original rpart object
      rtree = NULL,
      #' @field tree_sc the tree structure scheme
      tree_sc = NULL,
      #' @field n_nodes the total number of nodes in the tree
      n_nodes = NULL,
      #' @field node_names the names of the nodes in the tree
      node_names = NULL,

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

      #' @description
      #' Given a dataset and the corresponding node, recursively classify it.
      #' This is the private and recursive part of the function.
      #' @param dt the dataset that is being classified
      #' @param idx the specific rows of the dataset that are being processed
      #' @param node the environment of the tree node that the algorithm is traversing
      #' @return the classified dataset of the whole subtree
      classify_dt_rec = function(dt, idx, node){

        dt[idx, eval(as.character(node$name)) := 1] # Always note the instances of a node

        if(!is.null(node$l_node)){ # If not in a root node
          idx_1 <- dt[idx, which(eval(node$l_node$cut_e))] # Apply the cuts. More readable than in 1 line
          idx_1 <- idx[idx_1] # Notice that the previous line returns indexes from 1 to len(idx) due to the 'which'
          dt <- private$classify_dt_rec(dt, idx_1, node$l_node)

          idx_2 <- idx[dt[idx, which(eval(node$r_node$cut_e))]]
          dt <- private$classify_dt_rec(dt, idx_2, node$r_node)
        }

        return(dt)
      }
      )

    )
