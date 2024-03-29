% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtDBN.R
\name{mtDBN}
\alias{mtDBN}
\title{R6 class that defines the tree + DBN model with the rpart package}
\description{
R6 class that defines the tree + DBN model with the rpart package

R6 class that defines the tree + DBN model with the rpart package
}
\details{
The model is defined as a tree with several DBNs learned for each node.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{tree_sc}}{the tree structure scheme}

\item{\code{models}}{list of DBN models}

\item{\code{cl}}{Total number of DBN models}

\item{\code{f_vars}}{names of the t_0 variables}

\item{\code{size}}{the size of the networks learned}

\item{\code{homogen}}{whether the DBN structure is the same in all leaves or not}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{tree_sc}}{the tree structure scheme}

\item{\code{models}}{list of DBN models}

\item{\code{cl}}{Total number of DBN models}

\item{\code{f_vars}}{names of the t_0 variables}

\item{\code{size}}{the size of the networks learned}

\item{\code{homogen}}{whether the DBN structure is the same in all leaves or not}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-fit_model}{\code{mtDBN$fit_model()}}
\item \href{#method-forecast_ts}{\code{mtDBN$forecast_ts()}}
\item \href{#method-export_tree}{\code{mtDBN$export_tree()}}
\item \href{#method-get_models}{\code{mtDBN$get_models()}}
\item \href{#method-get_tree}{\code{mtDBN$get_tree()}}
\item \href{#method-clone}{\code{mtDBN$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit_model"></a>}}
\if{latex}{\out{\hypertarget{method-fit_model}{}}}
\subsection{Method \code{fit_model()}}{
Fit the model tree and each leaf's DBN
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mtDBN$fit_model(
  dt_train,
  size,
  method,
  obj_var,
  mv = FALSE,
  f_dt = NULL,
  homogen = TRUE,
  prune_val = 0.03,
  min_ind = 160,
  inc = 0.005,
  max_depth = 3,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dt_train}}{a data.table with the training dataset}

\item{\code{size}}{the size of the networks learned}

\item{\code{method}}{the structure learning method used}

\item{\code{obj_var}}{the objective variable for the model tree construction}

\item{\code{mv}}{if TRUE, a multivariate tree will be made. If FALSE, it will be univariate}

\item{\code{f_dt}}{a previously folded dataset, in case some rows had to be deleted beforehand}

\item{\code{prune_val}}{complexity parameter for the rpart prune function}

\item{\code{min_ind}}{the minimum number of instances per leaf node}

\item{\code{inc}}{the increment added to prune_val each time an invalid tree is generated}

\item{\code{max_depth}}{maximum depth of the tree}

\item{\code{...}}{additional parameters for the structure learning}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new 'causlist' object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-forecast_ts"></a>}}
\if{latex}{\out{\hypertarget{method-forecast_ts}{}}}
\subsection{Method \code{forecast_ts()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mtDBN$forecast_ts(
  f_dt,
  obj_vars,
  ini,
  len,
  prov_ev,
  print_res = TRUE,
  plot_res = TRUE,
  debug_m = TRUE
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-export_tree"></a>}}
\if{latex}{\out{\hypertarget{method-export_tree}{}}}
\subsection{Method \code{export_tree()}}{
Export the tree model into a .ps file. To convert it into pdf, something
like 'ps2pdf -dDEVICEWIDTHPOINTS=1479 -dDEVICEHEIGHTPOINTS=1598 rtree.ps rtree.pdf'
in linux could be used.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mtDBN$export_tree(
  exp_dir = NULL,
  width = 20,
  height = 20,
  paper = "special",
  horizontal = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{exp_dir}}{the path in which to save the file. No path will put the file in the current directory}

\item{\code{width}}{width of the exported image}

\item{\code{height}}{height of the exported image}

\item{\code{paper}}{size of the image. Set to special so that width and height can be used}

\item{\code{horizontal}}{the orientation of the image. Horizontal if true, vertical otherwise.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_models"></a>}}
\if{latex}{\out{\hypertarget{method-get_models}{}}}
\subsection{Method \code{get_models()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mtDBN$get_models()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_tree"></a>}}
\if{latex}{\out{\hypertarget{method-get_tree}{}}}
\subsection{Method \code{get_tree()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mtDBN$get_tree()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{mtDBN$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
