#' Normal To Anything (NorTA) core
#'
#' Generate correlated samples with arbitrary marginal distributions
#' using the NorTA approach.
#'
#' @param n Integer, number of observations (samples).
#' @param cor Correlation structure: either a scalar correlation
#'   (for a 2x2 case) or a correlation matrix (D x D).
#' @param qdist Quantile function for the target marginal distribution.
#'   Must have an argument named \code{p} for the probabilities.
#' @param param Optional named vector or matrix/data.frame of parameters
#'   for \code{qdist}. If a vector, parameters are shared across all
#'   dimensions. If a matrix/data.frame, rows correspond to dimensions
#'   and columns to parameter names (which must match arguments of \code{qdist}).
#' @param check_pd (default TRUE) Logical, whether to check that \code{cor} is positive-definite.
#' @param tol (default 1e-8) Numerical tolerance for positive-definiteness.
#'
#' @return A numeric matrix of size \code{n x D}.
#'
#' @export
norta <- function(
    n,
    cor,
    qdist,
    param    = NULL,
    check_pd = TRUE,
    tol      = 1e-8
) {
  # --- basic checks ---------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1L || n <= 0 || n != round(n)) {
    cli::cli_abort("`n` must be a positive integer scalar.")
  }
  
  if (!is.function(qdist)) {
    cli::cli_abort("`qdist` must be a function (a quantile function).")
  }
  
  qformals <- names(formals(qdist))
  if (!("p" %in% qformals)) {
    cli::cli_abort(
      "`qdist` must have an argument named {.arg p} for the probabilities."
    )
  }
  
  # --- handle correlation input --------------------------------------------
  if (length(cor) == 0L) {
    cli::cli_abort("`cor` has length 0. Please provide a scalar or a matrix.")
  }
  
  if (is.null(dim(cor))) {
    # scalar -> 2x2 correlation matrix
    if (!is.numeric(cor) || length(cor) != 1L || abs(cor) > 1) {
      cli::cli_abort(
        "`cor` as scalar must be a single numeric value in [-1, 1]."
      )
    }
    cor <- matrix(c(1, cor, cor, 1), nrow = 2L)
  } else {
    if (!is.numeric(cor)) {
      cli::cli_abort("`cor` matrix must be numeric.")
    }
    if (nrow(cor) != ncol(cor)) {
      cli::cli_abort("`cor` must be a square matrix.")
    }
    if (!isSymmetric(unname(cor), tol = 1e-10)) {
      cli::cli_abort("`cor` must be symmetric.")
    }
    if (any(abs(cor) > 1 + 1e-10)) {
      cli::cli_abort("All entries of `cor` must lie in [-1, 1].")
    }
    if (any(abs(diag(cor) - 1) > 1e-10)) {
      cli::cli_abort("Diagonal entries of `cor` must all be equal to 1.")
    }
  }
  
  D <- ncol(cor)
  
  if (is.null(colnames(cor))) {
    colnames(cor) <- paste0("V", seq_len(D))
    cli::cli_inform(
      "Column names of `cor` were missing; default names {colnames(cor)} have been assigned."
    )
  }
  
  # --- positive-definiteness check -----------------------------------------
  if (check_pd) {
    ev <- eigen(cor, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < tol) {
      cli::cli_abort(c(
        "The correlation matrix `cor` is not positive-definite.",
        "i" = "Minimum eigenvalue is {signif(min(ev), 3)}, which is < tol = {tol}.",
        "i" = "Consider adjusting or regularizing `cor` before calling `norta()`."
      ))
    }
  }
  
  # --- handle param --------------------------------------------------------
  if (!is.null(param)) {
    
    # data.frame -> matrix
    if (is.data.frame(param)) {
      param <- as.matrix(param)
    }
    
    # vector -> 1-row matrix with column names
    if (is.null(dim(param))) {
      if (is.null(names(param))) {
        cli::cli_abort(
          "`param` is a vector but has no names. Parameter names must match arguments of `qdist`."
        )
      }
      param <- matrix(param, nrow = 1L)
      colnames(param) <- names(param)
    }
    
    if (is.null(colnames(param))) {
      cli::cli_abort(
        "Column names of `param` are missing. They must match arguments of `qdist`."
      )
    }
    
    # Check that all param names are valid arguments of qdist
    if (!all(colnames(param) %in% qformals)) {
      bad <- setdiff(colnames(param), qformals)
      cli::cli_abort(c(
        "Some parameter names in `param` are not arguments of `qdist`.",
        "x" = "Invalid names: {paste(bad, collapse = ', ')}"
      ))
    }
    
    # Recycle single row to all dimensions
    if (nrow(param) == 1L && D > 1L) {
      param <- param[rep(1L, D), , drop = FALSE]
      cli::cli_inform(
        "Single-row `param` has been recycled to all {D} dimensions."
      )
    }
    
    if (nrow(param) != D) {
      cli::cli_abort(c(
        "`param` must have either 1 row or as many rows as dimensions D.",
        "x" = "nrow(param) = {nrow(param)}, D = {D}."
      ))
    }
  }
  
  # --- NorTA core ----------------------------------------------------------
  # 1. multivariate normal
  z <- mvtnorm::rmvnorm(n = n, sigma = cor)
  
  # 2. Gaussian CDF -> uniforms
  u <- stats::pnorm(z)
  
  # 3. Apply marginal quantile
  ans <- matrix(NA_real_, nrow = n, ncol = D)
  colnames(ans) <- colnames(cor)
  
  other_args <- setdiff(qformals, c("p", "..."))
  
  if (is.null(param)) {
    # Use default parameters of qdist
    for (d in seq_len(D)) {
      ans[, d] <- qdist(p = u[, d])
    }
  } else {
    for (d in seq_len(D)) {
      arglist <- list(p = u[, d])
      
      for (par_name in other_args) {
        if (par_name %in% colnames(param)) {
          arglist[[par_name]] <- param[d, par_name]
        }
      }
      
      ans[, d] <- do.call(qdist, arglist)
    }
  }
  
  ans
}
