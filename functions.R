#' Fonctions pour l'implémentation de l'ACP sparse
#' @author R. Bordas A. Hollands

#' racine carrée d'une matrice x
sqrt.matrix <- function (x)
{
  x.eigen <- eigen(x)
  d <- x.eigen$values
  d <- (d + abs(d)) / 2
  v <- x.eigen$vectors
  return(v %*% diag(sqrt(d)) %*% t(v))
}

compute.pev <- function(beta_matrix, X_data) {
  # calcule la variance expliquée suivant les Zou et al. 2006
  # les PC dans l'ACP sparse ne sont plus orthogonaux => la variance expliquée
  # tient compte et compense cette covariance non nulle (astuce d'une 
  # décomposition QR, voir section 3.4 de Zou et al. 2006)
  svdobj <- svd(X_data)
  totalvariance <- sum((svdobj$d) ^ 2)
  u <- X_data %*% beta_matrix
  R <- qr.R(qr(u))
  explained.var <- diag(R ^ 2) / totalvariance
  
  return(explained.var)
}

#' Implémentation personnelle de l'ACP sparse d'après Zou et al. 2006
#' 
#' 
#' Retourne une liste avec les matrices A et B optimisées ainsi que le nombre
#' d'itérations effectivement réalisées
#' 
#' @param X les données centrées
#' @param k le nombre de composantes
#' @param lambda les lambda_1j pour la pénalisation L1
#' @param alpha le paramètre du Ridge (voir doc glmnet pour détails)
#' @param max_iteration 
#' @param epsilon 
#' @param stop.criteria norm pour la norm, toutes autres chaines de caractères pour le critère du package elasticnet
#' @param use.cov pour utiliser soit la covariance soit la matrice des données dans l'elastic net
standard.spca <- function(X,
                          k,
                          lambda,
                          alpha,
                          max_iteration,
                          epsilon = 1e-3,
                          stop.criteria = "norm",
                          use.cov = FALSE) {
  # X is assumed to be centered
  # k is the number of components
  # alpha is the ridge param
  # lambda is a vector for lasso param
  
  # initialisation
  A <- svd(X)$v[, 1:k]
  B <- matrix(0L, nrow = dim(X)[2], ncol = k)
  
  # sample covariance to avoid recomputation at each iteration
  covX <- 1 / dim(X)[1] * t(X) %*% X
  if (use.cov) {
    sigma <- sqrt.matrix(covX) # too slow for p >> n ??
  } else {
    sigma <- X
  }
  iter <- 1
  
  c_B <- Inf
  c_A <- Inf

  pb <- txtProgressBar(min = 1, max = max_iteration, style = 3)
  while (iter < max_iteration & c_B > epsilon & c_A > epsilon) {
    b_tmp <- B
    a_tmp <- A
    Y <- sigma %*% A
    # compute B
    for (j in 1:k) {
      # here alpha = lambda_2 / (lambda_1 + lambda_2)
      pc.mdl <-
        glmnet(
          sigma,
          (sigma %*% A)[, j],
          lambda = lambda[j],
          alpha = alpha,
          family = "gaussian"
        )
      B[, j] <- pc.mdl$beta[, 1]
    }
    
    # compute A
    svd.A <- svd(covX %*% B)
    A <- svd.A$u %*% t(svd.A$v)
    
    # checking convergence
    if (stop.criteria == "norm") {
      c_A <- norm(A - a_tmp, type = "F")
      c_B <- norm(B - b_tmp, type = "F")
    } else {
      # this code is identical to convcheck
      # in elasticnet package
      a <- apply(abs(B + b_tmp), 2, max)
      b <- apply(abs(B - b_tmp), 2, max)
      d <- length(a)
      x <- rep(1, d)
      for (i in 1:d) {
        x[i] <- min(a[i], b[i])
      }
      c_B <- max(x)
    }
    iter <- iter + 1
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  
  return(list(A = A,
              B = B,
              n_iter = iter))
}
