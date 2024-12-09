#' Implémentation de l'ACP sparse
#' 
#' Algorithme d'après Zou et al. 2006
#' 
#' Dépendances de ce script : 
#' - utils.R => pré-traitement des données génomiques
#' - functions.R => implémentation de l'ACP sparse
#' 
#' Ce script génère les figures du rapport : 
#' - simulation de données d'ACP
#' - exemple sur des données génomiques
#' 
#' @author R. Bordas A. Hollands
#' 
library(elasticnet)
library(glmnet)
library(mvtnorm)
library(tidyverse)
library(FactoMineR)
library(factoextra)

colors <- c(
  PCA = "#FDAE61",
  SPCA = "#ABDDA4",
  spca.crit2 = "#3288BD",
  norm = "#ABDDA4",
  sparse = "#3288BD"
)

source("utils.R")
source("functions.R")

##### simulations - inspiré de Chavent & Chavent (2020)

eigenval <- c(200, 100, 50, 50, 6, 5, 4, 3, 2, 1)

# True loading matrix
v1 <- c(1, 2, 0.5, 1, 0, 0, 0, 0, 0.9, 0.9)
v2 <- c(0, 0, 0, 0, 1, 4, 0.75, 1, -0.3, 0.3)
v1 <- v1 / norm(v1, type = "2")
v2 <- v2 / norm(v2, type = "2")

sim.PCA.data <- function(eigvects, eigvals, n_, p_)
{
  V_star <- matrix(0L, nrow = p_, ncol = p_)
  m <- length(eigvects)
  for (j in 1:m)
    V_star[, j] <- eigvects[[j]]

  set.seed(42)
  V_star[, (m + 1):p_] <- matrix(runif(p_ * (p_ - m)), nrow = p_, ncol = p_ - m)
  
  # matrice orthogonale selon Gram-Schmidt (i.e. décomposition QR)
  V <- qr.Q(qr(V_star))
  V <- round(V, digits = 12) # pour éviter erreurs d'arrondis
  
  C <- V %*% diag(eigvals) %*% t(V) / n_
  print(dim(C))
  set.seed(42)
  A <- mvtnorm::rmvnorm(n_,sigma=C)
  
  return(list(
    A = A,
    Vtrue = V
  ))
}

sim.data <- sim.PCA.data(list(v1, v2), eigenval, n_ = 100, p_ = 10)
print(dim(sim.data$A))
xtable(sim.data$Vtrue[, 1:2], digits = 3)
sim.pca.res <- PCA(sim.data$A, ncp = 2)
print(sim.pca.res$svd$V[, 1:2])
xtable(sim.pca.res$svd$V[, 1:2], digits = 3)
spca.sim <- standard.spca(
  sim.data$A,
  k = 2,
  lambda = c(0.1, 0.1),
  stop.criteria = "norm",
  max_iteration = 200,
  alpha = 0.5
)
print(sim.data$V[, 1:2])
print(spca.sim$B)
xtable(spca.sim$B, digits = 3)


##### données génomiques (Golub et al. 1999) #####

data.train <- load_gene_data(
  filename = "../data/all_aml/all_aml_train.csv", 
  output_name = "../data/all_aml/all_aml_train_preprocessed.csv"
)
print(dim(data.train$data)) # 38 x 7129

data.test <- load_gene_data(
  filename = "../data/all_aml/all_aml_test.csv", 
  output_name = "../data/all_aml/all_aml_test_preprocessed.csv"
)
print(dim(data.test$data)) # 35 x 7129

##### LARS-EN pour trouver des lambdas_1j optimaux sur les données génomiques
pca.1 <- PCA(X_train, ncp = 37)
V <- pca.1$svd$V

lambdas <- list()

for (i in 1:5) {
  # LARS-EN implémenté dans elasticnet::enet
  mdl.enet <- enet(V[, -i], V[, i], lambda = 4, trace = FALSE, intercept = TRUE)
  candidates <- mdl.enet$L1norm[mdl.enet$L1norm > 0]
  lambdas[[i]] <- list(c(min(candidates), max(candidates)))
}
lambdas

###### SPCA avec le critère d'arrêt de la norme de A et B
# using lambda_1,j from LARS-EN optim
start_time <- Sys.time()
testing.spca1 <- standard.spca(
  X_train,
  k = 3,
  # correcting for the way glmnet defines penalties : 
  # lambda = c(0.07813366, 0.02340993, 0.05277198) / 0.5,
  lambda = c(0.44479509, 0.75967155, 0.36226425),
  stop.criteria = "norm",
  max_iteration = 200,
  alpha = 0.5
)
end_time <- Sys.time()
print(end_time - start_time)
pev.1 <- compute.pev(testing.spca1$B, X_train)
print(pev1)
print(colSums(testing.spca1$B != 0))

pca.1 <- PCA(X_train, ncp = 6)
print(pca.1$eig[1:3, "percentage of variance"])
print(sum(pca.1$eig[1:3, "percentage of variance"]))
print(sum(pca.1$eig[1:2, "percentage of variance"]))

res.df <- data.frame(
  "method" = c("PCA", "SPCA"),
  "total_var" = c(sum(pca.1$eig[1:3, "percentage of variance"]), sum(pev))
)

res.df.per_pc <- data.frame(
  "method" = as.factor(c(rep("PCA", 3), rep("SPCA", 3))),
  "pc" = as.factor(c(rep(c("PC1", "PC2", "PC3"), 2))),
  "var" = c(pca.1$eig[1:3, "percentage of variance"], pev1 * 100)
)
print(res.df.per_pc)

res.df.per_pc %>%
  ggplot(aes(x = pc, y = var, fill = method)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = round(var, 2)),vjust = 1.5,
            position = position_dodge(.9), size = 2) +
  labs(x = "principal components", y = "explained variance (%)") +
  scale_fill_manual(values = colors) + 
  theme_classic() +
  theme(text = element_text(size = 8))
ggsave(
  "pev_spca_pca.pdf",
  path ="figures",
  width=3,
  height=1.5
)

data_to_plot <- data.frame(X_train %*% testing.spca1$B, Y_labels = as.factor(y_train))
data_to_plot$Y_labels <- factor(reduced_data$Y, levels = c(0, 1), labels = c("ALL", "AML"))
fig <- plot_ly(data_to_plot, x = ~X1, y = ~X2, z = ~X3, color = ~Y_labels,
               type = "scatter3d")
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

###### ACP sparse avec le critère d'arrêt sparse (voir elasticnet package)
# using lambda_1,j from LARS-EN optim
testing.spca1.crit <- standard.spca(
  X_train,
  k = 3,
  # lambda = c(0.5, 0.35, 0.85),
  # correcting for the way glmnet defines penalties : 
  # lambda = c(0.07813366, 0.02340993, 0.05277198) / 0.5, 
  lambda = c(0.44479509, 0.75967155, 0.36226425),
  stop.criteria = "other.stopping",
  max_iteration = 200,
  alpha = 0.5
)
svdobj <- svd(X_train)
totalvariance <- sum((svdobj$d) ^ 2)
u <- X_train %*% testing.spca1.crit$B
R <- qr.R(qr(u))
pev2 <- diag(R ^ 2) / totalvariance
print(pev2)
print(colSums(testing.spca1.crit$B != 0))

data_to_plot <- data.frame(X_train %*% testing.spca1.crit$B, Y_labels = as.factor(y_train))
fig <- plot_ly(data_to_plot, x = ~X1, y = ~X2, z = ~X3, color = ~Y_labels,
               type = "scatter3d")
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

res.sparse.per_pc <- data.frame(
  "stopping_criteria" = as.factor(c(rep("norm", 3), rep("sparse", 3))),
  "pc" = as.factor(c(rep(c("PC1", "PC2", "PC3"), 2))),
  "var" = c(pev1 * 100, pev2 * 100),
  "sparsity_lvl" = c(
    colSums(testing.spca1$B != 0),
    colSums(testing.spca1.crit$B != 0)
  )
)
print(res.sparse.per_pc)
res.sparse.per_pc %>%
  ggplot(aes(x = pc, y = sparsity_lvl, fill = stopping_criteria)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "principal components", y = "non-zero coordinates", fill = "stopping crit.") +
  scale_fill_manual(values = colors) + 
  theme_classic() +
  theme(text = element_text(size = 8))
ggsave(
  "sparsity_counts_spca.pdf",
  path ="figures",
  width=3,
  height=1.5
)

res.sparse.per_pc %>%
  ggplot(aes(x = pc, y = var, fill = stopping_criteria)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = round(var, 2)),vjust = 1.5,
            position = position_dodge(.9), size = 2) +
  labs(x = "principal components", y = "explained variance (%)", fill = "stopping crit.") +
  scale_fill_manual(values = colors) + 
  theme_classic() +
  theme(text = element_text(size = 8))
ggsave(
  "variance_stopping_crit_spca.pdf",
  path ="figures",
  width=3,
  height=1.5
)


##### ACP sparse vs ACP sparse pour données de puces à ADN
sparsity_lvls <- matrix(0L, nrow = 8, ncol = 3)
alphas <- c(0.5, 0.2, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6)
for (i in 1:8) {
  testing.spca1.crit <- standard.spca(
    X_train,
    k = 3,
    # correcting for the way glmnet defines penalties : 
    lambda = c(0.07813366, 0.02340993, 0.05277198) / 0.5, 
    stop.criteria = "norm",
    max_iteration = 200,
    alpha = alphas[i]
  )
  sparsity_lvls[i, ] <- colSums(testing.spca1.crit$B != 0)
}

pdf(file = "figures/sparsity_lvls_microarray_custom_code.pdf", width = 4, height = 3)
par(mar = c(4, 5, 1, 1))
plot(
  1 - alphas,
  sparsity_lvls[, 1] / 7129,
  log = "x",
  lty = 2,
  type = "b",
  xlab = TeX(r"($\alpha$)"),
  ylab = TeX(r"(|$\beta_{1}|_0$ / p)"),
  ylim = c(0, 1)
)
abline(v = 0.5, col = "red", lty = 1)
#arrows(0.65, 0.4, 0.75)
# text(x = c(0.7), y = 0.3, TeX(r"($\alpha \rightarrow 1$)"), adj = c(0.5, 0.5))
dev.off()


