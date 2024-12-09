#' @author R. Bordas A. Hollands

library(rpart)

source("utils.R")

##### loading genes data #####

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

##### preparing data #####

X_train <- scale(as.matrix(data.train$data))
y_train <- data.train$classes
table(y_train) # 0 is ALL ; 1 is AML
X_test <- scale(as.matrix(data.test$data))
y_test <- data.test$classes
table(y_test) # 0 is ALL ; 1 is AML


X_train %>%
  as_tibble() %>%
  rowid_to_column(var="X") %>%
  gather(key="Y", value="Z", -1) %>%
  mutate(Y=as.numeric(gsub("V","",Y))) %>%
  ggplot(aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral") +
  labs(x = "patients", y = "genes") +
  theme_classic() +
  theme(legend.position="none")
ggsave(
  "training_dataset.pdf",
  path = "figures",
  width = 3,
  height = 6
)

X_large_ds <- scale(as.matrix(large.dataset$data))

##### standard PCA #####

pca.1 <- PCA(X_train, ncp = 6)
print(pca.1$eig[1:3, "percentage of variance"])
print(sum(pca.1$eig[1:3, "percentage of variance"]))
print(sum(pca.1$eig[1:2, "percentage of variance"]))

reduced_data <- data.frame(pca.1$ind$coord[, 1:5], Y = factor(y_train))
table(reduced_data$Y)
reduced_data$Y_labels <- factor(reduced_data$Y, levels = c(0, 1), labels = c("ALL", "AML"))

fviz_eig(pca.1, ncp = 20)
ggsave(
  "eigenvalues_standard_pca.pdf",
  path ="figures",
  width=4,
  height=2.5
)
fviz_pca_ind(
  pca.1,
  axes = c(1, 2),
  col.ind = reduced_data$Y_labels,
  palette = c("red", "blue"),
  label = "none"
)
ggsave(
  "individuals_standard_pca.png",
  path ="figures",
  width=4,
  height=4
)

fig <- plot_ly(reduced_data, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, color = ~Y_labels,
               type = "scatter3d")
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                   yaxis = list(title = 'PC2'),
                                   zaxis = list(title = 'PC3')))

fig

as.data.frame(pca.1$svd$V)[,1:3] %>%
  pivot_longer(cols=1:3, names_to="PC") %>%
  mutate(PC = str_replace_all(PC, "V", "PC")) %>%
  ggplot(aes(x = value)) +
  geom_vline(xintercept = 0, colour="red") +
  geom_histogram(aes(y = after_stat(density)), bins = 40) +
  scale_x_continuous(breaks = c(-0.025, 0.025), limits = c(-0.04, 0.04)) + 
  scale_y_continuous(limits = c(0, 82)) +
  facet_wrap(PC ~ .) +
  xlab("PC loadings") +
  theme_classic()
ggsave(
  "loadings_standard_pca.pdf",
  path ="figures",
  width=4,
  height=1.5
)

##### classification #####

pos_weights = 1.0 / (nrow(subset(reduced_data, Y == 1)) / nrow(reduced_data))
neg_weights = 1.0 / (nrow(subset(reduced_data, Y == 0)) / nrow(reduced_data))

case_weights <- ifelse(reduced_data$Y == 1, pos_weights, neg_weights)

res.1 <- glm(Y ~ Dim.1 + Dim.2 + Dim.3, family = binomial(link = "logit"), data = reduced_data)
summary(res.1)


# classification on standard PCA dimension reduction
pcaTree <- rpart(Y ~ Dim.1 + Dim.2 + Dim.3,
                 data = reduced_data,
                 control = rpart.control(minsplit = 5, cp = 0),
                 weights = case_weights)

X_test_pca <- predict(pca.1, newdata = X_test)$coord[, 1:3]

predictions <- predict(pcaTree, newdata = data.frame(X_test_pca), type = "class")

confusion_matrix <- table(predictions, y_test)
print(confusion_matrix)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(accuracy)

###### cette section utilise des variables définies dans custom_sparsepca
# TODO : MAJ => la présentation utilisait les résultats du package elasticnet de Zou et al. 2006
# Le code ci-dessous utilise notre implémentation personnelle => résultats légèrement différents 
# (car paramètres des pénalisations différents)
data_to_plot <- data.frame(X_train %*% testing.spca1$B, Y_labels = as.factor(y_train))
data_to_plot$Y_labels <- factor(reduced_data$Y, levels = c(0, 1), labels = c("ALL", "AML"))
X_test_spca <- as.data.frame(X_test %*% testing.spca1$B[, 1:3])

pcaTree_sparse <- rpart(Y_labels ~ X1 + X2 + X3,
                 data = data_to_plot,
                 control = rpart.control(minsplit = 5, cp = 0),
                 weights = case_weights)

colnames(X_test_spca) <- c("X1", "X2", "X3")

predictions <- predict(pcaTree_sparse, newdata = X_test_spca, type = "class")

confusion_matrix <- table(predictions, y_test)
print(confusion_matrix)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(accuracy)


