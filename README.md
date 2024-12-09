# M2 Mathématiques et IA - Méthodes statistiques non supervisées

Auteurs : R. Bordas ; A. Hollands

Référence du projet : 

Zou, H., Hastie, T., & Tibshirani, R. (2006). Sparse Principal Component Analysis. Journal of Computational and Graphical Statistics, 15(2), 265–286. https://doi.org/10.1198/106186006X113430

Structure du code : 

- `functions.R` contient la fonction principale `standard.spca` qui implémente l'algorithme 1 du rapport (+ des fonctions utilitaires annexes).
- `example_sparse_pca.R` contient le code des figures du rapport, la simulation des données et les comparaisons des différentes SPCA.
- `classif_models.R` donne un exemple de classification de données génomiques (à titre illustratif, l'arbre de décision n'a pas fait l'objet de fine-tuning).
- `utils.R` contient le pré-traitement des données génomiques.
