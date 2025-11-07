# Génération de données aléatoires (10 individus, 5 variables)
set.seed(123)
X <- matrix(rnorm(50), ncol=5)

# Centrage des données
X_centered <- scale(X, center=TRUE, scale=FALSE)

# ACP
pca <- prcomp(X_centered, center=FALSE, scale.=FALSE)
pca<-PCA(X_centered,scale.unit = F)

# On garde les 2 premières composantes principales
r <- 2
U_r <- pca$svd$V[, 1:r]  # Vecteurs propres
X_mean <- attr(X_centered, "scaled:center")  # Moyenne

# Nouvelles données (3 individus, 5 variables)
X_new <- matrix(rnorm(15), ncol=5)

# Centrage des nouvelles données
X_new_centered <- scale(X_new, center=X_mean, scale=FALSE)

# Projection sur les composantes principales
F_new <- X_new_centered %*% U_k

# Reconstruction
X_new_reconstructed <- F_new %*% t(U_k) + rep(X_mean, each=nrow(X_new))

# Affichage
print("Données originales (nouvelles) :")
print(X_new)
print("Données reconstruites :")
print(X_new_reconstructed)
