library(SixSigma)
library(readxl)
predict_parafac <- function(X_new, model) {
  # X_new : tenseur (I_new x J x K)
  # model : liste {A, B, C} du modèle PARAFAC déjà ajusté
  # Retourne : A_new, X_hat, SPE, R2
  
  if (length(dim(X_new)) != 3) stop("X_new doit être un tenseur 3D.")
  
  A <- model$A
  B <- model$B
  C <- model$C
  
  I_new <- dim(X_new)[1]
  J <- dim(X_new)[2]
  K <- dim(X_new)[3]
  F <- ncol(B)
  
  # Produit de Khatri–Rao : C ⊙ B
  khatri_rao <- function(M1, M2) {
    do.call(cbind, lapply(1:ncol(M1), function(f) kronecker(M1[,f], M2[,f])))
  }
  
  # Déroulement mode-1
  X1 <- matrix(aperm(X_new, c(1,2,3)), nrow = I_new)
  # Scores A_new estimés par régression
  KR <- khatri_rao(C, B)
  A_new <- X1 %*% KR %*% solve(t(KR) %*% KR)
  # Reconstruction
  X_hat <- fitted.parafac(list(A = A_new, B = B, C = C))
  
  # Erreur de prédiction (SPE)
  SPE <- sapply(1:I_new, function(i) sum((X_new[i,,] - X_hat[i,,])^2))
  
  # Variance expliquée (R²)
  var_total <- apply(X_new, 1, function(x) sum((x - mean(x))^2))
  R2 <- 1 - SPE / var_total
  
  return(list(A_new = A_new, X_hat = X_hat, SPE = SPE, R2 = R2))
}

dta=read.csv("ex_qPCR.csv",sep=";",dec=",")
str(dta)

? ss.rr

# modele à 2 facteurs souche x operateur
table(dta$souche,dta$operateur)

res_RR=ss.rr(reponse,souche,operateur,data=dta,print_plot=FALSE)

# autres modèles d'anova

# modele complet
anova(lm(reponse~souche+operateur, data=dta))


# modele a 3 facteurs, sans interaction
anova(lm(reponse~souche+operateur+jour, data=dta))

model=aov(reponse~souche, data=dta)
modelm <- summary(model)
modelm[[1]]


library(multiSPC)
library(multiway)
X=carbon1
T2=T2_Hotelling_kn(X)
LSC=LSC_T2_Hotelling_kn(X)
p1=plot_chart(T2,LSC=LSC[1],Type="T2 Phase I")
plot_chart(T2,LSC=LSC[1])
X2=carbon2
T2=T2_Hotelling_kn(X2,MoyT = Moy_T(X),S=covariance_X(X))

p2=plot_chart(T2,LSC=LSC[2],Type="T2 Phase II")

k=dim(X)[1]
X1=matrix(NA,nrow=k,ncol=p)
for(t in 1:k){
  X1[t,]=t(apply(X[t,,],1,mean))
}
colnames(X1)=dimnames(X)[[2]]
Xc=scale(X1,scale=F)
library(FactoMineR)
pca=PCA(Xc,scale.unit=T)
r <- 2
scores <- pca$ind$coord[,1:r]
#S=cov(diff(scores))/2
S=cov(scores)
T2pca=T2_Hotelling_k1(scores,S=S)
LSCpca=LSC_T2_Hotelling_k1(scores)
p3=plot_chart(T2pca,LSC=LSCpca["PhaseI"],Type="T2 ACP Phase I")
k=dim(X2)[1]
X1=matrix(NA,nrow=k,ncol=p)

for(t in 1:k){
  X1[t,]=t(apply(X2[t,,],1,mean))
}
colnames(X1)=dimnames(X)[[2]]
Xc=scale(X1,scale=F)

scores2=predict.PCA(pca,newdata=Xc)$coord[,1:r]
T2pca=T2_Hotelling_k1(scores2,S=S)
p4=plot_chart(T2pca,LSC=LSCpca["PhaseII"],Type="T2 ACP Phase II")

gridExtra::grid.arrange(p1,p2,p3,p4,ncol=2)

model <- parafac(X, nfac = 2, nstart = 5)


res=predict_parafac(X2,model)
# Hotelling T²
T2 <- apply(model$A, 1, function(t) t(t) %*% solve(cov(model$A)) %*% t)

LIC=mean(T2)-3*sd(T2)
LSC=mean(T2)+3*sd(T2)
p5=plot_chart(T2,LIC=LIC,LSC=LSC)


T2 <- apply(res$A_new, 1, function(t) t(t) %*% solve(cov(model$A)) %*% t)

p6=plot_chart(T2,LIC=LIC,LSC=LSC)

gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
