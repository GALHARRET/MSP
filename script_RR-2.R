library(SixSigma)
library(readxl)

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
