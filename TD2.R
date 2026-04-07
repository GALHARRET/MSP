data <- read.csv("~/Downloads/exemple_carte_Shewhart.csv", 
                 row.names=1, stringsAsFactors=TRUE)

dim(data)

summary(data)

# Selection des 10 premiers prélèvements

data1=data[1:10,]

mu=apply(data1,1,mean)
M=mean(mu)

n=5
library(multiSPC)
sigma=apply(data1,1,sd)
S=mean(sigma)/c4(n)

LIC=M-3*S/sqrt(n)
LSC=M+3*S/sqrt(n)

plot_chart(mu,LIC=LIC,LSC=LSC)



### Sur les autres prélèv.
data2=data[11:24,]
mu=apply(data2,1,mean)
plot_chart(mu,LIC=LIC,LSC=LSC)

sig_S=(sqrt(1-c4(n)^2)/c4(n))*S

lic=max(S-3*sig_S,0)
lsc=S+3*sig_S

g1=plot_chart(sigma,LIC=lic,LSC=lsc)

### Sur les autres prélèv.
sig=apply(data2,1,sd)
g2=plot_chart(sig,LIC=lic,LSC=lsc,Type = "carte de l'écart type")

library(gridExtra)
grid.arrange(g1,g2,nrow=2)

load("~/Downloads/X.rda")
dim(X)

graph1=list()
graph2=list()

for(i in 1:dim(X)[2]){
  data1=X[1:15,i,]
  mu=apply(data1,1,mean)
  M=mean(mu)
  
  n=6
  sigma=apply(data1,1,sd)
  S=mean(sigma)/c4(n)
  
  LIC=M-3*S/sqrt(n)
  LSC=M+3*S/sqrt(n)
  
  graph1[[i]]=plot_chart(mu,LIC=LIC,LSC=LSC,
                         Type=paste("carte moy. (PhaseI) pour caract.",i))
  
  data2=X[16:40,i,]
  mu=apply(data2,1,mean)
  graph2[[i]]=plot_chart(mu,LIC=LIC,LSC=LSC,
                         Type=paste("carte moy. (PhaseII) pour caract.",i))
}

library(gridExtra)

p=grid.arrange(graph1[[1]],graph2[[1]],
             graph1[[2]],graph2[[2]],
             graph1[[3]],graph2[[3]],
             graph1[[4]],graph2[[4]],
             graph1[[5]],graph2[[5]],ncol=2)


Moy_T(X)

covariance_X(X)

## Construction de la carte
M=Moy_T(X[1:15,,])
S=covariance_X(X[1:15,,])
T2_I=T2_Hotelling_kn(X[1:15,,])
LSC=LSC_T2_Hotelling_kn(X[1:15,,])

plot_chart(T2_I,LSC=LSC[1],Type="carte de T2 (Phase I)")

T2_II=T2_Hotelling_kn(X[16:40,,],MoyT = M,S = S)

plot_chart(T2_II,LSC=LSC[2],Type="carte de T2 (Phase II)")


M=Moy_T(X[1:15,,])
S=covariance_X(X[1:15,,])
T2_I=T2_Hotelling_kn(X[1:15,,])
LSC=LSC_T2_Hotelling_kn(X[1:15,,])

plot_chart(T2_I,LSC=LSC[1],Type="carte de T2 (Phase I)")

T2_II=T2_Hotelling_kn(X[16:40,,],MoyT = M,S = S)

plot_chart(T2_II,LSC=LSC[2],Type="carte de T2 (Phase II)")

####
T2_Hotelling_kn_MYT(k_hc=5,k=15,select=1,X[16:40,,],MoyT = M)


T2_Hotelling_kn_MYT(k_hc=12,k=15,select=1,X[16:40,,],MoyT = M)
