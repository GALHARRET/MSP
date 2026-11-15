library(multiSPC)
library(ggplot2)
library(dplyr)


## Cartes univariées

### Cartes sur les moyennes

df<-read.csv("cap_data.csv",sep=",")
View(df)
Moy<-df |> group_by(sample) |> summarise("M_j"=mean(obs)) 
SD<- df|>group_by(sample) |> summarise("S_j"=mean(obs))

