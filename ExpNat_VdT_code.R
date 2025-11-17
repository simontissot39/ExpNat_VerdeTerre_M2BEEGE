library(tidyverse)
library(ggplot2)
library(readr)
library(vegan)
library(DescTools)
library(FSA)

Data = read.csv("C:/Users/PC/Desktop/M2/S1/ExpNat/ExpNat_VdT/Data terrain vdt - Feuille 1.csv")

Data_count = Data %>%
  group_by(Placette, Replicat) %>%
  summarise(count = n(),
            Shannon = {
              pi <- table(Famille) / n()
              -sum(ifelse(pi > 0, pi * log(pi), 0))  
            },
            Distance = mean(Distance),
            Prop_Matures = sum(Maturite == "oui", na.rm = TRUE) / n(),
            .groups = "drop")  %>%
  ungroup()


## Taille VS Placette ----

ggplot(data = Data)+
  geom_boxplot(aes(x = Placette,
                   y = Taille.mm.))+
  geom_jitter(aes(x = Placette,
                  y = Taille.mm.),
              color = "blue",
              width = 0.1)+
  labs(x = "Placette", y = "Taille (mm)", title = "Taille des vers selon les différentes placettes")+
  theme_minimal()

tst = aov(data = Data,
          Taille.mm.~Placette)

shapiro.test(tst$residuals)
#Test de Shapiro significatif (p<0.05) = données non normales

kruskal.test(data = Data,
             Taille.mm.~Placette)
#Test de Kruskall.Wallis significatif (p<0.05) = différences entre les moyennes des groupes

dunnTest(data = Data,
         Taille.mm.~Placette,
         method = "bonferroni")
#Test de Dunn significatif pour B>C et D>C

### T vs P - lm -----

ggplot(data = Data)+
  geom_jitter(aes(x = Distance,
                  y = Taille.mm.),
              color = "blue",
              width = 10)+
  geom_vline(aes(xintercept = 18))+
  geom_vline(aes(xintercept = 58))+
  geom_vline(aes(xintercept = 115))+
  geom_vline(aes(xintercept = 174))+
  geom_smooth(aes(x = Distance,
                  y = Taille.mm.),
              se = TRUE,
              method = "lm",
              linetype = "dashed",
              color = "black")+
  labs(x = "Distance", y = "Taille (mm)", title = "Taille des vers selon la distance")+
  theme_minimal()

model_TD = lm(data = Data, Taille.mm.~Distance)
par(mfrow=c(2,2))
plot(model_TD)

summary(model_TD)




## Count VS Placette ----

ggplot(data = Data_count)+
  geom_point(aes(x = Placette,
                 y = count),
             color = "blue")+
labs(x = "Placette", y = "Nombre d'individus", title = "Nombre d'individus aux différentes placettes")+
  theme_minimal()
#Comme décrit ds le rapport = pas de test


### Count VS Placette lm -----

ggplot(data = Data_count)+
  geom_point(aes(x = Distance,
                 y = log10(count)),
             color = "blue")+
  geom_smooth(aes(x = Distance,
                  y = log10(count)),
              se = TRUE,
              method = "lm",
              linetype = "dashed",
              color = "black")+
  labs(x = "Distance à la route (m)", y = "log10(Nombre d'individus)", title = "Nombre d'individus avec éloignement")+
  theme_minimal()

model_CvsP = lm(data = Data_count, log10(count)~Distance)

par(mfrow=c(2,2))
plot(model_CvsP)

summary(model_CvsP)




## Mature VS Placette ----

ggplot(data = Data_count)+
  geom_point(aes(x = Placette,
                 y = Prop_Matures),
             color = "blue")+
labs(x = "Placette", y = "Proportion d'individus matures", title = "Proportion d'individus matures aux différentes placettes")+
  theme_minimal()
#Comme décrit ds le rapport = pas de test


### Maturite vs Placette lm -----

ggplot(data = Data_count)+
  geom_point(aes(x = Distance,
                 y = Prop_Matures),
             color = "blue")+
  geom_smooth(aes(x = Distance,
                  y = Prop_Matures),
              se = TRUE,
              method = "lm",
              linetype = "dashed",
              color = "black")+
  labs(x = "Distance à la route (m)", y = "Proportion d'individus matures", title = "Proportion d'individus matures avec éloignement")+
  theme_minimal()

model_MvsD = lm(data = Data_count, Prop_Matures~Distance)

par(mfrow=c(2,2))
plot(model_MvsD)

summary(model_MvsD)

## Shannon VS Placette ----

ggplot(data = Data_count)+
  geom_point(aes(x = Placette,
                 y = Shannon),
             color = "blue")+
  labs(x = "Placette", y = "Diversité (indice de Shannon)", title = "Diversité aux différentes placettes")+
  theme_minimal()
#Comme décrit ds le rapport = pas de test


### Shannon vs Placette lm -----

ggplot(data = Data_count)+
  geom_point(aes(x = Distance,
                 y = Shannon),
             color = "blue")+
  geom_smooth(aes(x = Distance,
                  y = Shannon),
              se = TRUE,
              method = "lm",
              linetype = "dashed",
              color = "black")+
  labs(x = "Distance à la route (m)", y = "Diversité (indice de Shannon)", title = "Diversité avec éloignement")+
  theme_minimal()

model_SvsD = lm(data = Data_count, Shannon~Distance)

par(mfrow=c(2,2))
plot(model_SvsD)

summary(model_SvsD)

# PCoA communauté ----

