library(tidyverse)
library(ggplot2)
library(readr)
library(vegan)
library(DescTools)
library(FSA)

# Communautés selon distance ----

Data = read.csv("...path.../Data terrain vdt - Feuille 1.csv")

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

# Communauté ----

df_prop = Data


df_prop$espece <- replace(df_prop$Famille,
                          grepl("anec", df_prop$Famille, ignore.case = TRUE),
                          "anécique") 

df_prop <- df_prop %>%
  select(-Famille) %>%
  group_by(Placette, Replicat, espece) %>% #, Maturite) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Placette, Replicat) %>%
  mutate(prop = n / sum(n)) %>%
  filter(espece != "ind")

colnames(df_prop) = c("Placette", "replicat", "espece", "n", "prop" ) #, "maturite")

df_prop = left_join(
  df_prop,
  Data %>%
    select(Placette,Distance) %>%
    distinct(Placette, Distance),
  by = "Placette" 
)


### Communauté barplot ----

ggplot(df_prop,
       aes(x = replicat,
           y = prop,
           fill = espece #,       # couleur = espèce
          # alpha = maturite)
          )) + # luminosité = maturité
  geom_bar(stat = "identity", color = "black") +
  
  facet_wrap(~ Placette, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format()) +
  
  # ---- Échelle couleur par espèce (palette standard) ----
scale_fill_brewer(name = "Gp. fonctionel", palette = "Set2") +
  
  # ---- Échelle alpha séparée pour maturité ----
#scale_alpha_manual(
#  name = "Maturité",
 # values = c("non" = 0.5,  # juvénile = clair
  #           "oui" = 1),   # mature = foncé
  #labels = c("non" = "Juvénile",
   #          "oui" = "Mature")
#) +
  
  labs(
    x = "Réplicat",
    y = "Proportion",
    title = "Composition spécifique par réplicat et placette" #\n(clair = juvénile, foncé = mature)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold")
  )


### Comm. typique barplot ----

df_prop_mean = df_prop %>%
  group_by(Placette, espece) %>%
  mutate(prop_mean = sum(prop)/3) %>%
  ungroup %>%
  distinct(Placette, espece, .keep_all = TRUE) %>%
  select(Placette, espece, prop_mean)


ggplot(df_prop_mean,
       aes(x = Placette,
           y = prop_mean,
           fill = espece #,       # couleur = espèce
           # alpha = maturite)
       )) + # luminosité = maturité
  geom_bar(stat = "identity", color = "black") +
  labs(
    x = "Réplicat",
    y = "Proportion",
    title = "Composition spécifique moyenne par réplicat et placette" #\n(clair = juvénile, foncé = mature)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold")
  )


### Comm. modèle linéaire ----

ggplot(data = df_prop %>%
         filter(espece != "ind"))+
  geom_point(aes(x = Distance,
                 y = prop,
                 colour = espece))+
  geom_smooth(aes(x = Distance,
                  y = prop,
                  colour = espece,
                  fill = espece),
              se = TRUE,
              method = "lm")+ 
  facet_wrap(~ espece)+
  labs(
    x = "Distance",
    y = "Proportion",
    title = "Composition spécifique par réplicat et placette",
    color = "Groupe fonctionel",
    fill = "Groupe fonctionel"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold")
  )

mdl_prop_dist = lm(data = df_prop,
                   prop ~ Distance*espece)

par(mfrow = c(2,2))
plot(mdl_prop_dist)
summary(mdl_prop_dist)


#Vibrations dans le sol ----

folder = "...path.../Moutarde"
files = list.files(folder, pattern = "\\.spec$",
                   full.names = TRUE)
data_list = lapply(files, read.table, header = FALSE)
names(data_list) = tools::file_path_sans_ext(basename(files))

data_list = lapply(data_list, function(df) {
  colnames(df) = c("freq_Hz", "freq_mean", "freq_sd_inf", "freq_sd_sup")  
  df
}) 

data_list <- Map(function(df, nom) {
  df$nom_dataset <- nom  
  df
}, data_list, names(data_list))

data_list = lapply(data_list, function(df) {
  df$direction = substr(df$nom_dataset, nchar(df$nom_dataset), nchar(df$nom_dataset))  
  df$malette = substr(df$nom_dataset, nchar(df$nom_dataset)-3, nchar(df$nom_dataset)-2)
  df
}) 

data_vibrations = bind_rows(data_list)

malettes = data.frame(
  malette = c("43", "40", "42", "45"),
  Placette = c("A", "B", "C", "D")
)

data_vibrations = left_join(
  data_vibrations, malettes, by = "malette"
)


ggplot(data = data_vibrations)+
  geom_ribbon(aes(x = freq_Hz,
                  ymin = freq_sd_inf,
                  ymax = freq_sd_sup,
                  fill = Placette),
              alpha = 0.2)+
  geom_line(aes(x = freq_Hz,
                y = freq_mean,
                color = Placette))+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~ direction)+
  labs(
    x = "Frequence (Hz)",
    y = "Power spectral density ((m/s)²/Hz)",
    title = "Spectre de puissance des vibration terrestres des placettes"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold")
  )






