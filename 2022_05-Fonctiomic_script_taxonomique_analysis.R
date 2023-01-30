
#####################################################################
##   "FonctioMic - Analyse Diversite Taxonomique microbienne"      ##
##            __ Samuel Desquiet / Alphonse THIAW __               ##
#####################################################################

library(tidyverse)
library(phyloseq)
library(ape)
library(hillR)
library(DESeq2)
library(rstatix)
library(ggpubr)
library(kableExtra)
library(ggsci) 

###############################
##  LOADING AND PARSING DATA ##
###############################

## LOAD DES DONNEES BRUTES RAREFIES, ENREGISTREES SOUS FORMES DE RData
#Donnees taxonomiques de comptage des lectures brutes et de comptage des génes predicts
metadata = read_rds("./data/sols_physeq.rds")
load("./data/2022_05_Fonctiomic_taxonomicData.RData")

`%notin%` <- Negate(`%in%`) #Creer une fonction contraire de %in%

## SPLIT DES DONNEES SELON LES CONDITIONS D'ETUDES (culture/prairie et ou controle/amendé) POUR CHAQUE REGNE
#Archées
arch.c = subset_samples(archae.rarefied, MODE.USAGE=="Culture") %>% t()
arch.p = subset_samples(archae.rarefied, MODE.USAGE=="Prairie") %>% t()
arch.cc = subset_samples(archae.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="CONTROL") %>% t()
arch.cr = subset_samples(archae.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="AMENDE") %>% t()
arch.pc = subset_samples(archae.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="CONTROL") %>% t()
arch.pr = subset_samples(archae.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="AMENDE") %>% t()
#Bacteries
bact.c = subset_samples(bact.rarefied, MODE.USAGE=="Culture") %>% t()
bact.p = subset_samples(bact.rarefied, MODE.USAGE=="Prairie") %>% t()
bact.cc = subset_samples(bact.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="CONTROL") %>% t()
bact.cr = subset_samples(bact.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="AMENDE") %>% t()
bact.pc = subset_samples(bact.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="CONTROL") %>% t()
bact.pr = subset_samples(bact.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="AMENDE") %>% t()
#Eucaryotes
fungi.c = subset_samples(fungi.rarefied, MODE.USAGE=="Culture") %>% t()
fungi.p = subset_samples(fungi.rarefied, MODE.USAGE=="Prairie") %>% t()
fungi.cc = subset_samples(fungi.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="CONTROL") %>% t()
fungi.cr = subset_samples(fungi.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="AMENDE") %>% t()
fungi.pc = subset_samples(fungi.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="CONTROL") %>% t()
fungi.pr = subset_samples(fungi.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="AMENDE") %>% t()
#Protistes
prot.c = subset_samples(prot.rarefied, MODE.USAGE=="Culture") %>% t()
prot.p = subset_samples(prot.rarefied, MODE.USAGE=="Prairie") %>% t()
prot.cc = subset_samples(prot.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="CONTROL") %>% t()
prot.cr = subset_samples(prot.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="AMENDE") %>% t()
prot.pc = subset_samples(prot.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="CONTROL") %>% t()
prot.pr = subset_samples(prot.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="AMENDE") %>% t()
#Virus
virus.c = subset_samples(virus.rarefied, MODE.USAGE=="Culture") %>% t()
virus.p = subset_samples(virus.rarefied, MODE.USAGE=="Prairie") %>% t()
virus.cc = subset_samples(virus.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="CONTROL") %>% t()
virus.cr = subset_samples(virus.rarefied, MODE.USAGE=="Culture" & AMENDEMENT.BLE=="AMENDE") %>% t()
virus.pc = subset_samples(virus.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="CONTROL") %>% t()
virus.pr = subset_samples(virus.rarefied, MODE.USAGE=="Prairie" & AMENDEMENT.BLE=="AMENDE") %>% t()


###########################
## ANALYSES TQXONOMIQUES ##
###########################

## I- COURBES D'ACCUMULATION

sac_compute = function(data){
  data.acc = vegan::specaccum(otu_table(data), method = "exact", permutations = 100,
                              conditioned = TRUE, gamma("jack1"), w=NULL )
  return(data.acc)
}

list_phylum=list(arch.c,arch.p,bact.c,bact.p,fungi.c,fungi.p,prot.c,prot.p,virus.c,virus.p)
res=lapply(list_phylum, sac_compute)
names(res)=c("arch.c","arch.p","bact.c","bact.p","fungi.c","fungi.p","prot.c","prot.p","virus.c","virus.p")

plot_sac = function(maliste, Names, title, x, y) {
  png(paste("./results/courbeAccumulation/", title, ".png", sep=""))
  # calcule des limites de l'axe y requise
  l = maliste
  ylm = range(sapply(l, "[[", "richness") + 
                sapply(l, "[[", "sd") * c(-2, 2))
  
  # Fonction pour le plotting des courbes d'accumulation de chaque Règne 
  sapply(seq_along(l), function(i) {
    if (i==1) { #  Si 1ère élément de la liste, utiliser la fonction "plot"
      with(l[[i]], {
        plot(sites, richness, type='l', ylim=ylm, 
             xlab="Sites", ylab="exact", las=1)
        segments(seq_len(max(sites)), y0=richness - 2*sd, 
                 y1=richness + 2*sd)
      })    
    } else {
      with(l[[i]], { # Pour les éléments suivants, utilisés "lines()"
        lines(sites, richness, col=i)
        segments(seq_len(max(sites)), y0=richness - 2*sd, 
                 y1=richness + 2*sd, col=i)
      })     
    }
  })
  
  legend("bottomright", Names, lty=1, col = x:y, bty="n", inset=0.025, cex = 0.75)
  dev.off()
}

#Dossier pour les figures des courbes d'accumulation
dir.create("./results")
dir.create("./results/courbeAccumulation")
# Courbe d'accumulation des Archae
archee=plot_sac(list(res$arch.c, res$arch.p), c("archae.Culture", "archae.Prairie"), "Archees", 1, 2)
# Courbe d'accumulation des Bactéries
plot_sac(list(res$bact.c, res$bact.p),c("bactérie.Culture", "bactérie.Prairie"), "Bacterie", 1, 2)
# Courbe d'accumulation des Champignons, Protistes et des Virus
plot_sac(list(res$virus.c, res$virus.p, res$prot.c, res$prot.p, res$fungi.c, res$fungi.p),
         c("virus.C","virus.P","protiste.C","protiste.P","fungi.C","fungi.P"), "Eucaryote_Virus", 1,6) 



## II IMPACT DES PRATIQUES AGRICOLES SUR LA BIODIVERSITÉ MICROBIENNE DU SOL
### 1. Analyse de biomasse entre culture et prairie

dir.create("./results/modeUsageImpact")
#Abondance des Règnes au temps T0
data_regne.t0 = subset_samples(metadata.rarefied, TIME=="T0")%>%
  tax_table(.)%>%as.data.frame()%>%group_by(Kingdom)%>%
  summarise(count=n())%>%filter(count!=1 )

data_regne.t0%>%ggplot(aes(x=Kingdom, y=count))+
  geom_bar(stat="identity", color="grey20")
ggsave("./results/modeUsageImpact/T0_phylumAbondance.png")

### Difference d'abondance des genres entre prairie et culture
metadata.culture = subset_samples(metadata.rarefied, MODE.USAGE=="Culture" & TIME=="T0")
metadata.prairie = subset_samples(metadata.rarefied, MODE.USAGE=="Prairie" & TIME=="T0")

genus.culture=as.data.frame(tax_table(metadata.culture))$Genus #Récup des noms des OTU(especes) dans le table taxonomique

# On à des triplicats donc on prend les otus present dans au moins un replicat
data_genus.culture = as.data.frame(otu_table(metadata.culture))%>%cbind(Genus=genus.culture)%>%
  filter(`C1-C-T0`!=0 & `C3-C-T0`!=0 & `C5-C-T0`!=0 & Genus !="NA")#%>%select(Species)

genus.prairie=as.data.frame(tax_table(metadata.prairie))$Genus
data_genus.prairie = as.data.frame(otu_table(metadata.prairie))%>%cbind(Genus=genus.prairie)%>%
  filter(`P1-C-T0`!=0 & `P3-C-T0`!=0 & `P5-C-T0`!=0 & Genus!="NA")%>%select(Genus)

data_genus.list = list(Culture=data_genus.culture$Genus, Prairie=data_genus.prairie$Genus)

venn_data.taxo = ggVennDiagram::Venn(data_genus.list)
data.gvd = ggVennDiagram::process_data(venn_data.taxo)
ggplot()+
  geom_sf(aes(fill=count), data = ggVennDiagram::venn_region(data.gvd))+
  geom_sf_text(aes(label=name,), fontface="bold", size=10, data = ggVennDiagram::venn_setlabel(data.gvd))+
  geom_sf_label(aes(label=count), fontface="bold", size=8, data=ggVennDiagram::venn_region(data.gvd))+
  scale_fill_distiller(palette = "GnBu")+
  theme_void(base_size = 20)
ggsave("./results/modeUsageImpact/Venn_diagram.pdf", height = 8, width = 14)


### 2. Analyse d'abondance relative entre culture et prairie
#Identifier les Phyla majoritaire en culture et en prairie et leur proportion à TO
#Pour les cultures
abundance.culture = data.frame(as.data.frame(otu_table(metadata.culture))%>%
                                 mutate(mean_count_T0=rowMeans(select(as.data.frame(otu_table(metadata.culture)), c(`C1-C-T0`,`C3-C-T0`,`C5-C-T0`)))),
                               as.data.frame(tax_table(metadata.culture)))%>%
  group_by(Phylum)%>%
  summarise(culture=sum(mean_count_T0))%>%
  mutate(Culture=round(culture/sum(culture)*100, 2))%>%arrange(desc(Culture))

abundance.culture.top20=abundance.culture[1:22,]%>%filter(Phylum  %notin% c("Unclassified", "NA"))

abundance.prairie = data.frame(as.data.frame(otu_table(metadata.prairie))%>%
                                 mutate(mean_count_T0=rowMeans(select(as.data.frame(otu_table(metadata.prairie)), c(`P1-C-T0`,`P3-C-T0`,`P5-C-T0`)))),
                               as.data.frame(tax_table(metadata.prairie)))%>%
  group_by(Phylum)%>%
  summarise(prairie=sum(mean_count_T0))%>%
  mutate(Prairie=round(prairie/sum(prairie)*100,2))%>%arrange(desc(Prairie))

abundance.prairie.top20=abundance.prairie[1:22,]%>%filter(Phylum  %notin% c("Unclassified", "NA"))

data_abundance.top20=full_join(abundance.culture.top20,abundance.prairie.top20)
data_abundance.top20[is.na(data_abundance.top20)]=0
data_abundance.top20=data_abundance.top20%>%select(-c("culture","prairie"))%>%
  pivot_longer(cols = -Phylum)
colnames(data_abundance.top20)=c("Phylum","echantillons","abondance relative")

data_abundance.top20%>%
  ggplot(aes(x=Phylum,y=`abondance relative`, fill=echantillons))+
  geom_bar(stat = "identity", position = position_dodge(0.8),width = 0.7, show.legend = T)+
  scale_fill_manual(values = c("#fbbb3f","#01fa03"))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank() ,axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title = element_text(size=8.5), axis.text.x = element_text(angle = 70,  hjust = 1))

ggsave("./results/modeUsageImpact/TopPhylum_Abundace.png", height = 6.5, width = 14)


### 3 Analyse de la diversité microbienne en fonction des conditions d'etudes
#Extraction et transposition des tables d'abondances des OTUs
dir.create("./results/Diversites")
arch.data = archae.rarefied@otu_table %>% t()
bact.data = bact.rarefied@otu_table %>% t()
euk.data = euk.rarefied@otu_table %>% t()
fungi.data = fungi.rarefied@otu_table %>% t()
prot.data = prot.rarefied@otu_table %>% t()
virus.data = virus.rarefied@otu_table %>% t()

#Fonction hilNumber.compute pour calculer les nombres de hill dans chaque sous-ensemble de données spécifique à chaque Règne
regne.list=list(arch.data,bact.data,euk.data,fungi.data,prot.data,virus.data)
hillNumber.compute=function(data, indice){
  regne.q=hill_taxa(data, q=indice, MARGIN = 1, base = exp(1))
}

regne.q0 = sapply(regne.list, hillNumber.compute, indice=0)
colnames(regne.q0)=c("archae.q0","bacterie.q0","eukaryote.q0","fungi.q0","protiste.q0","virus.q0")
regne.q1 = sapply(regne.list, hillNumber.compute, indice=1)
colnames(regne.q1)=c("archae.q1","bacterie.q1","eukaryote.q1","fungi.q1","protiste.q1","virus.q1")
regne.q2 = sapply(regne.list, hillNumber.compute, indice=2)
colnames(regne.q2)=c("archae.q2","bacterie.q2","eukaryote.q2","fungi.q2","protiste.q2","virus.q2")

hillNumber.df = data.frame(regne.q0,regne.q1,regne.q2, sample_data(archae.subset)[, c("MODE.USAGE", "AMENDEMENT.BLE", "TIME", "PLACETTE")])
hillNumber.df=as_tibble(hillNumber.df)
hillNumber.df$PLACETTE = as.factor(hillNumber.df$PLACETTE)
#Réorganiser les données, crer un le vecteur id pour caracteriser les differentes mesures repetées pour les tests statistiques
hillNumber.df = arrange(hillNumber.df, AMENDEMENT.BLE, TIME)
hillNumber.df$id = rep(c(1,2,3), 16) 

## Test anova à trois facteur sur mesure répéte sur les nombres de hill
tab_length=length(colnames(hillNumber.df))-5
Y = colnames(hillNumber.df)[1:tab_length]
cutoff=0.05
for (i in Y) {
  arch.aov = anova_test(
    data = hillNumber.df, dv = i, wid = id,
    within = c(MODE.USAGE, AMENDEMENT.BLE, TIME))
  res=get_anova_table(arch.aov)
  df=data.frame(Effect=res$Effect, p_value=res$p)
  p_sing = filter(df, p_value<=cutoff)
  if (length(p_sing$Effect)>0){
    print(colnames(hillNumber.df[i]))
    print(p_sing)
  }
  
}


#### Représentation graphique des nombres de Hill
arch_boxplot0 = ggplot(data = hillNumber.df , aes(x=TIME, y=archae.q0, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = "Archaea", x=NULL, y="richness")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

arch_boxplot1 = ggplot(data = hillNumber.df , aes(x=TIME, y=archae.q1, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x=NULL, y="Hill shannon")+
  theme_bw()

arch_boxplot2 = ggplot(data = hillNumber.df , aes(x=TIME, y=archae.q2, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x="Date prelevement", y="Hill Simpson")+
  theme_bw()


bact_boxplot0 = ggplot(data = hillNumber.df , aes(x=TIME, y=bacterie.q0, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = "Bacterie", x=NULL, y=NULL)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


bact_boxplot1 = ggplot(data = hillNumber.df , aes(x=TIME, y=bacterie.q1, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x=NULL, y=NULL)+
  theme_bw()

bact_boxplot2 = ggplot(data = hillNumber.df , aes(x=TIME, y=bacterie.q2, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x="Date prelevement", y=NULL)+
  theme_bw()

virus_boxplot0 = ggplot(data = hillNumber.df , aes(x=TIME, y=virus.q0, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = "Virus", x=NULL, y=NULL)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

virus_boxplot1 = ggplot(data = hillNumber.df , aes(x=TIME, y=virus.q1, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x=NULL, y=NULL)+
  theme_bw()

virus_boxplot2 = ggplot(data = hillNumber.df , aes(x=TIME, y=virus.q2, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x="Date prelevement", y=NULL)+
  theme_bw()

fungi_boxplot0 = ggplot(data = hillNumber.df , aes(x=TIME, y=fungi.q0, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = "Champignons", x=NULL, y=NULL)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

fungi_boxplot1 = ggplot(data = hillNumber.df , aes(x=TIME, y=fungi.q1, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x=NULL, y=NULL)+
  #theme(legend.position = "none")+
  theme_bw()

fungi_boxplot2 = ggplot(data = hillNumber.df , aes(x=TIME, y=fungi.q2, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x="Date prelevement", y=NULL)+
  #theme(legend.position = "none")+
  theme_bw()

prot_boxplot0 = ggplot(data = hillNumber.df , aes(x=TIME, y=protiste.q0, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = "Protiste", x=NULL, y=NULL)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

prot_boxplot1 = ggplot(data = hillNumber.df , aes(x=TIME, y=protiste.q1, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x=NULL, y=NULL)+
  theme_bw()

prot_boxplot2 = ggplot(data = hillNumber.df , aes(x=TIME, y=protiste.q2, fill=AMENDEMENT.BLE)) +
  facet_wrap(~MODE.USAGE) +
  geom_boxplot(show.legend = FALSE)+
  labs(title = NULL, x="Date prelevement", y=NULL)+
  theme_bw()

cowplot::plot_grid(arch_boxplot0, bact_boxplot0,fungi_boxplot0, prot_boxplot0, virus_boxplot0,
                   arch_boxplot1, bact_boxplot1,fungi_boxplot1, prot_boxplot1, virus_boxplot1,
                   arch_boxplot2, bact_boxplot2,fungi_boxplot2, prot_boxplot2, virus_boxplot2,
                   #labels = c("A","D","G","B","E","H","C","F","I"),
                   ncol = 5,
                   nrow = 3)
ggsave("./results/Diversites/hillIndice_diversity.png", height = 8, width = 14)



## 3. Analyse de la structuration des populations
####  Analyse en Coordonnées principale PCoA
#fonction "pcoa_analysis" pour faire des analyses en coordoonées principales (PCoA) et test PREMANOVA
times = sample_data(arch.cc)$TIME
pcoa_analysis = function(data, titre) {
  dist = vegan::vegdist(otu_table(data), method = "bray")
  PCOA = pcoa(dist)
  #Veigv = barplot(PCOA$values$Relative_eig[1:10])
  var_ax1 = round(PCOA$values$Relative_eig[1] * 100)
  var_ax2 = round(PCOA$values$Relative_eig[2] * 100)
  sumEigval = round(sum(PCOA$values$Relative_eig)*100)
  
  test.peraov = vegan::adonis2(dist~times,
                               permutations = 999,
                               sqrt.dist = FALSE, 
                               add = FALSE, 
                               by = "terms")
  pv = test.peraov$`Pr(>F)`[1]
  
  dat = as.data.frame(PCOA$vectors) %>%
    select(c(Axis.1, Axis.2))
  dat$MODE.USAGE <- sample_data(data)$MODE.USAGE
  dat$AMENDEMENT.BLE <-sample_data(data)$AMENDEMENT.BLE
  dat$TIME <- sample_data(data)$TIME
  
  dat_centroid <- dat %>%
    group_by(TIME) %>%
    summarise(Axis.1 = mean(Axis.1), Axis.2 = mean(Axis.2)) %>%
    mutate(type = "centroid")
  
  pcoa.plot = dat %>%
    ggplot(aes(x=Axis.1, y=Axis.2, colour=TIME, fill=TIME, group=TIME, ))+
    labs(title=titre,
         x = paste("PCoA 1 [", var_ax1, "%",  "  ", "pvalue", "=", pv, "]", sep = ""),
         y=paste("PCoA 2 [", var_ax2, "%]", " ", sep = ""))+
    geom_polygon(show.legend = FALSE,alpha = 0.1, aes(linetype = NA)) +
    geom_point(shape = 20, show.legend = FALSE) +
    geom_point(data = dat_centroid, size = 3, show.legend = FALSE) +
    theme_bw(base_size = 8)+
    theme(axis.title = element_text(size=8.5), plot.title = element_text(margin=margin(t=40,b=-9), size = 9, face="bold", colour="black"))
  return(pcoa.plot)
}

##### **PCoA sur les données Archae**
set.seed(1234)
archpr.plot = pcoa_analysis(arch.pr, "PR")
archpc.plot = pcoa_analysis(arch.pc, "PC")
archcc.plot = pcoa_analysis(arch.cc, "CC")
archcr.plot = pcoa_analysis(arch.cr, "CR")

cowplot::plot_grid(archcc.plot, archcr.plot,
                   archpc.plot, archpr.plot,
                   ncol = 2,nrow = 2)
ggsave("./results/Diversites/PCoA_archee.png", height = 10, width = 14)

##### PCoA sur les données de Bactéries
set.seed(1234)
bactpr.plot = pcoa_analysis(bact.pr, "PR")
bactpc.plot = pcoa_analysis(bact.pc, "PC")
bactcc.plot = pcoa_analysis(bact.cc, "CC")
bactcr.plot = pcoa_analysis(bact.cr, "CR")

cowplot::plot_grid(bactcc.plot, bactcr.plot,
                   bactpc.plot, bactpr.plot,
                   ncol = 2,nrow = 2)
ggsave("./results/Diversites/PCoA_bacterie.png", height = 10, width = 14)

##### PCoA sur les données des Champignons
set.seed(1234)
fungipr.plot = pcoa_analysis(fungi.pr, "PR")
fungipc.plot = pcoa_analysis(fungi.pc, "PC")
fungicc.plot = pcoa_analysis(fungi.cc, "CC")
fungicr.plot = pcoa_analysis(fungi.cr, "CR")

cowplot::plot_grid(fungicc.plot, fungicr.plot,
                   fungipc.plot, fungipr.plot,
                   ncol = 2,nrow = 2)
ggsave("./results/Diversites/PCoA_fungis.png", height = 10, width = 14)


##### PCoA sur les données des protistes
set.seed(1234)
protpr.plot = pcoa_analysis(prot.pr, "PR")
protpc.plot = pcoa_analysis(prot.pc, "PC")
protcc.plot = pcoa_analysis(prot.cc, "CC")
protcr.plot = pcoa_analysis(prot.cr, "CR")

cowplot::plot_grid(protcc.plot, protcr.plot,
                   protpc.plot, protpr.plot,
                   ncol = 2,nrow = 2)
ggsave("./results/Diversites/PCoA_protiste.png", height = 10, width = 14)

##### PCoA sur les données des virus
set.seed(1234)
viruspr.plot = pcoa_analysis(virus.pr, "PR")
viruspc.plot = pcoa_analysis(virus.pc, "PC")
viruscc.plot = pcoa_analysis(virus.cc, "CC")
viruscr.plot = pcoa_analysis(virus.cr, "CR")

cowplot::plot_grid(viruscc.plot, viruscr.plot,
                   viruspc.plot, viruspr.plot,
                   ncol = 2,nrow = 2)
ggsave("./results/Diversites/PCoA_virus.png", height = 10, width = 14)

cowplot::plot_grid(archcc.plot, bactcc.plot, fungicc.plot, protcc.plot, viruscc.plot,
                   archcr.plot, bactcr.plot, fungicr.plot, protcr.plot, viruscr.plot,
                   archpc.plot, bactpc.plot, fungipc.plot, protpc.plot, viruspc.plot,
                   archpr.plot, bactpr.plot, fungipr.plot, protpr.plot, viruspr.plot,
                   labels=c("Archée","Bactrérie","Champignons","Protiste","Virus"),
                   hjust = -0.75,
                   ncol=5, nrow = 4)
ggsave("./results/Diversites/PCoA_taxo.png", height = 10, width = 14)


####  Analyse NMDS (Non-Métrique Multidimentionnal scaling)
# NMDS a partir de la matrice de bray-Curtis calculer sur les métadonnées de comptage rarefié
metadata.otu = otu_table(metadata.rarefied) %>% t()%>%as.data.frame() 
metadata.nmds = vegan::metaMDS(metadata.otu,  distance="bray", k=2, trymax=100, autotransform=F)
# NMDS à partir des gènes prédits
taxo=c("Kingdom","Phylum","Class","Order","Family","Genus","Species") # Pour les infos taxonomique
data_count = GENE_DATA%>% select(!taxo)%>%t()%>%as.data.frame() # Données de comptage sans la taxonomie transposer
data_taxa = GENE_DATA%>% select(taxo) # Données taxonomiques
nmds.gene = vegan::metaMDS(data_count, distance="bray", k=2, trymax=100, autotransform = FALSE)
# Function "plot_nmds" to plot NMDS
plot_nmds = function(nmdsRes, titre){
  Stress=round(nmdsRes$stress, 3)
  data.scores = as.data.frame(vegan::scores(nmdsRes)$sites)
  data.scores$MODE.USAGE <- sample_data(metadata)$MODE.USAGE
  data.scores$AMENDEMENT.BLE <-sample_data(metadata)$AMENDEMENT.BLE
  data.scores$TIME <- sample_data(metadata)$TIME
  
  data.scores = data.scores%>%mutate(sampleID=paste(str_sub(rownames(.), start=0,end=1), str_sub(rownames(.), start=3, end=7), sep = ""))
  
  data.scores.centroids = merge(data.scores, aggregate(cbind(mean.NMDS1=NMDS1, mean.NMDS2=NMDS2)~sampleID, 
                                                       data.scores,mean), by="sampleID")
  X=min(data.scores.centroids$NMDS1)/1.5
  Y=max(data.scores.centroids$NMDS1)
  
  nmds.plot = data.scores.centroids %>%
    ggplot()+
    geom_point(aes(x=NMDS1, y=NMDS2), size=0.5)+
    geom_point(aes(x=mean.NMDS1, y=mean.NMDS2, shape=AMENDEMENT.BLE, colour=MODE.USAGE), size=10, show.legend=FALSE)+
    scale_shape_manual(values = c(15, 16))+
    #ggsci::scale_color_jco()+
    scale_color_manual(values = c("#fbbb3f", "#01fa03"))+
    labs(title = titre)+
    geom_segment(aes(x=mean.NMDS1, y=mean.NMDS2, xend=NMDS1, yend=NMDS2), colour = "dimgray")+
    geom_text(aes(x=mean.NMDS1, y=mean.NMDS2, label=sampleID), alpha=2,cex=3.5, 
              colour="black", hjust=0.0, vjust=0.0, fontface=2)+
    annotate(geom="text", x=X, y=Y, label=paste("Stress = ", Stress, sep=""), color="red", size=5)+
    theme_bw(base_size = 16)+
    theme(plot.title = element_text(size = 20, face = "bold"), 
          axis.title = element_text(size = 14, face = "bold"))
  return(nmds.plot)
}

metadata.plot = plot_nmds(metadata.nmds, "A")
genes.plot = plot_nmds(nmds.gene, "B")

cowplot::plot_grid(metadata.plot+theme(legend.position = "none",),
                   genes.plot,
                   ncol = 2,nrow = 1)
ggsave("./results/Diversites/NMDS_taxo_gene.png", height = 8, width = 16)


## III. ANALYSE ABONDANCE TAXONOMIQUE
dir.create("./results/Abondances")
deseq2_viz = function(data, rawdata, t2, t1, groupe, temps) {
  sample_data(data)$TIME = as.factor(sample_data(data)$TIME)
  ds = phyloseq_to_deseq2(data, ~TIME)
  ds = DESeq(ds)
  
  cutoff=0.05
  res = results(ds, contrast = c("TIME", t2, t1), alpha = cutoff)
  res = res[order(res$padj, na.last = NA), ]
  res.sig = res[(res$padj<cutoff), ]
  if (length(res.sig$padj)>0){
    res.sig = cbind(as(res.sig, "data.frame"), as(tax_table(rawdata)[rownames(res.sig), ], "matrix"))
    res.sig$TEMPS = temps
    res.sig$GROUPE = groupe
    return(res.sig)
  }
}

arch.cr.T0_T1 = deseq2_viz(arch.cr, archGenus.subset, "T1", "T0", "CR", "T1-T0")
arch.cr.T4_T0 = deseq2_viz(arch.cr, archGenus.subset, "T4", "T0", "CR", "T4_T0")
arch.cr.T6_T0 = deseq2_viz(arch.cr, archGenus.subset, "T6", "T0", "CR", "T6_T0")
arch.pr.T0_T1 = deseq2_viz(arch.pr, archGenus.subset, "T1", "T0", "PR", "T1-T0")
arch.pr.T4_T0 = deseq2_viz(arch.pr, archGenus.subset, "T4", "T0", "PR", "T4-T0")
arch.pr.T6_T0 = deseq2_viz(arch.pr, archGenus.subset, "T6", "T0", "PR", "T6-T0")
arch.resDeseq = rbind(arch.cr.T4_T0, arch.cr.T6_T0, arch.cr.T0_T1,
                      arch.pr.T0_T1, arch.pr.T4_T0, arch.pr.T6_T0)

bact.cr.T0_T1 = deseq2_viz(bact.cr, bactGenus.subset, "T1", "T0", "CR", "T1_T0")
bact.cr.T4_T0 = deseq2_viz(bact.cr, bactGenus.subset, "T4", "T0", "CR", "T4_T0")
bact.cr.T6_T0 = deseq2_viz(bact.cr, bactGenus.subset, "T6", "T0", "CR", "T6_T0")
bact.pr.T0_T1 = deseq2_viz(bact.pr, bactGenus.subset, "T1", "T0", "PR", "T1_T0")
bact.pr.T4_T0 = deseq2_viz(bact.pr, bactGenus.subset, "T4", "T0", "PR", "T4_T0")
bact.pr.T6_T0 = deseq2_viz(bact.pr, bactGenus.subset, "T6", "T0", "PR", "T6_T0")
bact.resDeseq = rbind(bact.cr.T0_T1, bact.cr.T4_T0, bact.cr.T6_T0,
                      bact.pr.T0_T1, bact.pr.T4_T0, bact.pr.T6_T0)

fungi.cr.T0_T1 = deseq2_viz(fungi.cr, fungiGenus.subset, "T1", "T0", "CR", "T1_T0")
fungi.cr.T4_T0 = deseq2_viz(fungi.cr, fungiGenus.subset, "T4", "T0", "CR", "T4_T0")
fungi.cr.T6_T0 = deseq2_viz(fungi.cr, fungiGenus.subset, "T6", "T0", "CR", "T6_T0")
fungi.pr.T0_T1 = deseq2_viz(fungi.pr, fungiGenus.subset, "T1", "T0", "PR", "T1_T0")
fungi.pr.T4_T0 = deseq2_viz(fungi.pr, fungiGenus.subset, "T4", "T0", "PR", "T4_T0")
fungi.pr.T6_T0 = deseq2_viz(fungi.pr, fungiGenus.subset, "T6", "T0", "PR", "T6_T0")
fungi.resDeseq = rbind(fungi.cr.T0_T1, fungi.cr.T4_T0, fungi.cr.T6_T0, 
                       fungi.pr.T0_T1,fungi.pr.T4_T0, fungi.pr.T6_T0)


prot.cr.T0_T1 = deseq2_viz(prot.cr, protisteGenus.subset, "T1", "T0", "CR", "T1_T0")
prot.cr.T4_T0 = deseq2_viz(prot.cr, protisteGenus.subset, "T4", "T0", "CR", "T4_T0")
prot.cr.T6_T0 = deseq2_viz(prot.cr, protisteGenus.subset, "T6", "T0", "CR", "T6_T0")
prot.pr.T0_T1 = deseq2_viz(prot.pr, protisteGenus.subset, "T1", "T0", "PR", "T1_T0")
prot.pr.T4_T0 = deseq2_viz(prot.pr, protisteGenus.subset, "T4", "T0", "PR", "T4_T0")
prot.pr.T6_T0 = deseq2_viz(prot.pr, protisteGenus.subset, "T6", "T0", "PR", "T6_T0")
prot.resDeseq = rbind(prot.cr.T0_T1, prot.cr.T4_T0, prot.cr.T6_T0, 
                      prot.pr.T0_T1, prot.pr.T4_T0, prot.pr.T6_T0)


virus.cr.T0_T1 = deseq2_viz(virus.cr, virusGenus.subset, "T1", "T0", "CR", "T1_T0")
virus.cr.T4_T0 = deseq2_viz(virus.cr, virusGenus.subset, "T4", "T0", "CR", "T4_T0")
virus.cr.T6_T0 = deseq2_viz(virus.cr, virusGenus.subset, "T6", "T0", "CR", "T6-T0")
virus.pr.T0_T1 = deseq2_viz(virus.pr, virusGenus.subset, "T1", "T0", "PR", "T1_T0")
virus.pr.T4_T0 = deseq2_viz(virus.pr, virusGenus.subset, "T4", "T0", "PR", "T4-T0")
virus.pr.T6_T0 = deseq2_viz(virus.pr, virusGenus.subset, "T6", "T0", "PR", "T6_T0")
virus.resDeseq = rbind(virus.cr.T0_T1, virus.cr.T4_T0, virus.cr.T6_T0, 
                       virus.pr.T0_T1, virus.pr.T4_T0, virus.pr.T6_T0)

couleur = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#0213FF","#fdcdac","#cbd5e8","#f4cae4","#e6f5c9",
                     "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3",
                     "#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6",
                     "#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e", "#fc8d62","#000000","#8da0cb","#e78ac3","#a6d854",
                     "2f1b0c")

plot_deseq2 = function(deseqdata, PHYLUM, titres) {
  data = arrange(deseqdata, PHYLUM)
  data = mutate(data, Family = as.factor(Family))
  data = mutate(data, Family = factor(Family, levels = unique(data$Family)))
  if ("Bacteria" %in% deseqdata[1,7]){
    ploting = ggplot(data, aes(x = Family, y = log2FoldChange, color = Phylum)) + 
      #geom_point()+
      geom_jitter(size = 0.95) +
      geom_hline(yintercept = 0, linetype="solid", color="red", size=0.25)+
      labs(title = titres)+
      theme_bw() + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major  = element_blank(),
            legend.text = element_text(size = 14), legend.position = "bottom", plot.title = element_text(size=16, face = "bold")) + 
      #facet_wrap(~TEMPS)+
      scale_color_manual(values = couleur)+
      facet_grid(rows = vars(GROUPE), cols = vars(TEMPS), scales = "free")
    return(ploting)
  }
  else {ploting = ggplot(data, aes(x = Family, y = log2FoldChange, color = PHYLUM)) + 
    geom_jitter(size =1.75) +
    geom_hline(yintercept = 0, linetype="solid", color="red", size=0.25)+
    theme_bw() + 
    labs(title = titres)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major  = element_blank(),
          legend.text = element_text(size = 13), legend.position = "bottom", plot.title = element_text(size=16, face = "bold")) + 
    scale_color_manual(values = couleur)+
    facet_grid(rows = vars(GROUPE), cols = vars(TEMPS), scales = "free")
  return(ploting)
  }
} 

##### Archae
archPlot.deseq = plot_deseq2(arch.resDeseq, arch.resDeseq$Phylum, "Archées")
archPlot.deseq
ggsave("./results/Abondances/resDeseq2_archee.png", height = 8, width = 16)

##### Bacterie
bactPlot.deseq = plot_deseq2(bact.resDeseq, bact.resDeseq$Phylum, "Bactéries")
bactPlot.deseq
ggsave("./results/Abondances/resDeseq2_bacterrie.pdf", height = 10, width = 16)

##### Champignons
fungiPlot.deseq = plot_deseq2(fungi.resDeseq, fungi.resDeseq$Phylum, "Champignons")
fungiPlot.deseq
ggsave("./results/Abondances/resDeseq2_fungi.pdf", height = 8, width = 16)

##### Protiste
protPlot.deseq = plot_deseq2(prot.resDeseq, prot.resDeseq$Phylum, "Protistes")
protPlot.deseq
ggsave("./results/Abondances/resDeseq2_protiste.pdf", height = 8, width = 16)

##### **virus**
virusPlot.deseq = plot_deseq2(virus.resDeseq, virus.resDeseq$Class, "Virus")
virusPlot.deseq
ggsave("./results/Abondances/resDeseq2_virus.pdf", height = 8, width = 16)

cowplot::plot_grid(archPlot.deseq,fungiPlot.deseq,
                   protPlot.deseq,virusPlot.deseq,
                   ncol = 2, nrow = 2)
ggsave("./results/Abondances/resDeseq2_Arch_fungi_prot_virus.pdf", height = 8, width = 16)


#### Exploration des résultats DESeq2
#Fonction resDeseq2_parsed pour determiner les differences entre les genres selon les temps, le mode d'usage et dans chaque condition
resDeseq2_parsed = function(data, groupe){
  if ("Viruses" %in% data[1,7]){
    data.parsed = data %>%
      filter(GROUPE == groupe)%>%
      group_by(Class, Genus, TEMPS)%>%
      summarise(count=n(), L2FC=log2FoldChange)%>%
      filter(Genus!="NA")%>% # Supprimer les Genres non affiliées
      pivot_wider(id_cols = Genus, names_from = TEMPS, values_from =c(count, L2FC)) # Pour Répartir suivant le temps en fonction du l2FC et du count
    colnames(data.parsed)=c("Genus",paste(groupe,colnames(data.parsed)[-1], sep = "_"))
    #data.parsed2=cbind(data.parsed, Superingdom=rep(data[1,7], nrow(data.parsed)))
    #data.cr[is.na(data.cr)]=0
    return(data.parsed)
  }
  else{
    data.parsed = data %>%
      filter(GROUPE == groupe)%>%
      group_by(Phylum, Genus, TEMPS)%>%
      summarise(count=n(), L2FC=log2FoldChange)%>%
      filter(Genus!="NA")%>% # Supprimer les Genres non affiliées
      pivot_wider(id_cols = Genus, names_from = TEMPS, values_from =c(count, L2FC)) # Pour Répartir suivant le temps en fonction du l2FC et du count
    colnames(data.parsed)=c("Genus",paste(groupe,colnames(data.parsed)[-1], sep = "_"))
    #data.parsed=cbind(data.parsed, Superingdom=rep(data[1,7], nrow(data.parsed)))
    #data.cr[is.na(data.cr)]=0
    return(data.parsed)
  }
  
}

#Execution de la fonction parse3_resDesq2
# il n'ya pas d'archae PR
archGenus_table.resdeseq2 = resDeseq2_parsed(arch.resDeseq, "CR")
archGenus_table.resdeseq2[is.na(archGenus_table.resdeseq2)]=0
write.table(archGenus_table.resdeseq2, file="./results/Abondances/2022_04_resultats_archaea_analysetaxonomique_deseq2_L2FCsorted.csv", sep=";", row.names=FALSE, col.names = TRUE)

bactCR_table.resdeseq2 = resDeseq2_parsed(bact.resDeseq, "CR")
bactPR_table.resdeseq2 = resDeseq2_parsed(bact.resDeseq, "PR")
bactGenus_table.resdeseq2=full_join(bactCR_table.resdeseq2, bactPR_table.resdeseq2)
bactGenus_table.resdeseq2[is.na(bactGenus_table.resdeseq2)]=0
write.table(bactGenus_table.resdeseq2, file="./results/Abondances/2022_04_resultats_bacterie_analysetaxonomique_deseq2_L2FCsorted.csv", sep=";", row.names=FALSE, col.names = TRUE)

fungiCR_table.resdeseq2 = resDeseq2_parsed(fungi.resDeseq, "CR")
fungiPR_table.resdeseq2 = resDeseq2_parsed(fungi.resDeseq, "PR")
fungiGenus_table.resdeseq2 = full_join(fungiCR_table.resdeseq2, fungiPR_table.resdeseq2)
fungiGenus_table.resdeseq2[is.na(fungiGenus_table.resdeseq2)]=0
write.table(fungiGenus_table.resdeseq2, file="./results/Abondances/2022_04_resultats_fungi_analysetaxonomique_deseq2_L2FCsorted.csv", sep=";", row.names=FALSE, col.names = TRUE)

protCR_table.resdeseq2 = resDeseq2_parsed(prot.resDeseq, "CR")
protPR_table.resdeseq2 = resDeseq2_parsed(prot.resDeseq, "PR")
protGenus_table.resdeseq2 = full_join(protCR_table.resdeseq2, protPR_table.resdeseq2)
protGenus_table.resdeseq2[is.na(protGenus_table.resdeseq2)]=0
write.table(protGenus_table.resdeseq2, file="./results/Abondances/2022_04_resultats_protiste_analysetaxonomique_deseq2_L2FCsorted.csv", sep=";", row.names=FALSE, col.names = TRUE)

virusCR_table.resdeseq2 = resDeseq2_parsed(virus.resDeseq, "CR")
virusPR_table.resdeseq2 = resDeseq2_parsed(virus.resDeseq, "PR")
virusGenus_table.resdeseq2 = full_join(virusCR_table.resdeseq2, virusPR_table.resdeseq2)
virusGenus_table.resdeseq2[is.na(virusGenus_table.resdeseq2)]=0
write.table(virusGenus_table.resdeseq2, file="./results/Abondances/2022_04_resultats_virus_analysetaxonomique_deseq2_L2FCsorted.csv", sep=";", row.names=FALSE, col.names = TRUE)


# Définir un fonction pour la classification des genres valable pour les Bactéries, Champignons et les Protistes
select_taxon=function(data, taxons){
  tax_select = data %>%
    mutate(Statut=case_when(
      CR_count_T1_T0+CR_count_T4_T0+CR_count_T6_T0+PR_count_T1_T0+PR_count_T4_T0+PR_count_T6_T0==6 ~ "common_crop_grassland",
      (CR_count_T1_T0+CR_count_T4_T0+CR_count_T6_T0==3)&(PR_count_T1_T0+PR_count_T4_T0+PR_count_T6_T0==0) ~ "only_crop",
      (PR_count_T1_T0+PR_count_T4_T0+PR_count_T6_T0==3)&(CR_count_T1_T0+CR_count_T4_T0+CR_count_T6_T0==0) ~ "only_grassland",
      TRUE ~ "Statut_Not_defined"),
      Copio_Oligo=case_when(
        (CR_count_T1_T0+PR_count_T1_T0==2)&(CR_count_T4_T0+CR_count_T6_T0+PR_count_T4_T0+PR_count_T6_T0==0)&CR_L2FC_T1_T0>0&PR_L2FC_T1_T0>0&CR_L2FC_T1_T0>CR_L2FC_T6_T0&PR_L2FC_T1_T0>PR_L2FC_T6_T0 ~ "common_copio",
        (CR_count_T1_T0==1)&(CR_count_T4_T0+CR_count_T6_T0+PR_count_T1_T0+PR_count_T4_T0+PR_count_T6_T0==0)&CR_L2FC_T1_T0>0&PR_L2FC_T1_T0==0&CR_L2FC_T1_T0>CR_L2FC_T6_T0 ~ "only_copio_crop",
        (PR_count_T1_T0==1)&(CR_count_T1_T0+CR_count_T4_T0+CR_count_T6_T0+PR_count_T4_T0+PR_count_T6_T0==0)&CR_L2FC_T1_T0==0&PR_L2FC_T1_T0>0&PR_L2FC_T1_T0>PR_L2FC_T6_T0 ~ "only_copio_grassland",
        (CR_count_T6_T0+PR_count_T6_T0==2)&(CR_count_T1_T0+CR_count_T4_T0+PR_count_T1_T0+PR_count_T4_T0==0)&CR_L2FC_T6_T0>0&PR_L2FC_T6_T0>0&CR_L2FC_T6_T0>CR_L2FC_T1_T0&PR_L2FC_T6_T0>PR_L2FC_T1_T0 ~ "common_oligo",
        (CR_count_T6_T0==1)&(CR_count_T1_T0+CR_count_T4_T0+PR_count_T1_T0+PR_count_T4_T0+PR_count_T6_T0==0)&CR_L2FC_T6_T0>0&PR_L2FC_T6_T0==0&CR_L2FC_T6_T0>CR_L2FC_T1_T0 ~ "only_oligo_crop",
        (PR_count_T6_T0==1)&(CR_count_T1_T0+CR_count_T4_T0+CR_count_T6_T0+PR_count_T1_T0+PR_count_T4_T0==0)&CR_L2FC_T6_T0==0&PR_L2FC_T6_T0>0&PR_L2FC_T6_T0>PR_L2FC_T1_T0 ~ "only_oligo_grassland",
        (CR_count_T4_T0+PR_count_T4_T0==2)&(CR_count_T1_T0+CR_count_T6_T0+PR_count_T1_T0+PR_count_T6_T0==0)&CR_L2FC_T4_T0>0&PR_L2FC_T4_T0>0 ~ "uncopio_oligo",
        TRUE ~ "Copio_Oligo_Not_defined"),
      Dynamic=case_when(
        (CR_count_T1_T0+PR_count_T1_T0+CR_count_T4_T0+PR_count_T4_T0==4)&(CR_count_T6_T0+PR_count_T6_T0==0) ~ "T1_T4_common",
        (CR_count_T1_T0+CR_count_T4_T0==2)&(CR_count_T6_T0+PR_count_T6_T0+PR_count_T1_T0+PR_count_T4_T0==0) ~ "T1_T4_crops",
        (PR_count_T1_T0+PR_count_T4_T0==2)&(CR_count_T6_T0+PR_count_T6_T0+CR_count_T1_T0+CR_count_T4_T0==0) ~ "T1_T4_grasslands",
        TRUE ~ "T1_T4_Not_defined"))
  colname=c("STATUT",taxons)
  tax.table=table(tax_select$Statut)%>% as.data.frame() %>% `colnames<-`(colname) %>%
    rbind(table(tax_select$Copio_Oligo)%>% as.data.frame() %>% `colnames<-`(colname),table(tax_select$Dynamic)%>% as.data.frame() %>% `colnames<-`(colname))
  return(tax.table)
}


# fonction select_taxon2 pour le traitement des virus
# Pas de Dynamic ni de "common_oligo" et "only_oligo_crop" par ce qu'il ya pas les temps CR_T6_T0 et PR_T4_T0,
# Pas necessaire de prendre en compte le log2FoldChange pour déterminer les copiotrophes et les oligotrophes car la plus grand majorité s'exprime tous au temps T1 
select_taxon2=function(data, taxons){
  tax_select = data %>%
    mutate(Statut=case_when(
      CR_count_T1_T0+CR_count_T4_T0+PR_count_T1_T0+PR_count_T6_T0==4 ~ "common_crop_grassland",
      (CR_count_T1_T0+CR_count_T4_T0==2)&(PR_count_T1_T0+PR_count_T6_T0==0) ~ "only_crop",
      (PR_count_T1_T0+PR_count_T6_T0==2)&(CR_count_T1_T0+CR_count_T4_T0==0) ~ "only_grassland",
      (CR_count_T1_T0+PR_count_T1_T0==2)&(CR_count_T4_T0+PR_count_T6_T0==0) ~ "common_copio",
      (CR_count_T1_T0==1)&(CR_count_T4_T0+PR_count_T1_T0+PR_count_T6_T0==0) ~ "only_copio_crop",
      (PR_count_T1_T0==1)&(CR_count_T1_T0+CR_count_T4_T0+PR_count_T6_T0==0) ~ "only_copio_grassland",
      #(CR_T6_T0+PR_T6_T0==2)&(CR_T1_T0+CR_T4_T0+PR_T1_T0+PR_T4_T0==0) ~ "common_oligo",
      #(CR_T6_T0==1)&(CR_T1_T0+CR_T4_T0+PR_T1_T0+PR_T4_T0+PR_T6_T0==0) ~ "only_oligo_crop",
      (PR_count_T6_T0==1)&(CR_count_T1_T0+CR_count_T4_T0+PR_count_T1_T0==0) ~ "only_oligo_grassland",
      TRUE ~ "Statut_Not_defined"))
  #Dynamic=case_when(
  #(CR_T1_T0+PR_T1_T0+CR_T4_T0+PR_T4_T0==4)&(CR_T6_T0+PR_T6_T0==0) ~ "T1_T4_common",
  #(CR_T1_T0+CR_T4_T0==2)&(CR_T6_T0+PR_T6_T0+PR_T1_T0+PR_T4_T0==0) ~ "T1_T4_crops",
  #(PR_T1_T0+PR_T4_T0==2)&(CR_T6_T0+PR_T6_T0+CR_T1_T0+CR_T4_T0==0) ~ "T1_T4_grasslands",
  #TRUE ~ "T1_T4_Not_defined"))
  colname=c("STATUT",taxons)
  tax.table=table(tax_select$Statut)%>% as.data.frame() %>% `colnames<-`(colname)
  return(tax.table)
}

# Traitement des archaea ou nous avons que les temps CR_T4_T0 et CR_T6_T0 
arch_select = archGenus_table.resdeseq2 %>%
  mutate(Statut=case_when(
    CR_count_T4_T0+CR_count_T6_T0==2 ~ "common_crop_grassland",
    (CR_count_T6_T0==1)&(CR_count_T4_T0==0) ~ "only_oligo_crop",
    TRUE ~ "Statut_Not_defined"))
colname=c("STATUT","ARCHAEA")
archae_select=table(arch_select$Statut)%>% as.data.frame() %>% `colnames<-`(colname)

bact_select=select_taxon(bactGenus_table.resdeseq2, "BACTERIE")
fungi_select=select_taxon(fungiGenus_table.resdeseq2, "CHAMPIGNONS")
prot_select=select_taxon(protGenus_table.resdeseq2, "PROTISTE")
virus_select=select_taxon2(virusGenus_table.resdeseq2, "VIRUS")
taxon_select = full_join(bact_select,fungi_select)%>%full_join(prot_select)%>%full_join(virus_select)%>%full_join(archae_select)%>%as.data.frame()
taxon_select[is.na(taxon_select)]=0
taxon_select = column_to_rownames(taxon_select, "STATUT")%>%t()

#Déterminer le nombre de total de genre et les differences entre Culture et prairie

Nombre_Total_Genre=c(length(rownames(bactGenus_table.resdeseq2)),
                     length(rownames(fungiGenus_table.resdeseq2)),
                     length(rownames(protGenus_table.resdeseq2)),
                     length(rownames(virusGenus_table.resdeseq2)),
                     length(rownames(archGenus_table.resdeseq2)))
Nombre_Genres_Different=c(length(rownames(bactGenus_table.resdeseq2))-length(rownames(bactCR_table.resdeseq2)),
                          length(rownames(fungiGenus_table.resdeseq2))-length(rownames(fungiCR_table.resdeseq2)), 
                          length(rownames(protGenus_table.resdeseq2))-length(rownames(protCR_table.resdeseq2)), 
                          length(rownames(virusGenus_table.resdeseq2))-length(rownames(virusCR_table.resdeseq2)), 
                          length(rownames(archGenus_table.resdeseq2)))
only_grassland=rep(0, length(taxon_select)-1) # Pour rajouter les info spécifique a prairie (only_grassland) égale à 0 dans tous les phylum
taxon_select=cbind(taxon_select, only_grassland, Nombre_Total_Genre, Nombre_Genres_Different)%>%t()
taxon_select=rownames_to_column(as.data.frame(taxon_select), "STATUT")


#### Tableau des Résultats de sélection Phylum statut
select.satut=filter(taxon_select, STATUT %in% c("common_crop_grassland","only_grassland","only_crop"))
select.nbre_genre=filter(taxon_select, STATUT %in% c("Statut_Not_defined","Nombre_Total_Genre","Nombre_Genres_Different"))
select.stat=rbind(select.satut,select.nbre_genre)
write.table(select.stat, file="./results/Abondances/2022_04_PhylumStatut.csv", sep=";", row.names=FALSE, col.names = TRUE)
#### Tableau des Résultats de sélection copio et oligotrope
select.copio=filter(taxon_select, STATUT %in% c("common_copio","only_copio_crop","only_copio_grassland"))
select.oligo=filter(taxon_select, STATUT %in% c("common_oligo","only_oligo_cr op","only_oligo_grassland","uncopio_oligo"))
select_not_def=filter(taxon_select, STATUT %in% ("Copio_Oligo_Not_defined"))
select.copio.oligo=rbind(select.copio,select.oligo,select_not_def)
write.table(select.copio.oligo, file="./results/Abondances/2022_04_PhylumCopio_Oligo.csv", sep=";", row.names=FALSE, col.names = TRUE)
#### Tableau des Résultats de sélection non copio-oligotrophe
select.T1_T4=filter(taxon_select, STATUT %in% c("T1_T4_common","T1_T4_crops","T1_T4_grasslands","Nombre_Total_Genre"))
write.table(select.T1_T4, file="./results/Abondances/2022_04_Phylum_NoCopio_Oligo.csv", sep=";", row.names=FALSE, col.names = TRUE)



