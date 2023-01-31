
#####################################################################
##   "FonctioMic - préparation données pour analyse des Réseaux"   ##
##            __ Samuel Desquiet / Alphonse THIAW __               ##
#####################################################################

# I. Objectifs
#Ce script a pour objectif de mettre en forme les fichiers de matrice OTUs issue de Reclustor pour 
#les bactéries et champignons.


#Chargement des packages
library(tidyverse)
library(magrittr)
library(Hmisc)
library(data.table)
library(reshape)
library(reshape2)
library(GGally)
library(RColorBrewer)
library(ggplot2)
library(igraph)
library(ggraph)
library(graphlayouts)
library(phyloseq)


# Définir le répertoire de travail
#DOSSIER <- "/home/alphonse/Abt_docs/abtStage_M2/Stage_M2"
DOSSIER <- "/home/abthiaw/UCA/test"
#DOSSIER <- "/home/abthiaw/work"
setwd(DOSSIER)

#########################################################################################


# I. Chargement et mise en forme des données

metadata <- read_rds("./data/sols_physeq.rds")
load("./data/2022_05_Fonctiomic_taxonomicData.RData")
#metadata.phylum <- tax_glom(metadata, "Phylum", NArm = FALSE)
#Raréfier les données 
#metadata.phylum.raref <- rarefy_even_depth(metadata.phylum, rngseed = 1, sample.size = 0.95*min(sample_sums(metadata.phylum)),replace=FALSE)
# Remarque : le rarefaction n'apporte pas de changement sur le jeu de données

### CHARGEMENT DE LA TABLE D'OTU ET RABOUTAGE DES INFOS DE KINGDOM
taxa_table <- tax_table(metadata.phylum)%>%as.data.frame()%>%
  mutate(rep = rep(1:nrow(.))) %>% 
  mutate(Superkingdom2 = paste(Superkingdom, rep, sep = "_"))

matOTUs <- t(otu_table(metadata.phylum)) %>% `colnames<-` (taxa_table$Superkingdom2) %>% as.data.frame(.)


#### Bacteria ####
bacterie <- matOTUs %>% select(grep(colnames(.), pattern = "Bacteria"))
### Eukaryota ###
eukaryote <- matOTUs %>% select(grep(colnames(.), pattern = "Eukaryota"))
### Viruses ###
virus = matOTUs %>% select(grep(colnames(.), pattern = "Viruses")) 
### Archaea ###
archae = matOTUs %>% select(grep(colnames(.), pattern = "Archaea"))


### MODALITES ### 
modalites <- data.frame(row.names(matOTUs)) %>% 
  `colnames<-` ("code_ech") %>% 
  mutate(modalite = paste(str_sub(code_ech, start = 0, end = 1), str_sub(code_ech, start = 4, end = 4), sep = "_")) %>% 
  column_to_rownames("code_ech")

#########################################################################################

# III. Script Réseau 

### CHOIX MODALITE
nom_modalite = "modalite"
groupes_temp <- modalites
groupes_temp$modalites <- as.factor(as.character(groupes_temp[, colnames(groupes_temp) == nom_modalite]))
summary(groupes_temp$modalites)
# choisir le nombre d'echantillons par modalite (prendre un chiffre plus petit que le nombre d'echantillons 
# dans la modalite la moins representee pour pouvoir faire des repetitions)
NbEch = 12
#NbEch = 9
# choisir le nombre de repetitions de reseau
NbRep = 1
AbMinOTUs = 5 # dans combien d'echantillons je dois retrouver mon OTU, a fixer en fonction du nombre d'echantillons, prendre 50%
#AbMinOTUs = 1
seuil.p = 0.01 # definir le seuil du pvalue
#seuil.p = 0.01

setwd(DOSSIER)
results <- NULL
for (m in 1 : NbRep){
  
  temp <- NULL
  for (k in 1 : nlevels(groupes_temp$modalites)){
    temp <- c(temp, sample(row.names(groupes_temp)[groupes_temp$modalites == levels(groupes_temp$modalites) [k]], NbEch, replace = FALSE))
    #temp <- c(temp, row.names(groupes_temp)[groupes_temp$modalites == levels(groupes_temp$modalites) [k]])
  }
  
  groupes <- groupes_temp [temp,]
  groupes <- groupes[order(row.names(groupes)),]
  # selection des echantillons dans les OTUS
  OTU_bact <- bacterie[row.names(bacterie) %in% row.names(groupes), ]
  OTU_bact <- OTU_bact[order(row.names(OTU_bact)),]
  row.names(groupes) == row.names(OTU_bact)
  OTU_euk <- eukaryote[row.names(eukaryote) %in% row.names(groupes), ]
  OTU_euk <- OTU_euk[order(row.names(OTU_euk)),]
  row.names(groupes) == row.names(OTU_euk)
  OTU_virus <- virus[row.names(virus) %in% row.names(groupes), ]
  OTU_virus<- OTU_virus[order(row.names(OTU_virus)),]
  row.names(groupes) == row.names(OTU_virus)
  OTU_arch <- archae[row.names(archae) %in% row.names(groupes), ]
  OTU_arch<- OTU_arch[order(row.names(OTU_arch)),]
  row.names(groupes) == row.names(OTU_arch)
  
  for (i in 1 : length (levels(groupes$modalites))){
    dir.create(paste(DOSSIER, "Reseaux_Metagenomique2_Phylum", sep = "/"))
    dir.create(paste(DOSSIER, "Reseaux_Metagenomique2_Phylum", nom_modalite, sep = "/"))
    dir.create(paste(DOSSIER, "Reseaux_Metagenomique2_Phylum", nom_modalite,  paste(seuil.p, "_", NbEch, "ech_", NbRep, "rep_", AbMinOTUs, "GENRESmin_2021_2022", sep = ""), sep = "/"))
    dir.create(paste(DOSSIER, "Reseaux_Metagenomique2_Phylum", nom_modalite,  paste(seuil.p, "_", NbEch, "ech_", NbRep, "rep_", AbMinOTUs, "GENRESmin_2021_2022", sep = ""), levels(groupes$modalites) [i], sep = "/"))
    dir.create(paste(DOSSIER, "Reseaux_Metagenomique2_Phylum", nom_modalite,  paste(seuil.p, "_", NbEch, "ech_", NbRep, "rep_", AbMinOTUs, "GENRESmin_2021_2022", sep = ""), levels(groupes$modalites) [i], "rep",sep = "/"))
    dir.create(paste(DOSSIER, "Reseaux_Metagenomique2_Phylum", nom_modalite,  paste(seuil.p, "_", NbEch, "ech_", NbRep, "rep_", AbMinOTUs, "GENRESmin_2021_2022", sep = ""), levels(groupes$modalites) [i], "rep", m, sep = "/"))
    setwd(paste(DOSSIER, "Reseaux_Metagenomique2_Phylum", nom_modalite,  paste(seuil.p, "_", NbEch, "ech_", NbRep, "rep_", AbMinOTUs, "GENRESmin_2021_2022", sep = ""), levels(groupes$modalites) [i], "rep", m, sep = "/"))
    #dir.create("Data_used")
    #dir.create ("Resultats")
    
    # setwd(paste(DOSSIER, "Reseaux_16Set18S", levels(groupes$modalites) [i], "rep", m, "Data_used", sep = "/"))
    dat_bact <- OTU_bact[groupes$modalites == levels(groupes$modalites)[i], ]
    #dat_bact <- dat_bact[, colSums(dat_bact>0) > AbMinOTUs]
    dat_bact <- dat_bact[, colSums(dat_bact) > AbMinOTUs]
    NbCol_bact <- ncol(dat_bact)
    dat_euk <- OTU_euk[groupes$modalites == levels(groupes$modalites)[i], ]
    #dat_euk <- dat_euk[, colSums(dat_euk>0) > AbMinOTUs]
    dat_euk <- dat_euk[, colSums(dat_euk) > AbMinOTUs]
    NbCol_euk <- ncol(dat_euk)
    dat_virus <- OTU_virus[groupes$modalites == levels(groupes$modalites)[i], ]
    #dat_virus <- dat_virus[, colSums(dat_virus>0) > AbMinOTUs]
    dat_virus <- dat_virus[, colSums(dat_virus) > AbMinOTUs]
    NbCol_virus <- ncol(dat_virus)
    dat_arch <- OTU_arch[groupes$modalites == levels(groupes$modalites)[i], ]
    #dat_arch<- dat_arch[, colSums(dat_arch>0) > AbMinOTUs]
    dat_arch<- dat_arch[, colSums(dat_arch) > AbMinOTUs]
    NbCol_arch <- ncol(dat_arch)
    
    dat <- cbind(dat_bact, dat_euk, dat_virus, dat_arch)
    
    write.csv2(dat, file = paste ("Reseaux_Metagenomique2_Phylum", levels(groupes$modalites)[i], ".csv", sep = ""), row.names = TRUE)
    # RV = rcorr(as.matrix(rbind(dat, dat)), type = "spearman") # calcul des correlations en doublant les replicats
    RV = rcorr(as.matrix(dat), type = "spearman") # calcul des correlations sans doubler les replicats
    #RV$r <- ifelse(RV$r == "NaN", 0, RV$r) # RV$r est la table reprenant les coefficients de correlation
    #RV$r <- diag.remove(dat = RV$r, remove.val = 0) # remplacer les valeurs de la diagonale par 0
    RV$P.adj<-p.adjust(RV$P, method=c("fdr"))
    RV$P.adj[is.na(RV$P.adj)]<-1 # remplacement de la pvalue avec des NAs par 1
    RV$P.adj<-matrix(RV$P.adj, ncol=ncol(RV$P))
    ## Préparation des tableaux de résultats (correlation positive et negative et sens de la correlation)
    cor.sig<-matrix(nrow=nrow(RV$r), ncol=ncol(RV$r))
    colnames(cor.sig)<-colnames(RV$r)
    rownames(cor.sig)<-colnames(RV$r)
    cor.pos<-matrix(nrow=nrow(RV$r), ncol=ncol(RV$r))
    colnames(cor.pos)<-colnames(RV$r)
    rownames(cor.pos)<-colnames(RV$r)
    cor.neg<-matrix(nrow=nrow(RV$r), ncol=ncol(RV$r))
    colnames(cor.neg)<-colnames(RV$r)
    rownames(cor.neg)<-colnames(RV$r)
    weight<-matrix(nrow=nrow(RV$r), ncol=ncol(RV$r))
    colnames(weight)<-colnames(RV$r)
    rownames(weight)<-colnames(RV$r)
    dir<-matrix(nrow=nrow(RV$r), ncol=ncol(RV$r))
    colnames(dir)<-colnames(RV$r)
    rownames(dir)<-colnames(RV$r)
    
    ## Selection des liens significatifs a partir du r 
    
    for(k in 1:ncol(RV$r)){     # LOOP 3
      for(j in 1:nrow(RV$r)){     # LOOP 4
        if (as.numeric(RV$P.adj[k,j])<seuil.p){     # LOOP 5 ### pvalue fixee a 0,06 car fdr est deja tres stringent
          cor.sig[k,j]<-as.numeric(RV$r[k,j])
          weight[k,j]<-abs(as.numeric(RV$r[k,j]))
          if (as.numeric(RV$r[k,j])>0){     # LOOP 6
            cor.pos[k,j]<-abs(as.numeric(RV$r[k,j]))
            cor.neg[k,j]<-0
            dir[k,j]<-c("positive")
          }
          else {
            cor.neg[k,j]<-abs(as.numeric(RV$r[k,j]))
            cor.pos[k,j]<-0
            dir[k,j]<-c("negative")
          }
        }
        else {
          cor.sig[k,j]<-0
          weight[k,j]<-0
          dir[k,j]<-c("NS")
        }
      }
    }
    
  
    ## Elaborer la table des liens
    
    dir2 <- as.data.frame(cbind(row.names(dir),dir))
    colnames(dir2) [1] <- "Input"
    dir.long <- melt(dir2, id = "Input") # transformation de la matrice en une colonne
    weight2 <- as.data.frame(cbind(row.names(weight),weight))
    colnames(weight2) [1] <- "Input"
    weight.long <- melt(weight2, id = "Input") # transformation de la matrice en une colonne
    links <- cbind(dir.long, weight.long[,3]) 
    colnames(links) <- c("input", "output", "direction", "weight")
    links <- links[links$direction != "NS", ]
    links <- links[links$input != links$output, ]
    links <- links[!duplicated(t(apply(links[,1:2], 1, sort))), ]
    
    # recuperation des regnes des noeuds input et output 
    links <- links %>% 
      mutate(input_regne = colsplit(input, pattern = "_", names = c("regne", "numero"))[,1],
             output_regne = colsplit(output, pattern = "_", names = c("regne", "numero"))[,1],
             regnes_lien = paste(input_regne, output_regne, sep = "_")) %>% 
      mutate(regnes_lien = case_when(regnes_lien == "Bacteria_Archea" ~ "Archaea_Bacteria",
                                     regnes_lien == "Eukaryota_Archaea" ~ "Archaea_Eukaryota",
                                     regnes_lien == "Viruses_Archaea" ~ "Archaea_Viruses",
                                     regnes_lien == "Eukaryota_Bacteria" ~ "Bacteria_Eukaryota",
                                     regnes_lien == "Viruses_Bacteria" ~ "Bacteria_Viruses",
                                     regnes_lien == "Viruses_Eukaryota" ~ "Eukaryota_Viruses",
                                     TRUE ~ regnes_lien))%>%
      select(c("input","output","weight", everything()))
    
    links$weight=as.numeric(links$weight)
    ## Élaborer la  table des noeuds
    nodes <- taxa_table%>%select(c( "Superkingdom2", "Superkingdom","Kingdom","Phylum"))
    
    ## Implémentation du réseau
    #N <- as.network(as.matrix(cor.sig), directed=FALSE)
    Net <- graph_from_data_frame(links, nodes, directed=F)
    ## Calcul des indices de réseaux
    Nb_echant <- nrow(dat) #Nombre d'echantillon
    #N0 <- ncol(dat) #Nbre de noeuds
    
    #N1 <- N$gal$n #Nbre de noeud calculer
    Nb_noeuds <- length(V(Net)) #Nbre de noeud
    Nblien_Pos <- sum(dir=="positive")/2 #Nbre de lien positive
    Nblien_Neg <- sum(dir=="negative")/2 #Nbre de lien Négative
    Nb_liens <- Nblien_Pos+Nblien_Neg #Nbre de lien
    ifelse(Nblien_Pos/Nblien_Neg!= "NaN", Ratio_lien_Pos.Neg <- Nblien_Pos/Nblien_Neg, Ratio_lien_Pos.Neg <- "NA") #Rapport lien positive/negative
    #connectance <- gden(N)
    Densite <- edge_density(Net, loops = F)
    #transitivity <- gtrans(N) # ATTENTION DONNE 1 DANS LE CAS OU IL N'Y A QUE 3 REPLICATS
    Transitivite <- transitivity(Net, type="global")
    #connectivity <- connectedness(N) # ATTENTION EST EGALE A LA CONNECTANCE DANS LE CAS OU IL N'Y A QUE 3 REPLICATS
    Degres = degree (Net) # nombre de liens pour chaque noeud
    Nb_noeud_Connect <- sum(Degres>0) # Nbre de noeuds connecter
    Ratio_noeud_Connect <- Nb_noeud_Connect/Nb_noeuds #Proportion de noeud connecter
    Degres_Moyen <- mean(Degres) #Moyenne des degrès
    #dN.sd <- sd(dN) # Ecart type des degrès
    #dN.max <- max(dN) # 
    Betweenness <- betweenness(Net, directed = F, weights = NA) # ATTENTION DONNE 0 DANS LE CAS OU IL N'Y A QUE 3 REPLICATS - calcul du nombre de chemins qui passent par chaque noeuds
    Betweenness_Moyen <- mean(Betweenness ) # ATTENTION DONNE 0 DANS LE CAS OU IL N'Y A QUE 3 REPLICATS
    Hubs_score <- hub_score(Net)$vector
    hubs_score_Moyen <- mean(Hubs_score)
    modalite <- levels(groupes$modalites)[i]

    
    # creation d'un tableau vide et remplissage
    resume <- data.frame(regnes_lien = c("Archaea_Archaea", "Archaea_Bacteria", "Archaea_Eukaryota", "Archaea_Viruses",
                                         "Bacteria_Bacteria", "Bacteria_Eukaryota", "Bacteria_Viruses", 
                                         "Eukaryota_Eukaryota", "Eukaryota_Viruses", "Viruses_Viruses")) %>% 
      left_join(links %>% 
                  group_by(regnes_lien) %>% 
                  summarise(comptage = n()), 
                by = "regnes_lien") %>% 
      column_to_rownames("regnes_lien") %>% 
      t() %>% as.data.frame()
    
    
    results <- as.data.frame(rbind(results,cbind(modalite, Nb_echant, Nb_noeuds, Nb_noeud_Connect, Ratio_noeud_Connect, Nb_liens, Nblien_Pos, Nblien_Neg,
                                                 Ratio_lien_Pos.Neg, Transitivite, Densite, Degres_Moyen, Betweenness_Moyen, hubs_score_Moyen, resume)))
    
    
    ## Plotting du réseau
    
    # Filtrer et  eliminer les noeud avec moins de 2 liens
    #filter_vertice=filter_vertices=V(Net)[degree(Net)<2]
    filter_vertice=V(Net)[degree(Net)<2]
    Net <- Net%>%delete_vertices(filter_vertice)
    V(Net)$hubs=hubs_score_Moyen
    V(Net)$degres=degree(Net)
    V(Net)$Betwenness=betweenness(Net, directed = F, weights = NA)
    meshubs <- sort(V(Net)$Betwenness, decreasing = T)# Reordonner les valeurs de betweenness par ordre descroissante
    mes20hubs <- meshubs[1:20] #Récuperer les 20 premiers noeud avec les betweenness les plus grands
    
    colors=c("#1A5878", "#ffffb3","#8968CD","#C44237") 

    pdf(file = paste("Graph_reseau_", levels(groupes$modalites)[i], ".pdf", sep=""))
    #set.seed(123)
    print(ggraph(Net, layout = "stress")+
      geom_edge_link0(aes(edge_width=weight, color=direction), show.legend = F)+
      geom_node_point(aes(fill=Superkingdom, size=Betwenness), shape=21)+
      geom_node_text(aes(filter=Betwenness%in%mes20hubs, label=Phylum), family="serif", size=3, repel = T, show.legend = F)+
      scale_fill_manual(values = colors)+
      scale_edge_width(range = c(0.05, 1.5))+
      scale_edge_color_manual(values = c("negative"="red2", "positive"="gray70"))+
      scale_size(range = c(1,6)))
    graphics.off()
    
  }
}


write.csv2(x=results, file=paste(DOSSIER,"/Reseaux_Metagenomique2_Phylum/resultats_reseaux_", nom_modalite, "_", paste(seuil.p, "_", NbEch, "ech_", NbRep, "rep_", AbMinOTUs, "GENRESmin_2021_2022", sep = ""), ".csv", sep = ""), row.names = FALSE)
setwd(DOSSIER)

#############################################################################
#### NB LIENS
library(ggplot2)
library(agricolae)
tableau_donnees <- read.csv2(paste(DOSSIER,"/Reseaux_Metagenomique2_Phylum/resultats_reseaux_modalite_0.01_12ech_1rep_5GENRESmin_2021_2022.csv", sep=""), header = TRUE)
#tableau_donnees <- read.csv2("/home/alphonse/Abt_docs/abtStage_M2/Stage_M2/Reseaux_Metagenomique2_sansT0/resultats_reseaux_modalite_0.06_12ech_1rep_5GENRESmin_2019_2020.csv", header = TRUE)
variable <- "Nb_liens" 
modalite <- "modalite"
#couleur <- "localisation_origine"
nom_axe_x <- NULL
nom_axe_y <- "Nb_liens"
#nom_legende_couleur <- "Loc. Origine"
angle_label_axe_x <- 90

#tableau_donnees$couleur <- tableau_donnees[, colnames(tableau_donnees) == couleur]
tableau_donnees$variable <- as.numeric(as.character(tableau_donnees[, colnames(tableau_donnees) == variable]))
#tableau_donnees$modalite <- factor(tableau_donnees[, colnames(tableau_donnees) == modalite], c("Extensif", "Semi_intensif", "Intensif"))
tableau_donnees$modalite <- tableau_donnees[, colnames(tableau_donnees) == modalite]

ggplot(tableau_donnees, aes(x = modalite, y = variable))+
  geom_boxplot() +
  ylab(nom_axe_y) +
  xlab(nom_axe_x) +
  #ylim(0, 60)+
  #labs (color = nom_legende_couleur)+
  theme (axis.title = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(angle = angle_label_axe_x), 
         axis.text = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  #geom_text(data = kruskal$groups, aes(label = groups, y = max(tableau_donnees$variable)), colour="black", size=8, angle = 90, hjust = "right")+
  geom_jitter(shape=16, position=position_jitter(0.2), col = "purple")+
  stat_summary(fun=mean, geom="point", shape=18, size=4, col = "red")

ggsave(paste("Reseaux_Metagenomique2_Phylum/graphe_metric_network_fonctiomic2_Phylum",seuil.p,".png", sep = ""), width = 8, height = 8)

