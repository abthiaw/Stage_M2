# Stage_M2 : ETUDE METAGENOMIQUE DE LA DIVERSITE TAXONOMIQUE DES COMMUNAUTES MICROBIENNES DU SOLS
==================================================================================================
By Alphonse Birane Thiaw




------------
 1) RESUME |
------------

Ce projet consiste à caractériser l’influence de deux modes d’usage contrastés sur les communautés microbiennes de sols avec et, à déterminer le
rapport entre ces changements et le recyclage de la matière organique a partir de l'analyse de données de sequencge metagenomique sans a priori.
Il se deroule en 2 etapes: Analyse de la diversite taxonomiques microbiennes, Reseau de co-occurences pour évaluer les différences d’interactions entre
micro-organismes apres apport de matiere orgsniques dans les sols.

# Pour executer les scripts en local 
- git clone git@github.com:abthiaw/Stage_M2.git
- ou telecharger le repertoire

- les Donnees analyses sont des donnees rarefies et agglomerer au niveau du genre au prealable puis stocker au format RData.

----------------------------
2) DIVERSITES TAXONOMIQUES |
----------------------------

 - Execution du script <2022_05-Fonctiomic_script_taxonomique_analysis.R> donne les resultats suivants
     ==> Courbe d'accumulation evaluant l'effort d'echantillonnage
     ==> Analyse de l'abondances relative des communautes microbiennes
     ==> Etudes de la diversites taxonomiques microbiennes
         - Bases sur les nombres de Hill
     ==> Etudes de la structuration des communautes miicrobiennes
         - Analyse d'ordination PCoA et NMDS
     ==> Etude de la difference d'abondance microbiennes selon les conditions d'etudes
  
 <2022_05-Fonctiomic_script_taxonomique_analysis.rmd> est la version Rmarkdown qui donne en output un rapport de toutes les analyses au format pdf
     
         

------------------------------
2) RESEAUX DE CO-OCCURENCES  |
------------------------------
- Execution du script <2022_05-Fonctiomic_script_taxonomique_analysis.R> permet d'elaborer des reseaux de co-occurences entre communautes microbiennes
