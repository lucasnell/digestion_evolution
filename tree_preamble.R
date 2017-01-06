# Required packages
# Commented lines are packages I may need later, but not yet.

library(dplyr)
library(readr)
# library(tidyr)
# library(ggplot2)
library(ape)
library(rotl)
library(phytools)
# library(phylolm)
# library(vegan)
# library(lme4)
# library(pez)
# library(lmerTest)



# Data frame of species needed:
    
sp_df <- 
    'species,type
Mus musculus,rodent
Rattus norvegicus,rodent
Thaptomys nigrita,rodent
Akodon montensis,rodent
Delomys sublineatus,rodent
Olygoryzomys nigripes,rodent
Sooretamys angouya,rodent
Euryoryzomys russatus,rodent
Peromyscus leucopus,rodent
Microtus pennsylvanicus,rodent
Artibeus lituratus,bat
Carollia perspicillata,bat
Desmodus rotundus,bat
Molossus rufus,bat
Molossus molossus,bat
Eumops glaucinus,bat
Tadarida brasiliensis,bat
Myotis lucifugus,bat
Eptesicus fuscus,bat
' %>% 
    read_csv(.)


