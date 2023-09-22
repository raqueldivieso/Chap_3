
# install.packages("phylosignal")
# library(phylosignal)
library(dplyr)
library(ape)
library(geiger)
library(picante)
library(parallel)
library(doParallel)


### Code for analysis in "Phylogenetic patterns of rarity in terrestrial mammals"
### In this code the phylogenetical signal of the seven forms of rarity will be tested. Our hypothesis is that considering the phylogenetical signal of the rarity axis, even combined we will see this pattern.  

source("code.R")


#Range size - PHYLACINE database (https://doi.org/10.5281/zenodo.3690867) 
range<-read.csv("data/Spatial_metadata_PHYLACINE.csv")
range<-select(range,Binomial.1.2, Number.Cells.Present.Natural.Range)
colnames(range)<-c("binomial", "range")
range$binomial<-gsub("_", " ", range$binomial)
table(duplicated(range$binomial))

#Density estimates - Santini et al 2022 (https://doi.org/10.1111/geb.13476) ##PredMd
dens<-read.csv("data/Santini_etal_density_estimates_2022.csv")
table(duplicated(dens$binomial))

# Body mass and habitat breadth - COMBINE database (https://doi.org/10.6084/m9.figshare.13028255.v4)
traits<-read.csv("data/trait_data_imputed_COMBINE.csv")
traits<-select(traits, iucn2020_binomial, order, adult_mass_g, habitat_breadth_n)
traits<-rename(traits, binomial=iucn2020_binomial)

table(duplicated(traits$binomial)) #cheking for duplicated names
dfd<- data.frame(traits %>% group_by(binomial) %>% filter(n() > 1))
#All duplicated species rows have the same mass and habitat breadth, let's keep just one of each.
traits <- traits %>% distinct(binomial, .keep_all = TRUE)

data1<-merge(traits, dens, by="binomial", all.x = T)
data2<-merge(data1, range, by="binomial", all.x = T)
data3<-na.omit(data2)
dim(data3) #We have 4455 to conduct the analisis.


###Calculating the percentil for crate tree classes of size

data3$log_adult_mass_g<-log(data3$adult_mass_g)
percentiles <- quantile(data3$log_adult_mass_g, probs = c(0.25, 0.50, 0.75))

# 25%      50%      75% 
# 2.993229 4.156599 6.159946 

data_ext_small <- subset(data3, log_adult_mass_g <= 2.993229)
data_small <- subset(data3, log_adult_mass_g > 2.993229 & log_adult_mass_g <=4.156599)
data_medium <- subset(data3, log_adult_mass_g > 4.156599 & log_adult_mass_g <=6.159946)
data_large <- subset(data3, log_adult_mass_g > 6.159946)

dim(data_ext_small) #1113     spp
dim(data_small) #1096     spp
dim(data_medium) #1132 spp
dim(data_large) #1114 spp


data_ext_small$dp_l<-log(data_ext_small$PredMd) # predicted density log
data_ext_small$dp<- 1-(data_ext_small$dp_l/max(data_ext_small$dp_l, na.rm = T)) # predicted
data_small$dp_l<-log(data_small$PredMd) # predicted density log
data_small$dp<- 1-(data_small$dp_l/max(data_small$dp_l, na.rm = T)) # predicted density index***
data_medium$dp_l<-log(data_medium$PredMd) # predicted density log
data_medium$dp<- 1-(data_medium$dp_l/max(data_medium$dp_l, na.rm = T))
##
data_large$dp_l<-log(data_large$PredMd) # predicted density log
data_large$dp<- 1-(data_large$dp_l/max(data_large$dp_l, na.rm = T))

#### The medians of each group for density index
median(data_ext_small$dp, na.rm = T) #extra_small = 0.6162482
median(data_small$dp, na.rm = T) #small =  0.3851455
median(data_medium$dp, na.rm = T) #medium =  0.4787729
median(data_large$dp, na.rm = T) #large =  0.7355161

data_ext_small$dp_rare<-ifelse(data_ext_small$dp>=0.6162482, TRUE, FALSE)
data_small$dp_rare<-ifelse(data_small$dp>=0.3851455, TRUE, FALSE)
data_medium$dp_rare<-ifelse(data_medium$dp>=0.4787729, TRUE, FALSE)
data_large$dp_rare<-ifelse(data_large$dp>=0.7355161, TRUE, FALSE)

data_dp<-rbind(data_ext_small, data_small, data_medium, data_large)
head(data_dp)
data_dp<-select(data_dp,binomial, dp, dp_rare )
table(data_dp$dp_rare)

data3$rs_l<-log(data3$range) # range size log
data3$rs<- 1-(data3$rs_l/max(data3$rs_l, na.rm = T)) # range size index***
# hist(data3$rs)

median(data3$rs, na.rm = T)  #0.5486034
data3$rs_rare<- ifelse(data3$rs>=0.5486034, TRUE, FALSE)

#All species with habitat breadth = 1 are considering rare:
data3$hb_rare<- ifelse(data3$habitat_breadth_n==1, TRUE, FALSE)
# table(data3$hb_rare)

dat<-merge(data3, data_dp, by="binomial", all = F, no.dups = TRUE)

dat$rank<- ifelse(dat$rs_rare==TRUE & dat$hb_rare==TRUE & dat$dp_rare==TRUE, "H", ifelse(dat$rs_rare==TRUE & dat$hb_rare==FALSE & dat$dp_rare==TRUE, "G",ifelse(dat$rs_rare==TRUE & dat$hb_rare==TRUE & dat$dp_rare==FALSE, "F", ifelse(dat$rs_rare==TRUE & dat$hb_rare==FALSE & dat$dp_rare==FALSE, "E", ifelse(dat$rs_rare==FALSE & dat$hb_rare==TRUE & dat$dp_rare==TRUE, "D",ifelse(dat$rs_rare==FALSE & dat$hb_rare==FALSE & dat$dp_rare==TRUE, "C", ifelse(dat$rs_rare==FALSE & dat$hb_rare==TRUE & dat$dp_rare==FALSE, "B", "A")))))))

dat$rar_levels<-dat$rank
dat$rar_levels<-sub('A','0', dat$rar_levels)
dat$rar_levels<-sub('B','1', dat$rar_levels)
dat$rar_levels<-sub('C','1', dat$rar_levels)
dat$rar_levels<-sub('D','2', dat$rar_levels)
dat$rar_levels<-sub('E','1', dat$rar_levels)
dat$rar_levels<-sub('F','2', dat$rar_levels)
dat$rar_levels<-sub('G','2', dat$rar_levels)
dat$rar_levels<-sub('H','3', dat$rar_levels)

write.csv(dat, "mammals_rarity.csv") #This data will be used in the next analysis and also for the structural model.

####### Analysis
### Calculating the phylogenetical signal of Rarity (forms and levels)

data<-read.csv("data/mammals_rarity.csv")
data$binomial<-gsub(" ", "_", data$binomial)

#Opening first the consensus tree from Uphan et al., 2019 obteined in: https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_MCC_v2_target.tre
tree<-read.nexus("data/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_namesok.nex") ### The code for the names configuration of this tree are in the following code.
# The final parte of this code is the calculation of the phylogenetical signal for the Uphan et al., 2019 multiple trees.


# The phylogenetical metric used here is one analog of the Shannon entropy (for categorical traits, as the Rarity forms and levels)
# Reference: Rui Borges, João Paulo Machado, Cidália Gomes, Ana Paula Rocha, Agostinho Antunes; Measuring phylogenetic signal between categorical traits and phylogenies, Bioinformatics, https://doi.org/10.1093/bioinformatics/bty800

#First let's guarantee that the names in data and phylogeny are matching:
#Keewping just the levels and forms of rarity:
dat<-select(data, rank, rar_levels)
rownames(dat)<-data$binomial
length(dat$rank)

d<- treedata(tree, dat, sort=TRUE)
rar_forms<-d$data[,1]
rar_levels<-d$data[,2]

delta_forms<- delta(rar_forms,d$phy,0.1,0.0589,1000,10,100)
res1<-paste("Delta - Rarity Forms = ", delta_forms)
write.table(res1, "Results/Result_Rar_forms_PhyloSignal.txt") #2.652576

delta_levels<- delta(rar_levels,d$phy,0.1,0.0589,10000,10,100)
res2<-paste("Delta - Rarity Levels = ", delta_levels)
write.table(res2, "Results/Result_Rar_levels_PhyloSignal") #1.32


### 100 ramdom datasets to calculate the ramdom deltas

# For the forms
random_delta_f <- rep(NA,10)
names(rar_forms)<-NULL
for(i in 1:10){
  rar_forms_ram<-sample(rar_forms)
  random_delta_f[i]<- delta(rar_forms_ram, d$phy,0.1,0.0589,1000,10,100)
  write.csv(random_delta_f, "random_delta_forms.csv")}

# For the levels
random_delta_l <- rep(NA,10)
names(rar_levels)<-NULL
for(i in 1:10){
  rar_levels_ram<-sample(rar_levels)
  random_delta_l[i]<- delta(rar_levels_ram,d$phy,0.1,0.0589,1000,10,100)
  write.csv(random_delta_l, "random_delta_levels.csv")}


### Now let's calculate the phylogenetic signal for each order

### Function for the forms
org_ord_f<-function(ord){
dat_o<- subset(data, order== "Diprotodontia") #selecting the order
rownames(dat_o)<-dat_o$binomial   #set the rownames
td<- treedata(tree, dat_o, sort=TRUE)        #comparing taxa in data and tree
td$data<- td$data  %>%            #Selecting the specific columns
    as.data.frame() %>%
    select(rank, rar_levels)
rar_forms<-td$data[,1]
names(rar_forms)<-rownames(td$data)
delta_forms_dip<- delta(rar_forms,td$phy,0.1,0.0589,1000,10,100) #calc the delta

 
res2<-paste(ord, round(delta_forms, 3), sep=",")
res2}

### Function for the levels
org_ord_l<-function(ord){
dat_o<- subset(data, order== ord) #selecting the order
rownames(dat_o)<-dat_o$binomial   #set the rownames
td<- treedata(tree, dat_o, sort=TRUE)        #comparing taxa in data and tree
td$data<- td$data  %>%            #Selecting the specific columns
    as.data.frame() %>%
    select(rank, rar_levels)
rar_levels<-td$data[,2]
names(rar_levels)<-rownames(td$data)
delta_levels<- delta(rar_levels,td$phy,0.1,0.0589,10000,10,100) #calc the delta
res2<-paste(ord, round(delta_levels, 3), sep=",")
res2}

list_o<-list(as.character(subset(data.frame(table(data$order)), Freq>100)$Var1))[[1]]

#Applying the functions 
res1<-lapply(list_o, org_ord_f)
x<-as.data.frame(unlist(res1))
colnames(x)<-c("order, PhyloSig_forms")
write.csv(x, "Results/PhyloSignal_all_orders_Rar_forms.csv")

delta_f_dip<-org_ord_f("Diprotodontia")
delta_f_eul<-org_ord_f("Eulipotyphla")
delta_f_pri<-org_ord_f("Primates")
delta_f_rod<-org_ord_f("Rodentia")

res2<-lapply(list_o, org_ord_l)
y<-as.data.frame(unlist(res2))
colnames(y)<-c("order, PhyloSig_forms")
write.csv(y, "Results/PhyloSignal_all_orders_Rar_levels.csv")


### Now let's calculate the 100 ramdom datasets and its respective delta for each specie

#Carnivora
dat_o <- subset(data, order== "Carnivora")
rownames(dat_o)<-dat_o$binomial
dat_o <- select(dat_o, rank, rar_levels)
dt_o<-treedata(tree, dat_o)
rar_forms<-dt_o$data[,1]
rar_levels<-dt_o$data[,2]

# For the forms
random_delta_f <- rep(NA,10)
names(rar_forms)<-NULL
for(i in 1:16){
  rar_forms_ram<-sample(rar_forms)
  random_delta_f[i]<- delta(rar_forms_ram, dt_o$phy,0.1,0.0589,1000,10,100)
  write.csv(random_delta_f, "Results/random_delta_forms_car.csv")}

# For the levels
random_delta_l <- rep(NA,10)
names(rar_levels)<-NULL
for(i in 1){
  rar_levels_ram<-sample(rar_levels)
  random_delta_l[i]<- delta(rar_levels_ram, dt_o$phy,0.1,0.0589,100,10,100)
  write.csv(random_delta_f, "Results/random_delta_levels_car.csv")}



setwd("C:/Users/raque/OneDrive - ufpr.br/Área de Trabalho/Ramd_deltas")

p_value <- sum(random_delta>deltaA)/length(random_delta)
boxplot(random_delta)
abline(h=deltaA,col="red")





