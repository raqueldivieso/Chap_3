library(data.table)
library(ape)
library(ggplot2)
library(dplyr)
library(phytools)
library(parallel)
library(doParallel)

# Code to calculate λDR values for the 1000 trees possibilities (Only DNA - bird death node dated) with 4098 species: 

#The tree reference:  tree-pruner-9c8a7b1a-48fa-4e6c-8acf-b6f77c675a95
trees<-read.nexus("data/PhyMammals_DNA_only_node-dated_1000.nex")

### Once we have 1000 different trees to calculate the diversification rates (λDR) between all species, the computacional demand is too high. Lets perform the analysis for 500 trees at time and combine the results:

### λDR statistic calculation available in: https://github.com/singhal/brazil_IBD/blob/main/analyses/speciesSpecificDivRate.R
#Calculates species-specific diversification rate, as per Jetz et al. 2012
#returns a vector of diversification rates, each element named by tip label.
#author of the following function: Pascal Title
#date: 31 Jan 2013

jetzDivRates <- function(tree) {
  
  spRate <- function(sp, tree) {
    #get branch lengths from root to tip
    edges <- vector()
    daughterNode <- match(sp, tree$tip.label)
    while (daughterNode != (length(tree$tip.label) + 1)) {
      parentNode <- tree$edge[which(tree$edge[,2] == daughterNode), 1]
      edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentNode & tree$edge[,2] == daughterNode)])
      daughterNode <- parentNode
    }
    
    res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
    res <- res ^ (-1)
    
    return(res)
  }
  
  rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
  names(rates) <- tree$tip.label
  
  return(rates)
  
}

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
clusterExport(cl, varlist = c("jetzDivRates","trees"))
clusterEvalQ(cl, library("phytools"))


Sys.time()
result2 <- parLapply(cl, fun = jetzDivRates, trees)
Sys.time()
length(result)
tab_res<-data.frame(result[[1]])
colnames(tab_res)<-"tree_1"
rownames(tab_res)

for (i in 2:1000){
  
  dr_tr<-data.frame(result[[i]])
  colnames(dr_tr)<-(paste("tree", i, sep="_"))
  tab_res<-merge(tab_res,dr_tr, by=0)
  rownames(tab_res)<-tab_res$Row.names
  tab_res<-tab_res[,-1]
}
# write.csv(tab_res, "LambdaDR_tr_1-500.csv")
stopCluster(cl)

# Now, for account for the phylogenetic uncertancy, we need to calculating the variance for 1000 lDR values and comparing the correlation between the consensus phylogeny rates and the average between the 1000 values.
spec_rates<-read.csv("data/LambdaDR_tr_1-1000.csv")
rownames(spec_rates)<-spec_rates$X
head(spec_rates[,1])
spec_rates <- spec_rates[,-1]

#transpose the table
df_t <- data.frame(transpose(spec_rates))
rownames(df_t) <- colnames(spec_rates)
colnames(df_t) <- rownames(spec_rates)

tab_var_rates<- data.frame(matrix(ncol = 3, nrow = 4099))
colnames(tab_var_rates) <- c("specie", "variance", "mean")
for(i in 1:4099){
  tab_var_rates$variance[i]<-var(df_t[,i], na.rm = T)
  tab_var_rates$mean[i]<-mean(df_t[,i], na.rm = T)
  tab_var_rates$specie[i]<-colnames(df_t)[i]}

# Now, we need to open the consensus tree (DNA-only) and match the names that are a quite different format.
# Upham, N. S., J. A. Esselstyn, and W. Jetz. 2019 - Fossilized birth-death, 4098 species + 76 fossil tips, backbone topology as in Zhou et al. (2013)

tr<-read.nexus("data/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_FBDasZhouEtAl_MCC_v2_target.nex")

phy_names<-tr$tip.label
#The names are like "Sorex_minutus_SORICIDAE_EULIPOTYPHLA" - except for the first five rows

phy_names <- as.data.frame(do.call('rbind', strsplit(as.character(phy_names),'_', fixed=TRUE)))
class(phy_names)

phy_names[1, 1]<- "Anolis"
phy_names[1, 2]<- "carolinensis"
phy_names[1, 3:4]<- "X"
head(phy_names)
phy_names_ok<-paste(phy_names$V1, phy_names$V2, sep="_")
length(phy_names_ok)
length(tr$tip.label)

table(phy_names$V1=="X") #76 tips are for node dated, this lets us 4100 species, and removing the external group specie (Anolis carolinensis) we have exacly 4099 species.

# Just checking:
phy_names_ok[c(30, 990, 1477, 2222, 3799, 4176)]
tr$tip.label[c(30, 990, 1477, 2222, 3799, 4176)]
tr$tip.label<-phy_names_ok

# Now, lets calculate the evolutionary rates for the consensus tree:

conc_tr_rates<-jetzDivRates(tr)
conc_tr_ratesdf<-data.frame(conc_tr_rates)
conc_tr_ratesdf$specie<-rownames(conc_tr_ratesdf)
dat<-merge(conc_tr_ratesdf,tab_var_rates, by="specie", all=T )
head(dat)
plot(dat$conc_tr_rates, dat$mean) #we have a R² of 88% in the relashionship/
mod<-lm(dat$mean~dat$conc_tr_rates)
abline(mod)
summary(mod)

min(as.matrix(df_t[,1:4099]))
max(as.matrix(df_t[,1:4099]))

spec_rates$specie<-rownames(spec_rates)

dat2<-merge(conc_tr_ratesdf, spec_rates, by="specie")
tab_res<-data.frame(matrix(ncol = 5, nrow = 1000))
colnames(tab_res) <- c("tree", "p_value", "intercept", "slope", "R2")

for (i in 1:1000){
  reg<-lm(dat2$conc_tr_rates~dat2[,i])
  res<-summary(reg)
  tab_res$tree[i]<- colnames(dat2)[i]
  tab_res$p_value[i] <- res$coefficients[2,4] 
  tab_res$intercept <- res$coefficients[1,1] 
  tab_res$slope[i] <- res$coefficients[2,1]
  tab_res$R2[i] <- res$r.squared} 

### Fig for Material suppementar
png("FigS2.png",  width=6, height=6, units="in", res=300)
plot(1, type='n',xlim=c(min(dat2$conc_tr_rates), max(dat2$conc_tr_rates-0.4)), ylim=c(min(dat2$conc_tr_rates), max(dat2$conc_tr_rates)), main=NULL, cex.axis=1,cex.lab=1.2, xlab="Speciation rates (λDR) for each of the 1000 possible phylogenies",ylab="Speciation rates (λDR) for consensus phylogeny")


for (i in 2:1000
     ) { 
  abline(a=0.0, b = tab_res$slope[i], col = alpha("#73aae0", 0.1))}
points(dat$mean, dat$conc_tr_rates, pch=20, col="#2b77c2", cex = 0.5)
abline(a=0, b = mean(tab_res$slope), col="#2b77c2", lwd = 1.5)
dev.off()

min(dat2$conc_tr_rates)
max(as.matrix(dat2[3:502]))
min((as.matrix(dat2[3:502])))

write.nexus(tr, file = "data/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_namesok.nex", translate = TRUE) ###This tree consensus will be used in the other analysis.
write.csv(conc_tr_ratesdf, "data/conc_tr_rates_speciation.csv") #the speciation rates also will be used in the subsequencial analysis. 





