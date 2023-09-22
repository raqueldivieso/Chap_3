### Phylopath analysis
# https://peerj.com/articles/4718/
# https://cran.r-project.org/web/packages/phylopath/vignettes/intro_to_phylopath.html

#### Phylogenetic path analysis
### In this code a structural model will be proposed to investigate the relashionship between rarity traits and speciation rates

library(phylopath)
library(ape)
library(dplyr)
library(ggm)
library(geiger)

setwd("C:/Users/raque/OneDrive - ufpr.br/Área de Trabalho/3_cap_final")

### Opening the tree and data (traits and speciation rates)
tree <-read.nexus("data/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_namesok.nex") 
rates<-read.csv("data/conc_tr_rates_speciation.csv")
rates<-dplyr::select(rates, specie, conc_tr_rates)
rates<-rename(rates, binomial=specie)
data<-read.csv("data/mammals_rarity.csv")

lat<-read.csv("data/hig_low_midpoint_by_spp_all.csv")
lat<-dplyr::select(lat, binomial, mp_lat)
dat<-dplyr::select(data, binomial, order, adult_mass_g, habitat_breadth_n, PredMd, range)
dat$binomial<-gsub(" ", "_", dat$binomial)
dat<-merge(dat, lat, by="binomial")
dat<-merge(dat, rates, by="binomial")

colnames(dat)<- c("binomial", "order", "body", "hab", "dens", "range", "lat", "speciation")

dat$body<-log(dat$body)
dat$dens<-log(dat$dens)
dat$range<-log(dat$range)
dat$lat<-abs(dat$lat)
dat$speciation<-log(dat$speciation)
rownames(dat)<-dat$binomial
dat<-na.omit(dat)
dim(dat)

# body = adult mass (log(g))
# hab = habitat breadht (numb of diff enviromn)
# dens = population density (log(ind/km²))
# range = range size (log(pixels = 95 km²))
# biome = predominant biome of the distribution (factors: desert, boreal, temper, dry tropics, humid tropics)
# speciation= tip rates by lambdaDR statistic (consensus tree)

# We have 3173 species to perform the analysis.

# m <- define_model_set(
#   null = c(),
#   direct = c(Status~Br),
#   indirect = c(L~Br, G~Br, W~Br),
#   both = c(Status~Br, L~Br, G~Br, W~Br),
#   .common = c(Br~B, P~B, L~B+G, W~G, Status~P+L+G+W+B)
# )

dt<-treedata(tree, dat, sort=T)
dat<-data.frame(dt$data)
dat <- dat %>% mutate_at(c("body", "hab", "dens", "range", "lat", "speciation"), as.numeric)
str(dat)

conc <- define_model_set(
  null = c(),
  direct = c(speciation~range, speciation~dens, speciation~hab, speciation~lat),
  interrel = c(range~hab, range~lat, range~dens, dens~body, dens~lat),
  both = c(speciation~range, speciation~dens, speciation~hab, speciation~lat,hab~range, range~lat, range~dens, dens~body),
  .common = c(speciation~range+dens+hab+lat, speciation~range*dens, speciation~range*hab))

result <- phylo_path(conc, data = dat, tree = dt$phy, model = 'lambda')
summary(result)
result$d_sep
result$model_set
best_model <- best(result)
best_model$se

plot(best_model)

coef_plot(best_model, error_bar = "se", order_by = "strength", to = "Status") + ggplot2::coord_flip()

positions <- data.frame(
  name = c('range', 'dens', 'body', 'lat', 'hab', 'speciation'),
  x = c(2, c(1, 1.75, 3.25, 4), 2.5),
  y = c(3, 3, 2, 2, 2, 2)
)






