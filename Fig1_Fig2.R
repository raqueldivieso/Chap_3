library(dplyr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(showtext)


deltaA<-read.csv("forms_levels_delta.csv")
deltaA
deltaA$order[1]<-"All"
deltaA$order[5]<-"Diprotodontia"

ord<- list(deltaA$order)[[1]]

# Table 1

Phylo_signal<-data.frame(matrix(nrow=8, ncol=7))
colnames(Phylo_signal)<-c("Order", "Delta_real_Forms", "Mean_ramd_delta_Forms", "p-value_Forms", "Delta_real_Levels", "Mean_ramd_delta_Levels", "p-value_Levels")

for (i in 1:8){
  xx<-paste("Ramd_deltas/100_deltas/", ord[i], "_random_delta_forms_100.csv", sep="")
  random_delta<-read.csv(xx)$x
  Phylo_signal$Order[i]<-ord[i]
  Phylo_signal$`p-value_Forms`[i] <- sum(random_delta>deltaA$Forms[i])/length(random_delta)
  Phylo_signal$Mean_ramd_delta_Forms[i]<- mean(random_delta)
  Phylo_signal$Delta_real_Forms[i]<-deltaA$Forms[i]}

for (i in 1:8){
  xx<-paste("C:/Users/raque/OneDrive - ufpr.br/Área de Trabalho/Ramd_deltas/100_deltas/", ord[i], "_random_delta_levels_100.csv", sep="")
  random_delta<-read.csv(xx)$x
  Phylo_signal$`p-value_Levels`[i] <- sum(random_delta>deltaA$Levels[i])/length(random_delta)
  Phylo_signal$Mean_ramd_delta_Levels[i]<- mean(random_delta)
  Phylo_signal$Delta_real_Levels[i]<-deltaA$Levels[i]}

write.csv(Phylo_signal, "Tab.1_Phylo_signal_Levels_CERTA.csv")

###### Fig 1
showtext_opts(dpi = 1500)
showtext_auto(enable = TRUE)

datam <- read.csv("C:/Users/raque/OneDrive - ufpr.br/Área de Trabalho/Ramd_deltas/100_deltas/All_orders_random_delta_forms_100.csv", header=TRUE, sep=",")
datam$X<-NULL

# Get 8 colors from the viridis palette
colors <- viridis_pal()(8)

# Fig. 1

# plot
p1 <- datam %>%
  ggplot( aes(x=All.species)) +
  geom_histogram(fill="#1b0022", color="#e9ecef", width = 0.5) +
  ggtitle("All species") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("Frequency")+
  xlab("Random λDR")+
  geom_vline(xintercept=0.5, col="red", linewidth=1.2, linetype = "dashed")
p1
# 
# png("Fig1_Forms_delta_teste.tiff",  width=18, height=26, units="in", res=500)
# p1
# dev.off()


p2 <- datam %>%
  ggplot( aes(x=Carnivora)) +
  geom_histogram(fill="#46337EFF", color="#e9ecef", width = 0.5) +
  ggtitle("Carnivora") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.748, col="red", linewidth=1.2, linetype = "dashed")
p2

p3 <- datam %>%
  ggplot( aes(x=Cetartiodactyla)) +
  geom_histogram(fill=colors[3], color="#e9ecef", width = 0.5) +
  ggtitle("Cetartiodactyla") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.292, col="red", linewidth=1.2, linetype = "dashed")
p3

p4 <- datam %>%
  ggplot( aes(x=Chiroptera)) +
  geom_histogram(fill=colors[4], color="#e9ecef", width = 0.5) +
  ggtitle("Chiroptera") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.3, col="red", linewidth=1.2, linetype = "dashed")
p4

p5 <- datam %>%
  ggplot( aes(x=Diprotodontia)) +
  geom_histogram(fill=colors[5], color="#e9ecef", width = 0.5) +
  ggtitle("Diprotodontia") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.51, col="red", linewidth=1.2, linetype = "dashed")
p5

p6 <- datam %>%
  ggplot( aes(x=Eulipotyphla)) +
  geom_histogram(fill=colors[6], color="#e9ecef", width = 0.5) +
  ggtitle("Eulipotyphla") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.6, col="red", linewidth=1.2, linetype = "dashed")
p6


p7 <- datam %>%
  ggplot( aes(x=Primates)) +
  geom_histogram(fill=colors[7], color="#e9ecef", width = 0.5) +
  ggtitle("Primates") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=2, col="red", linewidth=1.2, linetype = "dashed")
p7


p8 <- datam %>%
  ggplot( aes(x=Rodentia)) +
  geom_histogram(fill=colors[8], color="#e9ecef", width = 0.5) +
  ggtitle("Rodentia") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1, col="red", linewidth=1.2, linetype = "dashed")
p8


tiff("Fig1_Forms_delta_lambda_copiar a legenda.tiff",  width=18, height=26, units="in", res=500)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
dev.off()


###################
datal <- read.csv("C:/Users/raque/OneDrive - ufpr.br/Área de Trabalho/Ramd_deltas/100_deltas/All_orders_random_delta_levels_100.csv")
datal$X<-NULL

# plot
p1 <- datal %>%
  ggplot( aes(x=All.species)) +
  geom_histogram(fill="#1b0022", color="#e9ecef", width = 0.5) +
  ggtitle("All species") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.32, col="red", linewidth=1.2, linetype = "dashed")
p1


p2 <- datal %>%
  ggplot( aes(x=Carnivora)) +
  geom_histogram(fill="#46337EFF", color="#e9ecef", width = 0.5) +
  ggtitle("Carnivora") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.88, col="red", linewidth=1.2, linetype = "dashed")
p2

p3 <- datal %>%
  ggplot( aes(x=Cetartiodactyla)) +
  geom_histogram(fill=colors[3], color="#e9ecef", width = 0.5) +
  ggtitle("Cetartiodactyla") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=0.731, col="red", linewidth=1.2, linetype = "dashed")
p3

p4 <- datal %>%
  ggplot( aes(x=Chiroptera)) +
  geom_histogram(fill=colors[4], color="#e9ecef", width = 0.5) +
  ggtitle("Chiroptera") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.404, col="red", linewidth=1.2, linetype = "dashed")
p4

p5 <- datal %>%
  ggplot( aes(x=Diprotodontia)) +
  geom_histogram(fill=colors[5], color="#e9ecef", width = 0.5) +
  ggtitle("Diprotodontia") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=0.741, col="red", linewidth=1.2, linetype = "dashed")
p5

p6 <- datal %>%
  ggplot( aes(x=Eulipotyphla)) +
  geom_histogram(fill=colors[6], color="#e9ecef",width = 0.5) +
  ggtitle("Eulipotyphla") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=0.684, col="red", linewidth=1.2, linetype = "dashed")
p6


p7 <- datal %>%
  ggplot( aes(x=Primates)) +
  geom_histogram(fill=colors[7], color="#e9ecef", width = 0.5) +
  ggtitle("Primates") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=1.336, col="red", linewidth=1.2, linetype = "dashed")
p7


p8 <- datal %>%
  ggplot( aes(x=Rodentia)) +
  geom_histogram(fill=colors[8], color="#e9ecef", width = 0.5) +
  ggtitle("Rodentia") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5 , size = 12))+
  ylab("")+
  xlab("")+
  geom_vline(xintercept=0.5, col="red", linewidth=1.2, linetype = "dashed")
p8


tiff("Fig2_Forms_delta_lambda_arrumarridentia.tiff",  width=18, height=26, units="in", res=500)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
dev.off()





