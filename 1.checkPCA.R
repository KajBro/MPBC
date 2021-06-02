## ANCESTRY PLOT ##

hap3 <- read.table("hapMapPopulationsWithNamesChanged.txt", header = T)
# Here you put you output from --pca in Plink
pcsTemp <- read.table("pca.eigenvec", header = T)
# Here you put you output after merging your data with HapMap3 populations (see ancestry script)
fam <- read.table("hapmap3_bin_snplis.fam", header = F)
famTemp <- fam[,c("V1","V2","V6")]
colnames(famTemp) <- c("FID","IID","CASE")
famTemp$CASE[famTemp$CASE == -9] <- NA
colnames(pcsTemp) <- c ("FID", "IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
famTemp$index <- paste(famTemp$FID,famTemp$IID,sep = "_")
pcsTemp$index <- paste(pcsTemp$FID,pcsTemp$IID,sep = "_")
pcs <- merge(pcsTemp, famTemp, by = "index", all.x = T)
hap3$index <- paste(hap3$FID,hap3$IID,sep = "_")
pcs$CASE[pcs$CASE == 1] <- "CONTROL"
pcs$CASE[pcs$CASE == 2] <- "CASE"
pcs$CASE[is.na(pcs$CASE)] <- "HAPMAP"
library(stringr)
# pcs$tempIID <-  gsub("_.+$", "", pcs$IID)
# pcs$index <- paste(pcs$FID,pcs$tempIID,sep = "_")
data <- merge(pcs,hap3, by = "index", all.x = T)
data$POPULATION <- as.character(data$population)
data$POPULATION[data$CASE == "CASE"] <- "CASE"
data$POPULATION[data$CASE == "CONTROL"] <- "CONTROL"
data$POPULATION2 <- ifelse(data$CASE != "CASE" & data$CASE != "CONTROL", as.character(data$population), "PD")
#pop2 in others is FID
dat <- subset(data, !is.na(data$POPULATION))
write.table(dat, "AnnotatedPopulationCLusters", quote = F,sep = "\t",row.names = F)
library(ggplot2)
library(grid)
data <- read.table("AnnotatedPopulationCLusters", header = T)
plot1 <- ggplot(data, aes(PC1, PC2)) + labs(title="COHORT with HapMap3 Populations") + geom_point(aes(colour = POPULATION2),alpha = 1/4, jitter = "position") + theme_bw()
dat <- subset(data, POPULATION == "CEU" | POPULATION == "TSI" | POPULATION == "ASW" | POPULATION == "LWK" | POPULATION == "MKK" | POPULATION == "YRI" | POPULATION == "CASE" | POPULATION == "CONTROL")
#cut off below based on visual inspection of clusters in plot1 and 6 SD of European mean
euros <- subset(data, POPULATION == "CEU" | POPULATION == "TSI")
eurosLowC1 <- mean(euros$PC1) - (6*sd(euros$PC1)) #Have changed from 9 SD to 6
eurosHighC1 <- mean(euros$PC1) + (6*sd(euros$PC1))
eurosLowC2 <- mean(euros$PC2) - (6*sd(euros$PC2))
eurosHighC2 <- mean(euros$PC2) + (6*sd(euros$PC2))
data$CONTINENT <- "ASIA"
data$CONTINENT[data$PC1 >= eurosLowC1 & data$PC1 <= eurosHighC1 & data$PC2 >= eurosLowC2 & data$PC2 <= eurosHighC2] <- "EUROPE"
data$CONTINENT[data$CONTINENT != "EUROPE" & data$PC2 < mean(euros$PC2)] <- "AFRICA"
onlyGwas <- subset(data, POPULATION2 == "PD")
africa <- subset(data, CONTINENT == "AFRICA")
asia <- subset(data, CONTINENT == "ASIA")
europe <- subset(data, CONTINENT == "EUROPE")
summary(africa)
summary(asia)
summary(europe)
plot2 <- ggplot(onlyGwas, aes(PC1, PC2)) + labs(title="Case and Control only") + geom_point(aes(colour = POPULATION),alpha = 1/4, jitter = "position") + theme_bw()
plot3 <- ggplot(asia, aes(PC1, PC2)) + labs(title="Asian ancestry only") + geom_point(aes(colour = POPULATION),alpha = 1/4, jitter = "position") + theme_bw()
plot4 <- ggplot(africa, aes(PC1, PC2)) + labs(title="African Ancestry only") + geom_point(aes(colour = POPULATION),alpha = 1/4, jitter = "position") + theme_bw()
plot5 <- ggplot(europe, aes(PC1, PC2)) + labs(title="European ancestry only") + geom_point(aes(colour = POPULATION),alpha = 1/4, jitter = "position") + theme_bw()
#write.table(paste(dat4$FID.x,dat4$IID.x, sep = " "), "pdbpHbs_EurAncIdDs.txt",quote = F,sep = "\t", row.names = F,col.names = F)
library(grid)
library(gridExtra)
pdf(file = "popstrat.pdf",width=8,height=16)
grid.arrange(plot1,plot2,plot3,plot4,plot5,nrow=5)
dev.off()

africaNeuroX <- subset(onlyGwas, CONTINENT == "AFRICA")
asiaNeuroX <- subset(onlyGwas, CONTINENT == "ASIA")
europeNeuroX <- subset(onlyGwas, CONTINENT == "EUROPE")
write.table(paste(africaNeuroX$FID.x,africaNeuroX$IID.x,sep = " "),"africanAncestrySampleIds.txt",quote = F,row.names = F,col.names = F)
write.table(paste(asiaNeuroX$FID.x,asiaNeuroX$IID.x,sep = " "),"asianAncestrySampleIds.txt",quote = F,row.names = F,col.names = F)
write.table(paste(europeNeuroX$FID.x,europeNeuroX$IID.x,sep = " "),"europeanAncestrySampleIds.txt",quote = F,row.names = F,col.names = F)
