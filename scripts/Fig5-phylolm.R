######### Pagels lambda calculations for Nif data #########
# Last Edited 05/20/2025 by Morgan Sobol


setwd("/Users/zaj33/Library/Mobile Documents/com~apple~CloudDocs/Work/PostDoc/nitrogenase/nif_phylogenomics/phylolm/")

library(seecolor)
library(ape)
library(dplyr)
library(phytools)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(ggplotify)
library(phylolm)
library(ggpmisc)
library(ggpubr)

# Import tree
tree <- read.newick("../GenomeData/genomeTrees/tree.nwk")
is.rooted(tree)

# Load metadata
tree_metadata <- read.csv("../metadata/Nif_genome_metadata_final.csv", row.names = 1)

# Replace blank phylum values with NA
tree_metadata$phylum[tree_metadata$phylum == ""] <- NA

# Phylum color key
palette_breaks <-c("Euryarchaeota","Crenarchaeota", "Thaumarchaeota",
                   "Acidobacteria", "Actinobacteria","Aquificae", 
                   "Armatimonadetes","Bacteroidetes", "Balneolaeota", 
                   "Caldiserica", "Calditrichaeota", "Candidate phylum",
                   "Chlamydiae","Chlorobi", "Chloroflexi",
                   "Chrysiogenetes","Cyanobacteria", "Deferribacteres",
                   "Deinococcus-Thermus", "Dictyoglomi","Elusimicrobia",
                   "Fibrobacteres","Firmicutes","Fusobacteria",
                   "Gemmatimonadetes","Ignavibacteriae","Kiritimatiellaeota",
                   "Lentisphaerae","Nitrospinae","Nitrospirae",
                   "Planctomycetes","Deltaproteobacteria","Alphaproteobacteria", 
                   "Acidithiobacillia", "Betaproteobacteria", "Gammaproteobacteria",
                   "Epsilonproteobacteria","Zetaproteobacteria", "Spirochaetes",
                   "Synergistetes","Tenericutes","Thermodesulfobacteria",
                   "Thermotogae","Verrucomicrobia")



palette_tax <-c("#864B61","#9A577E","#C78DAA", 
                "#4B81A3","#8AA86A","#7AA6B1", 
                "#72AC70","#D9E0BE","#9D8F7E",
                "#CFD08E","#6845A1","#D3D3D3", 
                "#CBBD9D","#E6D19D","#51906B",
                "#1A615C","#497970","#3D4066", 
                "#598A55","#E1E0D1","#675974", 
                "#B8862F","#94C1C1","#E9C073", 
                "#D7B7AF","#DFE5EE","#3C0D82", 
                "#CAD5F2","#8CB1B3","#9BC4EC", 
                "#7C7799","#6564B3","#786F95",
                "#BC89C1","#CAA8E4","#D2C7D4",
                "#DAD1DC","#8B32B1","#AA856B",
                "#F1F3B9","#9CBBA2","#7496E0",
                "#D7E0E5","#AAA8AF")

### Habitat coverage ~ Gene richness

# Ensure Habitat is still a data frame and contains both habitat_coverage and Gene
# Assuming `Gene` is another column in `tree_metadata`
Habitat <- data.frame(habitat_coverage = tree_metadata[["habitat_coverage"]],
                      gene_richness = tree_metadata[["gene_richness"]],
                      phylum = tree_metadata[["phylum"]],
                      genome_ID = rownames(tree_metadata))

# Remove rows with NA values 
Habitat <- Habitat[ !is.na(Habitat$habitat_coverage) &
                      Habitat$habitat_coverage != 0 &
                      !is.na(Habitat$gene_richness) &
                      Habitat$gene_richness != 0, ]

# Remove tip branches not in trait dataframe
treeH <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% Habitat$genome_ID])

treeH <- ladderize(treeH, right=TRUE)
is.binary(treeH)
is.rooted(treeH)
treeH <-multi2di(treeH)

# make sure metadata is sorted in order of tip labels in the tree
order <- cbind(treeH$tip.label)
order <- as.data.frame(order, seq=seq(1:length(order)))

# Ensure the row order matches the order in the tree
Habitat <- Habitat[order(match(Habitat$genome_ID, treeH$tip.label)), ]
rownames(Habitat) <- Habitat$genome_ID

# model with no transformation
habitat.pagels = phylolm(habitat_coverage~gene_richness, data=Habitat, phy=treeH, model="lambda")
summary(habitat.pagels)

# check assumptions
plot(fitted(habitat.pagels), residuals(habitat.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

# log transformed dependent
log.habitat.pagels = phylolm(log(habitat_coverage)~gene_richness, data=Habitat, phy=treeH, model="lambda")
summary(log.habitat.pagels)

# check assumptions
plot(fitted(log.habitat.pagels), residuals(log.habitat.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted - Log")
abline(h = 0, col = "red")

# sqrt transformed dependent
sqrt.habitat.pagels = phylolm(sqrt(habitat_coverage)~gene_richness, data=Habitat, phy=treeH, model="lambda")
summary(sqrt.habitat.pagels)

# check assumptions
plot(fitted(sqrt.habitat.pagels), residuals(sqrt.habitat.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted - Sqrt")
abline(h = 0, col = "red")

# standard linear model

habitat_lm = lm(formula = sqrt(habitat_coverage)~gene_richness, data=Habitat)
summary(habitat_lm)

# final pub plot
p1 = ggplot(Habitat, aes(x = gene_richness, y = sqrt(habitat_coverage))) + geom_point(aes(colour = phylum)) + theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +scale_x_continuous(
  limits = c(0.0, 0.8), 
  breaks = seq(0.0, 0.8, 0.1), 
  expand = c(0, 0))+  
  xlab("N2-Fixation Gene Richness") + ylab("Habitat Richness") +
  scale_color_manual(values=palette_tax,
                     breaks=palette_breaks,
                     labels=palette_breaks,
                     name="phylum", na.value="black")
p1 = p1 + geom_abline(intercept=coef(sqrt.habitat.pagels)[1], slope=coef(sqrt.habitat.pagels)[2], linetype = "dashed", color = "blue" ) 
p1 = p1 + geom_abline(intercept=coef(habitat_lm)[1], slope=coef(habitat_lm)[2], linetype = "solid", color = "blue" ) 

p1


### Metabolism coverage ~ Gene richness

Metabolism <- data.frame(metabolic_coverage = tree_metadata[["metabolic_coverage"]],
                      gene_richness = tree_metadata[["gene_richness"]],
                      phylum = tree_metadata[["phylum"]],
                      genome_ID = rownames(tree_metadata))

# Remove rows with NA values 
Metabolism <- Metabolism[!is.na(Metabolism$metabolic_coverage) & !is.na(Metabolism$gene_richness), ]

# Remove tip branches not in trait dataframe
treeM <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% Metabolism$genome_ID])

treeM <- ladderize(treeM, right=TRUE)
is.binary(treeM)
is.rooted(treeM)
treeH <-multi2di(treeM)

# make sure metadata is sorted in order of tip labels in the tree
order <- cbind(treeM$tip.label)
order <- as.data.frame(order, seq=seq(1:length(order)))

# Ensure the row order matches the order in the tree
Metabolism <- Metabolism[order(match(Metabolism$genome_ID, treeM$tip.label)), ]
rownames(Metabolism) <- Metabolism$genome_ID

Metabolism=Metabolism[order(match(Metabolism[,2],order[,1])),]
rownames(Metabolism)=Metabolism$genome_ID

# no transformation model
metabo.pagels = phylolm(metabolic_coverage~gene_richness, data=Metabolism, phy=treeM, model="lambda")
summary(metabo.pagels)

# check assumptions
plot(fitted(metabo.pagels), residuals(metabo.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

# log transformed dependent 
log.metabo.pagels = phylolm(log(metabolic_coverage)~gene_richness, data=Metabolism, phy=treeM, model="lambda")
summary(log.metabo.pagels)

# check assumptions
plot(fitted(log.metabo.pagels), residuals(log.metabo.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted - Log")
abline(h = 0, col = "red")


# sqrt transformed dependent
sqrt.metabo.pagels = phylolm(sqrt(metabolic_coverage)~gene_richness, data=Metabolism, phy=treeM, model="lambda")
summary(sqrt.metabo.pagels)

# check assumptions
plot(fitted(sqrt.metabo.pagels), residuals(sqrt.metabo.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted - Sqrt")
abline(h = 0, col = "red")

#standard linear model

metabo_lm = lm(formula = sqrt(metabolic_coverage)~gene_richness, data=Metabolism)
summary(metabo_lm)

# final pub figure
p2 = ggplot(Metabolism, aes(x = gene_richness, y = sqrt(metabolic_coverage))) + geom_point(aes(colour = phylum)) + theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +scale_x_continuous(
  limits = c(0.0, 0.8), 
  breaks = seq(0.0, 0.8, 0.1), 
  expand = c(0, 0))+
  xlab("N2-Fixation Gene Richness") + ylab("Metabolic Richness") +
  scale_color_manual(values=palette_tax,
                     breaks=palette_breaks,
                     labels=palette_breaks,
                     name="phylum", na.value="black")
p2 = p2 + geom_abline(intercept=coef(sqrt.metabo.pagels)[1], slope=coef(sqrt.metabo.pagels)[2], linetype = "dashed", color = "blue" ) 
p2 = p2 + geom_abline(intercept=coef(metabo_lm)[1], slope=coef(metabo_lm)[2], linetype = "solid", color = "blue" ) 

p2


### Genome size ~ Gene richness

Genome <- data.frame(genome_size = tree_metadata[["genome_size"]],
                     gene_richness = tree_metadata[["gene_richness"]],
                     phylum = tree_metadata[["phylum"]],
                     genome_ID = rownames(tree_metadata))

# Remove rows with NA values in genome size
Genome <- Genome[!is.na(Genome$genome_size) & !is.na(Genome$gene_richness), ]

# Remove tip branches not in trait dataframe
treeG <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% Genome$genome_ID])

treeG <- ladderize(treeG, right=TRUE)
is.binary(treeG)
is.rooted(treeG)
treeH <-multi2di(treeG)

# make sure metadata is sorted in order of tip labels in the tree
order <- cbind(treeG$tip.label)
order <- as.data.frame(order, seq=seq(1:length(order)))

# Ensure the row order matches the order in the tree
Genome <- Genome[order(match(Genome$genome_ID, treeG$tip.label)), ]
rownames(Genome) <- Genome$genome_ID

Genome=Genome[order(match(Genome[,2],order[,1])),]
rownames(Genome)=Genome$genome_ID

# non-transformed model
genome.pagels = phylolm(genome_size~gene_richness, data=Genome, phy=treeG, model="lambda")
summary(genome.pagels)

plot(fitted(genome.pagels), residuals(genome.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

# log transformed dependent 
log.genome.pagels = phylolm(log(genome_size)~gene_richness, data=Genome, phy=treeG, model="lambda")
summary(log.genome.pagels)

# check assumptions
plot(fitted(log.genome.pagels), residuals(log.genome.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted - Log")
abline(h = 0, col = "red")

# sqrt transformed dependent 
sqrt.genome.pagels = phylolm(sqrt(genome_size)~gene_richness, data=Genome, phy=treeG, model="lambda")
summary(sqrt.genome.pagels)

plot(fitted(sqrt.genome.pagels), residuals(sqrt.genome.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted - Sqrt")
abline(h = 0, col = "red")

# standard linear model

genome_lm = lm(formula = log(genome_size)~gene_richness, data=Genome)
summary(genome_lm)


# final pub fig
p3 = ggplot(Genome, aes(x = gene_richness, y = log(genome_size))) + geom_point(aes(colour = phylum)) + theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +scale_x_continuous(
  limits = c(0.0, 0.8), 
  breaks = seq(0.0, 0.8, 0.1), 
  expand = c(0, 0))+ 
  xlab("N2-Fixation Gene Richness") + ylab("Genome Size") +
  scale_color_manual(values=palette_tax,
                     breaks=palette_breaks,
                     labels=palette_breaks,
                     name="phylum", na.value="black")
p3 = p3 + geom_abline(intercept=coef(log.genome.pagels)[1], slope=coef(log.genome.pagels)[2], linetype = "dashed", color = "blue" ) 
p3 = p3 + geom_abline(intercept=coef(genome_lm)[1], slope=coef(genome_lm)[2], linetype = "solid", color = "blue" ) 

p3

### Temperature ~ Gene richness
Temp <- data.frame(temp = tree_metadata[["temperature"]],
                   gene_richness = tree_metadata[["gene_richness"]],
                   phylum = tree_metadata[["phylum"]],
                   genome_ID = rownames(tree_metadata))

# Remove rows with NA values in genome size
Temp <- Temp[!is.na(Temp$temp) & !is.na(Temp$gene_richness), ]

# Remove tip branches not in trait dataframe
treeT <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% Temp$genome_ID])

treeT <- ladderize(treeT, right=TRUE)
is.binary(treeT)
is.rooted(treeT)
treeT <-multi2di(treeT)
#tree2<-midpoint.root(tree2)

# make sure metadata is sorted in order of tip labels in the tree
order <- cbind(treeT$tip.label)
order <- as.data.frame(order, seq=seq(1:length(order)))

# Ensure the row order matches the order in the tree
Temp <- Temp[order(match(Temp$genome_ID, treeT$tip.label)), ]
rownames(Temp) <- Temp$genome_ID

Temp=Temp[order(match(Temp[,2],order[,1])),]
rownames(Temp)=Temp$genome_ID

# non-transformed model
temp.pagels = phylolm(temp~gene_richness, data=Temp, phy=treeT, model="lambda")
summary(temp.pagels)

plot(fitted(temp.pagels), residuals(temp.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")


# log-transformed model
log.temp.pagels = phylolm(log(temp)~gene_richness, data=Temp, phy=treeT, model="lambda")
summary(log.temp.pagels)


plot(fitted(log.temp.pagels), residuals(log.temp.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted - Log")
abline(h = 0, col = "red")

# sqrt-transformed model
sqrt.temp.pagels = phylolm(sqrt(temp)~gene_richness, data=Temp, phy=treeT, model="lambda")
summary(sqrt.temp.pagels)


plot(fitted(sqrt.temp.pagels), residuals(sqrt.temp.pagels), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted - Sqrt")
abline(h = 0, col = "red")

# standard linear model 

temp_lm = lm(formula = log(temp)~gene_richness, data=Temp)
summary(temp_lm)

# final pub fig

p4 = ggplot(Temp, aes(x = gene_richness, y = log(temp))) + geom_point(aes(colour = phylum)) + theme_bw() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +scale_x_continuous(
  limits = c(0.0, 0.8), 
  breaks = seq(0.0, 0.8, 0.1), 
  expand = c(0, 0))+
  xlab("N2-Fixation Gene Richness") + ylab("Temperature") +
  scale_color_manual(values=palette_tax,
                     breaks=palette_breaks,
                     labels=palette_breaks,
                     name="phylum", na.value="black")
p4 = p4 + geom_abline(intercept=coef(log.temp.pagels)[1], slope=coef(log.temp.pagels)[2], linetype = "dashed", color = "blue" ) 
p4 = p4 + geom_abline(intercept=coef(temp_lm)[1], slope=coef(temp_lm)[2], linetype = "solid", color = "blue" ) 


p4


# combines all scatter plots into one figure

setwd("../../nif_phylogenomics/figures/scatters/")

figure <- ggarrange(p3, p2, p4, p1,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2, widths = c(2, 2, 2, 2))
figure
