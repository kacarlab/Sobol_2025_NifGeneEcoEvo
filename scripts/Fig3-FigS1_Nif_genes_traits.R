######### Nif-gene and functional group matrices mapped to phylogeny and trait data #########
# Last Edited 05/20/2025 by Morgan Sobol, adapted from Amanda Garcia

setwd("/figures/")

library(seecolor)
library(ape)
library(dplyr)
library(reshape2)
library(phytools)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(ggplotify)

# Import tree
tree <- read.newick("../GenomeData/genomeTrees/tree.nwk")

# Import metadata, genome label in 1st column = tree tip label
tree_metadata <- read.csv("../metadata/Nif_genome_metadata_final.csv")
View(tree_metadata)

taxonomy_table <- tree_metadata[, c("genome_ID", "phylum")]

# Replace empty taxonomy names with NA
phylum <- tree_metadata$phylum
phylum["phylum"][phylum["phylum"] == ""] <- NA
phylum <- as.data.frame(phylum)


# Prune tree, keep only those with nif 
tree <- drop.tip(tree,tree$tip.label[-match(tree_metadata$genome_ID, tree$tip.label)])

# get order of tree labels to be used later
species_tree_labels=as.data.frame(tree$tip.label)
write.csv(species_tree_labels, file="202502_ladderized_Zhu_tiplabels.csv")

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
                   "Epsilonproteobacteria","Spirochaetes","Synergistetes",
                   "Tenericutes","Thermodesulfobacteria","Thermotogae",
                   "Verrucomicrobia")



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
                "#DAD1DC","#AA856B","#F1F3B9",
                "#9CBBA2","#7496E0","#D7E0E5",
                "#AAA8AF")




# Plot first tree

tree1 <- ggtree(tree, ladderize=TRUE, right=TRUE, size=1)  
tree1                    

# Plot tree with taxonomy annotations
tree2 <- tree1 %<+% tree_metadata +
  geom_tippoint(aes(color=phylum), size=0.75) +
  scale_color_manual(values=palette_tax,
                     breaks=palette_breaks,
                     labels=palette_breaks,
                     name="Taxonomy", na.value="black") +
  geom_treescale(x=0,offset=3,fontsize=1,linesize=0.25) +
  theme(legend.position = "none")
tree2


# Load the presence/absence gene matrix
gene_matrix <- read.csv("../matrix/Nif_binary_matrix_final.csv", row.names=1)
gene_matrix <- as.data.frame(gene_matrix)
gene_matrix[sapply(gene_matrix, is.numeric)] <- lapply(gene_matrix[sapply(gene_matrix, is.numeric)], as.factor)
gene_matrix$genome_ID <- rownames(gene_matrix)  

# Convert the gene matrix to long format
gene_long <- melt(gene_matrix, id.vars = "genome_ID")
colnames(gene_long) <- c("genome_ID", "gene", "presence")

# Ensure presence is numeric
gene_long$presence <- as.numeric(gene_long$presence)

# Merge with taxonomy information
gene_long <- merge(gene_long, taxonomy_table, by = "genome_ID", all.x = TRUE)

# Assign colors: presence (1) gets phylum color, absence (0) remains grey96
gene_long$fill_color <- ifelse(gene_long$presence == "1", gene_long$phylum, "grey96")

# Convert back to wide format for gheatmap
gene_wide <- reshape2::dcast(gene_long, genome_ID ~ gene, value.var = "fill_color")

# Ensure row names match tree tips
rownames(gene_wide) <- gene_wide$genome_ID
gene_wide$genome_ID <- NULL  

palette_combo <- setNames(
  c("#864B61","#9A577E","#C78DAA", 
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
    "#DAD1DC","#AA856B","#F1F3B9",
    "#9CBBA2","#7496E0","#D7E0E5",
    "#AAA8AF"),
  c("Euryarchaeota","Crenarchaeota", "Thaumarchaeota",
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
    "Epsilonproteobacteria","Spirochaetes","Synergistetes",
    "Tenericutes","Thermodesulfobacteria","Thermotogae",
    "Verrucomicrobia")
)


# Define color mapping 
color_map <- c(palette_combo, "grey96" = "grey96") 

# Generate heatmap with correct presence/absence coloring
tree2.5 <- tree2 + new_scale_fill()
tree3 <- gheatmap(tree2.5, gene_wide, offset = 0.05, width = 3.5, 
                  colnames_angle = 90, colnames_offset_y = -30, font.size = 2, 
                  color = NA) +
  scale_fill_manual(values = color_map, name = "N2-fixation Gene \nPresence/Absence") +
  theme(legend.position = "none")

tree3

# extract trait data to plot with tree and matrix
temp <- data.frame(row.names=c(tree_metadata$genome_name), temperature=c(tree_metadata$temperature))
temp$temperature <- as.numeric(as.character(temp$temperature))

o2 <- data.frame(row.names=c(tree_metadata$genome_name), oxygen=c(tree_metadata$oxygen_discrete))

pH <- data.frame(row.names=c(tree_metadata$genome_name), pH=c(tree_metadata$pH))
tree_metadata$pH <- as.numeric(as.character(tree_metadata$pH))

habitat <- data.frame(row.names=c(tree_metadata$genome_name), habitat=c(tree_metadata$habitat))

genomesize <- data.frame(row.names=c(tree_metadata$genome_name), genomesize=c(tree_metadata$genome_size))


## add genome data
tree3.5 <- tree3 + new_scale_fill()
tree4 <- gheatmap(tree3.5, genomesize, offset=4.75, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2 ,color=NA) +  scale_fill_gradient(low = "white", high = "black", name="Genome Size") + theme(legend.position = "none")
tree4


## add temp data
tree4.5 <- tree4 + new_scale_fill()
tree5 <- gheatmap(tree4.5, temp, offset=4.85, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2 ,color=NA)  + scale_fill_gradient2(
  low = "#F5E596", 
  mid = "#DB7D1F", 
  high = "#B2182B", 
  midpoint = 50, name="Optimal Temperature") + theme(legend.position = "none")
tree5

## add oxygen data
tree5.5 <- tree5 + new_scale_fill()
tree6 <- gheatmap(tree5.5, o2, offset=4.95, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2, color=NA) +  scale_fill_manual(breaks=c("N/A", "aerobe", "anaerobe", "microaerophile", "facultative"), values=c("#f5f5f5", "#2b83ba", "#EE5C42", "#abdda4", "#fdae61"), name="Oxygen") + theme(legend.position = "none")
tree6

## add pH data
tree6.5 <- tree6 + new_scale_fill()
tree7 <- gheatmap(tree6.5, pH, offset=5.05, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2, color=NA) +  scale_fill_gradient2(
  low = "#DCE319FF", 
  mid = "#29AF7FFF", 
  high = "#404788FF", 
  midpoint = 7, name="pH") +theme(legend.position = "none")
tree7

## add habitat data 
tree7.5 <- tree7 + new_scale_fill()
tree8 <- gheatmap(tree7.5, habitat1, offset=5.15, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2, color=NA)+ scale_fill_manual(breaks=c("Freshwater","Saline","Hydrothermal", "Terrestrial", "Plant","Host-associated"), values=c("#75d9de","#3a6c6f", "#d1d388", "#e8b18d","#95dca9","#d4afdf"), name="Habitat") + theme(legend.position = "TRUE")
tree8 + expand_limits(y = -50)


############ Functional cluster matrix #####
# Load the relative abundance matrix
abundance_matrix <- read.csv("../matrix/Nif_clusters_matrix_final.csv", row.names=1)
abundance_matrix <- as.data.frame(abundance_matrix)

# Ensure values are numeric for scaling
abundance_matrix[] <- lapply(abundance_matrix, as.numeric)
abundance_matrix$genome_ID <- rownames(abundance_matrix)  # Add genome names as a column

# Convert the matrix to long format
abundance_long <- reshape2::melt(abundance_matrix, id.vars = "genome_ID")
colnames(abundance_long) <- c("genome_ID", "gene", "abundance")

# Convert back to wide format for gheatmap
abundance_wide <- reshape2::dcast(abundance_long, genome_ID ~ gene, value.var = "abundance")

# Ensure row names match tree tips
rownames(abundance_wide) <- abundance_wide$genome_ID
abundance_wide$genome_ID <- NULL  # Remove genome_ID column


library("RColorBrewer")


# Generate heatmap with relative abundance scaling (no phylum colors)
tree2.5 <- tree2 + new_scale_fill()
tree3 <- gheatmap(tree2.5, abundance_wide, offset = 0.0, width = 0.5, 
                  colnames_angle = 45, colnames_offset_y = -50, font.size = 2, 
                  color = NA) +
  scale_fill_gradient2(
    low = "#Eceff4", 
    mid = "#91A8D0", 
    high = "#2d4f80", 
    midpoint = 0.5, name="Gene Abundance")+  
  theme(legend.position = "right")

tree3


## add genome data
tree3.5 <- tree3 + new_scale_fill()
tree4 <- gheatmap(tree3.5, genomesize, offset=0.7, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2 ,color=NA) +  scale_fill_gradient(low = "white", high = "black", name="Genome Size") + theme(legend.position = "none")
tree4

## add temp data
tree4.5 <- tree4 + new_scale_fill()
tree5 <- gheatmap(tree4.5, temp, offset=0.785, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2 ,color=NA) + scale_fill_gradient2(
  low = "#F5E596", 
  mid = "#DB7D1F", 
  high = "#B2182B", 
  midpoint = 50, name="Optimal Temperature") +  theme(legend.position = "none")
tree5

## add oxygen data
tree5.5 <- tree5 + new_scale_fill()
tree6 <- gheatmap(tree5.5, o2, offset=0.87, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2, color=NA) +  scale_fill_manual(breaks=c("N/A", "aerobe", "anaerobe", "microaerophile", "facultative"), values=c("#f5f5f5", "#2b83ba", "#EE5C42", "#abdda4", "#fdae61"), name="Oxygen") + theme(legend.position = "none")
tree6

## add pH data
tree6.5 <- tree6 + new_scale_fill()
tree7 <- gheatmap(tree6.5, pH, offset=0.955, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2, color=NA) +  scale_fill_gradient2(
  low = "#DCE319FF", 
  mid = "#29AF7FFF", 
  high = "#404788FF", 
  midpoint = 7, name="pH") +theme(legend.position = "none")
tree7


## add habitat data 
tree7.5 <- tree7 + new_scale_fill()
tree8 <- gheatmap(tree7.5, habitat1, offset=1.04, width=.05,colnames_angle=90, colnames_offset_y = -30, font.size = 2, color=NA) + scale_fill_manual(breaks=c("Freshwater","Saline","Hydrothermal", "Terrestrial", "Plant","Host-associated"), values=c("#75d9de","#3a6c6f", "#d1d388", "#e8b18d","#95dca9","#d4afdf"), name="Habitat") + theme(legend.position = "TRUE")
tree8 + expand_limits(y = -50)

