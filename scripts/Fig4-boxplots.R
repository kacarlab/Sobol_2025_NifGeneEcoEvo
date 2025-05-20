######### Boxplots and stats for Nif genes #########
# Last Edited 05/20/2025 by Morgan Sobol

setwd("/boxplots")

library(ggpubr)
library(ggplot2)
library(ggnewscale)
library(ggplotify)

## Import metadata

tree_metadata <- read.csv("../../metadata/Nif_genome_metadata_final.csv")
View(tree_metadata)

## Extract columns for each trait

O2_extract <- c("gene_richness", "oxygen", "phylum") 
O2vGene <- tree_metadata[, O2_extract, drop = FALSE]
#remove NAs
O2vGene <- O2vGene[complete.cases(O2vGene), ]

hab_extract <- c("gene_richness", "habitat", "phylum") 
HabitatvGene <- tree_metadata[, hab_extract, drop = FALSE]

temp_extract <- c("gene_richness", "temp_discrete", "phylum") 
TempvGene <- tree_metadata[, temp_extract, drop=FALSE]

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


#### Oxygen ################################################

# First remap actual oxygen categories to placeholder names
oxygen_map <- c("anaerobe" = "Group1",
                "microaerophile" = "Group2",
                "facultative" = "Group3",
                "aerobe" = "Group4")

O2vGene$oxygen_group <- oxygen_map[as.character(O2vGene$oxygen)]

# Add 2 dummy groups to reach 6 total
dummy_rows <- data.frame(
  oxygen_group = c("Group5", "Group6"),
  gene_richness = NA, 
  phylum = NA
)

oxygen_padded <- rbind(O2vGene[, c("oxygen_group", "gene_richness", "phylum")], dummy_rows)

# Set correct factor levels
oxygen_padded$oxygen_group <- factor(oxygen_padded$oxygen_group,
                                     levels = c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6"))


# Perform pairwise comparisons
compare_means(gene_richness ~ oxygen,  data = O2vGene)

p <- ggboxplot(O2vGene, x = "oxygen", y = "gene_richness",add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")


# Specify the comparisons you want
my_comparisons <- list( c("anaerobe", "aerobe"), c("anaerobe", "facultative"), c("anaerobe", "microaerophile"))

oxygen_labels <- c("anaerobe", "microaerophile", "facultative", "aerobe", "", "")

o2_box <- ggboxplot(oxygen_padded, x = "oxygen_group", y = "gene_richness", na.rm = TRUE) + geom_point(aes(color = phylum),
             position = position_jitter(width = 0.2),
             size = 2, na.rm = TRUE)+
  scale_color_manual(values=palette_tax,
                     breaks=palette_breaks,
                     labels=palette_breaks,
                     name="phylum")+
  ylab("N2 fixation Gene Richness") +
  scale_x_discrete(labels = oxygen_labels, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "none")

o2_box

#### Habitat ################################################

# Perform pairwise comparisons
compare_means(gene_richness ~ habitat,  data = HabitatvGene)


p <- ggboxplot(HabitatvGene, x = "habitat", y = "gene_richness",add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")

#Specify the comparisons you want
my_comparisons <- list(c("Hydrothermal", "Freshwater"),c("Hydrothermal", "Terrestrial"), c("Hydrothermal", "Saline"),c("Hydrothermal", "Host-associated"), c("Hydrothermal", "Plant"), c("Freshwater", "Saline"), c("Freshwater", "Terrestrial"), c("Freshwater","Host-associated"), c("Freshwater", "Plant"))

habitat_box=ggboxplot(HabitatvGene, x = "habitat", y = "gene_richness", color = "habitat" ,add="jitter")+ xlab("") + ylab("N2-fixation Gene Richness")  + stat_compare_means(label.y = 1.5, label.x=0.75) +theme_classic()+ ylim(0,1.0)+ theme(legend.position="none") 
habitat_box=habitat_box+stat_compare_means(comparisons = my_comparisons, aes(label = ..p.signif..))
habitat_box

habitat_box <- ggboxplot(HabitatvGene, x = "habitat", y = "gene_richness", width = 0.8, na.rm = TRUE) +
  geom_point(aes(color = phylum),
               position = position_jitter(width = 0.2),
               size = 2, na.rm = TRUE) +
  scale_color_manual(values=palette_tax,
                     breaks=palette_breaks,
                     labels=palette_breaks,
                     name="phylum")+
  ylab("N2-fixation Gene Richness") +
  scale_x_discrete(labels = HabitatvGene$habitat, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "none")

habitat_box


#### Temp ################################################
temp_map <- c("Mesophile" = "Group1",
              "Thermophile" = "Group2")

TempvGene$temp_group <- temp_map[as.character(TempvGene$temp_discrete)]

dummy_temp <- data.frame(
  temp_group = c("Group3", "Group4", "Group5", "Group6"),
  gene_richness = NA,
  phylum = NA
)

temp_padded <- rbind(TempvGene[, c("temp_group", "gene_richness", "phylum")], dummy_temp)

temp_padded$temp_group <- factor(temp_padded$temp_group,
                                 levels = c("Group1", "Group2", "Group3", "Group4", "Group5", "Group6"))


# Perform pairwise comparisons
compare_means(gene_richness ~ temp_discrete,  data = TempvGene)


p <- ggboxplot(TempvGene, x = "temp_discrete", y = "gene_richness",add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")


# Specify the comparisons you want
my_comparisons <- list(c("mesophile", "thermophile"))


temp_labels <- c("Mesophile", "Thermophile", "", "", "", "")

temp_box <- ggboxplot(temp_padded, x = "temp_group", y = "gene_richness", width = 0.8, na.rm = TRUE) +
  geom_point(aes(color = phylum),
             position = position_jitter(width = 0.2),
             size = 2, na.rm = TRUE) +
  scale_color_manual(values=palette_tax,
                     breaks=palette_breaks,
                     labels=palette_breaks,
                     name="phylum")+
  ylab("N2-fixation Gene Richness") +
  scale_x_discrete(labels = temp_labels, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "none")

temp_box


## plot all together

figure <- ggarrange(temp_box, o2_box, habitat_box,
                    labels = c("A", "B", "C"),
                    ncol = 3, widths = c(2, 2, 2))
figure
    



