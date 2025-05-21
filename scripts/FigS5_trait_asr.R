##### Nif Trait ASR #####
## Last edited 5/20/2025 by Morgan Sobol

setwd("/ASR/")

library("ape")
library("phangorn")
library("phytools")
library("geiger")


# Import tree
tree <- read.newick("../GenomeData/genomeTrees/tree.nwk")

# Import metadata, genome label in 1st column = tree tip label
tree_metadata <- read.csv("../metadata/Nif_genome_metadata_final.csv")

# Prune tree, keep only those with nif 
tree <- drop.tip(tree,tree$tip.label[-match(tree_metadata$genome_ID, tree$tip.label)])

tree1 <- ladderize(tree, right=TRUE)
tree1 <-multi2di(tree1, random=FALSE)
tree1_ultra <- chronos(tree1, model="relaxed")
is.ultrametric((tree1_ultra))
is.binary(tree1)

plotTree(tree1)
nodelabels(frame="none")
tiplabels()

#########  habitat coverage asr ##################

hab_df = tree_metadata[, c("genome_ID", "habitat_coverage")]

hab_df = hab_df[!is.na(hab_df$habitat_coverage), ]

tree2 <-drop.tip(tree1,tree1$tip.label[-match(hab_df$genome_ID, tree1_ultra$tip.label)])

# make sure hab data and tree tips are same order 

order <- cbind(tree2$tip.label)
order <- as.data.frame(order, seq=seq(1:length(order)))

hab_df <- hab_df[order(match(hab_df$genome_ID, tree2$tip.label)), ]

rownames(hab_df) = hab_df$genome_ID
hab_df$genome_ID = NULL

hab<-as.matrix(hab_df)[,1]


########

fit<-fastAnc(tree2,hab,vars=TRUE,CI=TRUE)
fit


# Since true ancestral states are unknown, plot tip values against immediate ancestor estimates
ancestors <- sapply(1:length(hab), function(i) getParent(tree2, i))

# Compute correlation of predicted to true values
cor.test(hab, fit$ace[ancestors - length(tree2$tip.label)])

plot(hab, fit$ace[ancestors - length(tree2$tip.label)],
     xlab="Observed Habitat Values at Tips",
     ylab="Estimated Ancestral Habitat Values",
     main="Strong Correlation between Observed and Ancestral Values")
abline(a=0, b=1, col="red", lty=2) 


# Check confidence estimates
estimates <- fit$ace
CI95 <- fit$CI95  

# Get node numbers for tips immediate ancestors
ancestors <- sapply(1:length(tree2$tip.label), function(x) {
  Ancestors(tree2, x, type = "parent")
})

observed_in_CI <- sapply(1:length(hab), function(i) {
  ancestor_node <- ancestors[i]
  observed_value <- hab[i]
  ci_lower <- CI95[ancestor_node - length(tree2$tip.label), 1]
  ci_upper <- CI95[ancestor_node - length(tree2$tip.label), 2]
  
  observed_value >= ci_lower & observed_value <= ci_upper
})

# Determine the proportion of data fitting within the ancestral CIs
mean(observed_in_CI)


#plot 

obj<-contMap(tree2,hab,plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(tree2)), ftype="off", lwd=1)

plot(obj, ftype = "off", leg.txt = "Habitat Coverage", lwd = 1)

viridis.cMap<-setMap(obj,viridisLite::viridis(n=8))
plot(viridis.cMap,ftype="off",lwd = 1.5,outline=FALSE, leg.txt="Habitat Coverage")



############# gene richness asr ###########
gene_df = tree_metadata[, c("genome_ID", "gene_richness")]
gene_df = gene_df[!is.na(gene_df$gene_richness), ]

tree2 <-drop.tip(tree1,tree1$tip.label[-match(gene_df$genome_ID, tree1$tip.label)])

order <- cbind(tree2$tip.label)
order <- as.data.frame(order, seq=seq(1:length(order)))

gene_df <- gene_df[order(match(gene_df$genome_ID, tree2$tip.label)), ]

rownames(gene_df) = gene_df$genome_ID
gene_df$genome_ID = NULL

gene<-as.matrix(gene_df)[,1]


####

fit_u<-fastAnc(tree2,gene,vars=TRUE,CI=TRUE)
fit_u

# Since true ancestral states are unknown, plot tip values against immediate ancestor estimates
ancestors <- sapply(1:length(gene), function(i) getParent(tree2, i))

# Compute correlation
cor.test(gene, fit_u$ace[ancestors - length(tree2$tip.label)])


plot(gene, fit_u$ace[ancestors - length(tree2$tip.label)],
     xlab="Observed Metabolism Values at Tips",
     ylab="Estimated Ancestral Metabolism Values",
     main="Strong Correlation between Observed and Ancestral Values")
abline(a=0, b=1, col="red", lty=2) # 1:1 line


# Check confidence estimates
estimates <- fit_u$ace
CI95 <- fit_u$CI95  

# Get node numbers for tips’ immediate ancestors
ancestors <- sapply(1:length(tree2$tip.label), function(x) {
  Ancestors(tree2, x, type = "parent")
})

observed_in_CI <- sapply(1:length(gene), function(i) {
  ancestor_node <- ancestors[i]
  observed_value <- gene[i]
  ci_lower <- CI95[ancestor_node - length(tree2$tip.label), 1]
  ci_upper <- CI95[ancestor_node - length(tree2$tip.label), 2]
  
  observed_value >= ci_lower & observed_value <= ci_upper
})

# Proportion of data fitting within the ancestral CIs
mean(observed_in_CI)

#plot 

obj<-contMap(tree2,gene,plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(tree2)), ftype="off", lwd=1)

viridis.cMap<-setMap(obj,viridisLite::viridis(n=8))
plot(viridis.cMap,ftype="off",lwd = 1.5,outline=FALSE, leg.txt="Gene Richness")

