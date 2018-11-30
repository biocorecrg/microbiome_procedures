# analysis_16S_v1.R
# By: Julia Ponomarenko, CRG
# Date: Nov 28, 2018
# This is an interactive script to play with data rather than a finalized one.
# Therefore, run it in R Studio.

# https://rpubs.com/dillmcfarlan/R_microbiotaSOP
# https://f1000research.com/articles/5-1492/v2

#Updated R to 3.5.1 (2018-07-02) and R Studio to the latest.  
#Got phyloseq 1.24.2 and it was installed smoothly. 
#With old version of R only old version of phyloseq can be installed, which didn't work.

# install microbiome
library(BiocInstaller)
source("http://www.bioconductor.org/biocLite.R")
biocLite("microbiome")
####

library(phyloseq)
library(microbiome)

library(caret)
library(reshape2)
library(ggplot2)
library(dplyr) 
library(ape)
library(gplots)
library(Matrix)
library(lme4)
library(phangorn)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(plyr)
library(RColorBrewer)

########## Get input data #################

#OTU table (shared file)
df = read.table("16S_OTU_shared.txt", header=F, sep="\t", as.is = TRUE)
df <- df[,-c(1,3)]

df <- t(df)

d <- df
colnames(d) <- d[1,]
d <- d[-1,]
rownames(d) <- d[,1]
d <- d[,-1]
d <- data.frame(d, stringsAsFactors=F)
d1 = data.frame(lapply(d, function(x) as.numeric(x)),
                check.names=F, row.names = rownames(d)) # to preserve column names
otumat <- as.matrix(d1)
class(otumat)

#Taxonomy of each OTU
df = read.table("16S_OTU_taxonomy.txt", header=TRUE, sep="\t", as.is = T)
df[df$OTU == "Otu002",]
df <- as.data.frame(sapply(df, function(x) {gsub("\\s*\\([^\\)]+\\)", "", x)}))
df[df$OTU == "Otu002",]

#genus is the lowest taxa level for 16S
df <- separate(df, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep=";")
df[df$OTU == "Otu002",]
row.names(df) <- df[,1]
df <- df[,-c(1,2)]
taxmat <- as.matrix(df)

# metadata with diversity indexes
meta = read.table("16S_metadata.txt", header=TRUE, sep="\t", as.is = T)
meta <- meta[,-1]

meta$group <- substr(meta$Sample, 1, 1)

meta <- cbind.data.frame(Sample=meta$group, Group = "", meta[,-1])
meta$Group <- substr(meta$Sample, 1, 1)
row.names(meta) <- meta[,1]
meta <- meta[,-1]
class(meta)
meta[rownames(meta) == "D6306",]$Group <- "D6306"
meta[rownames(meta) == "HM782",]$Group <- "HM782"
meta[rownames(meta) == "HM783",]$Group <- "HM783"

# read the tree
tree <- read_tree("16S_OTU.tre")
#plot(tree)

#Create phyloseq object
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)

sample_data <- sample_data(meta, errorIfNULL=TRUE)
ps = merge_phyloseq(physeq, sample_data, tree)
ps
######### End of input data #############

### Collapse OTUs by similar genus (the lowest taxa level for 16S)
ntaxa(ps)
ps2 <- tax_glom(ps, taxrank="Genus")
ntaxa(ps2); tax_table(ps2)[1:3, c("Phylum", "Class", "Order", "Family", "Genus")]
otumat <- otu_table(ps2)
Object <- ps2

##### Taxonomic filtering (see https://f1000research.com/articles/5-1492/v2)
#
#One of the reasons to filter in this way is to avoid spending much time 
#analyzing taxa that were seen only rarely among samples. This also turns 
# out to be a useful filter of noise (taxa that are actually just artifacts
#of the data collection process), a step that should probably be considered 
#essential for datasets constructed via heuristic OTU-clustering methods, 
#which are notoriously prone to generating spurious taxa.
# Considering that this is really a noise, we just remove those taxa, instead of aggregating them to other.

# Compute prevalence of each phylum, which is defined as 
# the number of samples in which the phylum appears
prevdf = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf, TotalCounts = taxa_sums(ps2), tax_table(ps2))
plyr::ddply(prevdf, "Phylum", function(df1){
  cbind(mean(df1$Prevalence),sum(df1$Prevalence), sum(df1$TotalCounts))})
# It can be seen that some phyla appear only in one or few samples
# and some have very low total counts

# Let's remove those phyla manually
remove_phyla <- c("candidate_division_WPS-2","Elusimicrobia","Latescibacteria","SR1","Tenericutes",
                  "Bacteria_unclassified")
ps2 <- subset_taxa(ps2, !Phylum %in% remove_phyla)

# Which phyla?
unique(factor(tax_table(ps2)[, "Phylum"]))
Object <- ps2

###### END of filtering taxa #######



###### Explore MOCK samples ############

ps2 <- Object
pmock <- subset_samples(ps2, sample_names(ps2) %in% c("D6306","HM782","HM783"))

# filter rare taxa
ps2 <- pmock
prevdf = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf, TotalCounts = taxa_sums(ps2), tax_table(ps2))
plyr::ddply(prevdf, "Phylum", function(df1){
  cbind(mean(df1$Prevalence),sum(df1$Prevalence), sum(df1$TotalCounts))})
keepTaxa = rownames(prevdf)[(prevdf$TotalCounts >= 100)]
ps2 = prune_taxa(keepTaxa, ps2)
tax_table(ps2)
otu_table(ps2)

# make a table of genera and counts 

fac = factor(tax_table(ps2)[, "Genus"])
tab = apply(otu_table(ps2), MARGIN = 2, function(x) {
  tapply(x, INDEX = fac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
tab
df <- as.data.frame(tab)
df <- cbind("Genus" = rownames(df),df)
Family <- microbiome::map_levels(as.vector(df$Genus), "Genus", "Family",  ps2)
df <- cbind(Family,df)
Order <- microbiome::map_levels(as.vector(df$Family), "Family", "Order", ps2)
df <- cbind(Order,df)
Class <- microbiome::map_levels(as.vector(df$Order), "Order", "Class", ps2)
df <- cbind(Class,df)
Phylum <- microbiome::map_levels(as.vector(df$Class), from="Class", to="Phylum", ps2)
df <- cbind(Phylum,df)
df <- df[order(df[,1], df[,2], df[,3], df[,4], df[,5]), ]
head(df)

file_out <- "16S_counts_MOCK_samples.txt"
write.table(df, file_out, quote=F, sep="\t", row.names=F) #, col.names=NA)

## plot rel abundances at genus level
glom <- tax_glom(ps2, taxrank = 'Genus')
dat <- psmelt(glom)
dat$Genus <- as.character(dat$Genus)
p1 <- ggplot(dat, aes(x=Genus, y=Abundance)) + geom_boxplot() + coord_flip() 
c_count = length(unique(dat$Genus))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
p2 <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Genus))
p2 <- p2 + geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values=getPalette(c_count)) + 
  theme(legend.position="bottom") + 
  guides(fill=guide_legend(nrow=5))
print(p1); print(p2)

pdf("plots_rel_abundances_MOCK_samples.pdf", onefile = TRUE, width = 8, height = 6 ) # size in cm
print(p2)
dev.off()


###### END of exploring MOCK samples ###################

#### Filter taxa with MOCK samples removed ###########

## Remove MOCK samples
ps2 <- Object
sample_names(ps2)
ps2 <- subset_samples(ps2, !sample_names(ps2) %in% c("D6306","HM782","HM783"))
sample_names(ps2)


prevdf = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf, TotalCounts = taxa_sums(ps2), tax_table(ps2))
plyr::ddply(prevdf, "Phylum", function(df1){
  cbind(mean(df1$Prevalence),sum(df1$Prevalence), sum(df1$TotalCounts))})
# It can be seen that some phyla appear only in one or few samples
# and some have very low total counts

# look at phyla prevalence versus total abundance
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
ggplot(prevdf1, aes(TotalCounts, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Let's remove those phyla manually
remove_phyla <- c("Nitrospinae")
ps2 <- subset_taxa(ps2, !Phylum %in% remove_phyla)

# Which phyla?
unique(factor(tax_table(ps2)[, "Phylum"]))

# Let's remove OTUs that present only in 1 sample
#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold 

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)

# Let's check how many phyla and features we have now
ps2
length(get_taxa_unique(ps2, taxonomic.rank = "Phylum"))

# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))

###### END of filtering taxa #######

###### Make files with counts by taxa levels ####

ps2 <- Object
# Create a factor corresponding to the taxa level
fac = factor(tax_table(ps2)[, "Phylum"])
# Tabulate the counts for each taxa in each sample
tab = apply(otu_table(ps2), MARGIN = 2, function(x) {
  tapply(x, INDEX = fac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(tab)[, 1:10]

df <- as.data.frame(tab)
df <- cbind(rownames(df),df)
head(df)[, 1:10]
colnames(df)[1] <- "Phylum"
file_out <- "counts_Phylum.txt"
write.table(df, file_out, quote=F, sep="\t", row.names=F) #, col.names=NA)


## Class level ####
fac = factor(tax_table(ps2)[, "Class"])
tab = apply(otu_table(ps2), MARGIN = 2, function(x) {
  tapply(x, INDEX = fac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(tab)[, 1:10]
df <- as.data.frame(tab)
df <- cbind("Class" = rownames(df),df)
head(df)[, 1:10]
Phylum <- microbiome::map_levels(as.vector(df$Class), from="Class", to="Phylum", ps2)
df <- cbind(Phylum,df)
df <- df[order(df[,1], df[,2]), ]

file_out <- "counts_Class.txt"
write.table(df, file_out, quote=F, sep="\t", row.names=F) #, col.names=NA)


## Order level ####
fac = factor(tax_table(ps2)[, "Order"])
tab = apply(otu_table(ps2), MARGIN = 2, function(x) {
  tapply(x, INDEX = fac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(tab)[, 1:10]
df <- as.data.frame(tab)
df <- cbind("Order" = rownames(df),df)
head(df)[, 1:10]
Class <- microbiome::map_levels(as.vector(df$Order), "Order", "Class", ps2)
df <- cbind(Class,df)
head(df)[, 1:10]
Phylum <- microbiome::map_levels(as.vector(df$Class), from="Class", to="Phylum", ps2)
df <- cbind(Phylum,df)
df <- df[order(df[,1], df[,2], df[,3]), ]
head(df)[, 1:10]

file_out <- "counts_Order.txt"
write.table(df, file_out, quote=F, sep="\t", row.names=F) #, col.names=NA)


## Family level ####
fac = factor(tax_table(ps2)[, "Family"])
tab = apply(otu_table(ps2), MARGIN = 2, function(x) {
  tapply(x, INDEX = fac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(tab)[, 1:10]
df <- as.data.frame(tab)
df <- cbind("Family" = rownames(df),df)
head(df)[, 1:10]
Order <- microbiome::map_levels(as.vector(df$Family), "Family", "Order", ps2)
df <- cbind(Order,df)
head(df)[, 1:10]
Class <- microbiome::map_levels(as.vector(df$Order), "Order", "Class", ps2)
df <- cbind(Class,df)
head(df)[, 1:10]
Phylum <- microbiome::map_levels(as.vector(df$Class), from="Class", to="Phylum", ps2)
df <- cbind(Phylum,df)
df <- df[order(df[,1], df[,2], df[,3], df[,4]), ]
head(df)[, 1:10]

file_out <- "counts_Family.txt"
write.table(df, file_out, quote=F, sep="\t", row.names=F) #, col.names=NA)


## Genus level ####
fac = factor(tax_table(ps2)[, "Genus"])
tab = apply(otu_table(ps2), MARGIN = 2, function(x) {
  tapply(x, INDEX = fac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(tab)[, 1:10]
df <- as.data.frame(tab)
df <- cbind("Genus" = rownames(df),df)
head(df)[, 1:10]
Family <- microbiome::map_levels(as.vector(df$Genus), "Genus", "Family",  ps2)
df <- cbind(Family,df)
head(df)[, 1:10]
Order <- microbiome::map_levels(as.vector(df$Family), "Family", "Order", ps2)
df <- cbind(Order,df)
head(df)[, 1:10]
Class <- microbiome::map_levels(as.vector(df$Order), "Order", "Class", ps2)
df <- cbind(Class,df)
head(df)[, 1:10]
Phylum <- microbiome::map_levels(as.vector(df$Class), from="Class", to="Phylum", ps2)
df <- cbind(Phylum,df)
df <- df[order(df[,1], df[,2], df[,3], df[,4], df[,5]), ]
head(df)[, 1:10]

file_out <- "counts_Genus.txt"
write.table(df, file_out, quote=F, sep="\t", row.names=F) #, col.names=NA)

###### END of making count files #########


###### Plot relative abundances ########

ps3 <- ps2

# transform counts to relative abundances
X <- ps3
ps2_rel = transform_sample_counts(X, function(x){x / sum(x)})

#### let's plot relative abundances for each sample by phylum

# agglomerate taxa
glom <- tax_glom(ps2_rel, taxrank = 'Phylum')

# create dataframe from phyloseq object
dat <- psmelt(glom)

# convert Phylum to a character vector from a factor
dat$Phylum <- as.character(dat$Phylum)

# group dataframe by Phylum, calculate median rel. abundance
medians <- plyr::ddply(dat, ~Phylum, function(x) c(median=median(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
rare <- medians[medians$median <= 0.01,]$Phylum
# change their name to "Other"
dat[dat$Phylum %in% rare,]$Phylum <- 'Other'

# boxplot
p1 <- ggplot(dat, aes(x=Phylum, y=Abundance)) + geom_boxplot() + coord_flip()

#bar-plot
c_count = length(unique(dat$Phylum))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(c_count, "Paired"))

p2 <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Phylum))
p2 <- p2 + geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values=getPalette(c_count)) + theme(legend.position="bottom") + 
  guides(fill=guide_legend(nrow=5))


pdf("plots_rel_abundances_Phylum_level_collapsed.pdf", onefile = TRUE, width = 12, height = 6 ) # size in cm
print(p1); print(p2)
dev.off()

#### let's plot relative abundances for each phylum

x <- unique(dat$Phylum)
x <- x[x != "Other"]
phyla_list <- x

for (ph in phyla_list){
  print(ph)
  file_pdf = paste("plots_rel_abundances_Phylum_",ph,".pdf",sep="")
  pdf(file_pdf, onefile = TRUE, width = 12, height = 6 ) # size in cm
  #subset by Phylum
  xp <- subset_taxa(ps2_rel, Phylum == ph)
  glom <- tax_glom(xp, taxrank = 'Class')
  dat <- psmelt(glom)
  dat$Class <- as.character(dat$Class)
  print(length(unique(dat$Class)))
  #if (length(unique(dat$Class)) < 2) next; #nothing to plot for just one Class!
  # if (length(unique(dat$Class)) > 20) {
  #   medians <- plyr::ddply(dat, ~Class, function(x) c(median=median(x$Abundance)))
  #   rare <- medians[medians$median <= 0.01,]$Class
  #   if (length(rare) > 1) dat[dat$Class %in% rare,]$Class <- 'Other'
  # }  
  
  p1 <- ggplot(dat, aes(x=Class, y=Abundance)) + geom_boxplot() + coord_flip() + ggtitle(paste("Phylum = ", ph))
  c_count = length(unique(dat$Class))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  p2 <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Class))
  p2 <- p2 + geom_bar(aes(), stat="identity", position="stack") + 
    ggtitle(paste("Phylum = ", ph)) +
    scale_fill_manual(values=getPalette(c_count)) + theme(legend.position="bottom") + 
    guides(fill=guide_legend(nrow=5))
  print(p1); print(p2)

  # plot at Order level for each Class that other than unclassified
  x <- unique(dat$Class)
  x <- x[!grepl("unclassified", x)]
  x <- x[x != "Other"]
  class_list <- x
  for (cl in class_list){
    print(cl)
    #subset by Phylum
    xc <- subset_taxa(ps2_rel, Class == cl)
    glom <- tax_glom(xc, taxrank = 'Order')
    dat <- psmelt(glom)
    dat$Order <- as.character(dat$Order)
    print(length(unique(dat$Order)))
    #if (length(unique(dat$Order)) < 2) next; #nothing to plot for just one Class!
    # if (length(unique(dat$Order)) >20) {
    #   medians <- plyr::ddply(dat, ~Order, function(x) c(median=median(x$Abundance)))
    #   rare <- medians[medians$median <= 0.001,]$Order
    #   if (length(rare) > 1) dat[dat$Order %in% rare,]$Order <- 'Other'
    # }

    p1 <- ggplot(dat, aes(x=Order, y=Abundance)) + geom_boxplot() + coord_flip() + 
      ggtitle(paste("Phylum = ", ph, "; Class = ", cl))
    c_count = length(unique(dat$Order))
    getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    p2 <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Order))
    p2 <- p2 + geom_bar(aes(), stat="identity", position="stack") + 
      ggtitle(paste("Phylum = ", ph, "; Class = ", cl)) +
      scale_fill_manual(values=getPalette(c_count)) + theme(legend.position="bottom") + 
      guides(fill=guide_legend(nrow=5))
    print(p1); print(p2)

    # plot at Family level for each Order that other than unclassified
    x <- unique(dat$Order)
    x <- x[!grepl("unclassified", x)]
    x <- x[x != "Other"]
    order_list <- x
    for (od in order_list){
      print(od)
      #subset by Phylum
      xc <- subset_taxa(ps2_rel, Order == od)
      glom <- tax_glom(xc, taxrank = 'Family')
      dat <- psmelt(glom)
      dat$Family <- as.character(dat$Family)
      print(length(unique(dat$Family)))
      #if (length(unique(dat$Family)) < 2) next; #nothing to plot for just one Class!
      # if (length(unique(dat$Family)) >20) {
      #   medians <- plyr::ddply(dat, ~Family, function(x) c(median=median(x$Abundance)))
      #   rare <- medians[medians$median <= 0.001,]$Family
      #   if (length(rare) > 1) dat[dat$Family %in% rare,]$Family <- 'Other'
      # }
      p1 <- ggplot(dat, aes(x=Family, y=Abundance)) + geom_boxplot() + coord_flip() + 
        ggtitle(paste("Phylum = ", ph, "; Class = ", cl, "; Order = ", od))
      c_count = length(unique(dat$Family))
      getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
      p2 <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Family))
      p2 <- p2 + geom_bar(aes(), stat="identity", position="stack") + 
        ggtitle(paste("Phylum = ", ph, "; Class = ", cl, "; Order = ", od)) +
        scale_fill_manual(values=getPalette(c_count)) + 
        theme(legend.position="bottom") + 
        guides(fill=guide_legend(nrow=5))
      print(p1); print(p2)

      # plot at Genus level for each Family that other than unclassified
      x <- unique(dat$Family)
      x <- x[!grepl("unclassified", x)]
      x <- x[x != "Other"]
      family_list <- x
      for (fm in family_list){
        print(fm)
        xc <- subset_taxa(ps2_rel, Family == fm)
        glom <- tax_glom(xc, taxrank = 'Genus')
        dat <- psmelt(glom)
        dat$Genus <- as.character(dat$Genus)
        #if (length(unique(dat$Genus)) < 2) next; #nothing to plot for just one Class!
        # if (length(unique(dat$Genus)) >20) {
        #   medians <- plyr::ddply(dat, ~Genus, function(x) c(median=median(x$Abundance)))
        #   rare <- medians[medians$median <= 0.001,]$Genus
        #   if (length(rare) > 1) dat[dat$Genus %in% rare,]$Genus <- 'Other'
        # }
        p1 <- ggplot(dat, aes(x=Genus, y=Abundance)) + geom_boxplot() + coord_flip() + 
          ggtitle(paste("Phylum = ", ph, "; Class = ", cl, "; Order = ", od, "; Family = ", fm))
        c_count = length(unique(dat$Genus))
        getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
        p2 <- ggplot(data=dat, aes(x=Sample, y=Abundance, fill=Genus))
        p2 <- p2 + geom_bar(aes(), stat="identity", position="stack") + 
          ggtitle(paste("Phylum = ", ph, "; Class = ", cl, "; Order = ", od, "; Family = ", fm)) +
          scale_fill_manual(values=getPalette(c_count)) + 
          theme(legend.position="bottom") + 
          guides(fill=guide_legend(nrow=5))
        print(p1); print(p2)
      }
    }
  }
  dev.off()
}
######## END of plotting rel abundances ###########



########## Calculate beta-diversity ##############
# Must use counts!!!! (If needed, phyloseq converts to rel abundances)
#
# Beta-diversity shows how different every sample is from every other sample. 
# Some metrics take abundance into account (i.e. diversity: Bray-Curtis,
# weighted UniFrac) and some only calculate based on presence-absence 
# (i.e. richness: Jaccard, unweighted UniFrac).
# https://joey711.github.io/phyloseq/distance.html

theme_set(theme_bw()) # theme for ggplot with white background

dist_methods <- unlist(distanceMethodList)
print(dist_methods)
dist_methods <- c("unifrac", "wunifrac", "jsd", "bray", "canberra", "jaccard")

plist <- vector("list", length(dist_methods))
object <- ps3

level <- "Phylum"
object <- tax_glom(object, taxrank = level)
tax_table(object)[1:3, c("Phylum", "Class", "Order", "Family", "Genus")]

for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(object, method=i)
  #save it to the file
  x <- as.matrix(iDist)
  file <- paste("dist_",i, "_level_", level, ".txt", sep="")
  write.table(x, file, row.names = T, sep = "\t", quote = FALSE)
  
  # Calculate ordination
  iMDS  <- ordinate(object, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(object, iMDS, color="Group")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
length(unique(meta$Group)) # this is how many colors I need and where to place black or whatever color
my14colors <- c("red","brown","blue","grey50","magenta","cornflowerblue","cyan","green","black",
                "forestgreen","darkorange","bisque","gold", "pink")

df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p <- ggplot(df, aes(Axis.1, Axis.2, color=Group))
p <- p + geom_point(size=2, alpha=1.0) + facet_wrap(~distance, scales="free") + 
  scale_colour_manual(values=my14colors) +
  ggtitle(paste("MDS on various distance metrics. Taxonomic level = ", level, sep=""))
print(p)

file <- paste("plots_beta_diversity_level_",level,".pdf", sep="")
pdf(file, onefile = TRUE, width = 10, height = 7 ) # size in cm
print(p)
dev.off()

########## Explore alpha-diversity metrics ##########
# http://joey711.github.io/phyloseq/plot_richness-examples

pdf("plots_alpha_diversity.pdf", onefile = TRUE, width = 8, height = 6 ) # size in cm
#plot_richness(ps, x="Group", color="Group", measures=c("Observed"))
plot_richness(ps, x="Group", measures=c("Observed"))
plot_richness(ps, x="Group", measures=c("Chao1"))
plot_richness(ps, x="Group", measures=c("ACE"))
plot_richness(ps, x="Group", measures=c("Simpson"))
plot_richness(ps, x="Group", measures=c("InvSimpson"))
plot_richness(ps, x="Group", measures=c("Fisher"))
dev.off()



######### Statistical analysis #################

# Alpha-diversity: check the distributions 
par(mfrow = c(2, 2))
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(1/meta$simpson, main="Inverse Simpson diversity", xlab="", breaks=10)
hist(meta$chao, main="Chao richness", xlab="", breaks=15)
hist(meta$ace, main="ACE richness", xlab="", breaks=15)

#check these distribution on normality
shapiro.test(meta$shannon)
shapiro.test(1/meta$simpson)
shapiro.test(meta$chao)
shapiro.test(meta$ace)

# It can be seen that chao and ace are normally distributed, 
# therefore t-test can be run on comparing groups of samples by these metrics
# While for other, Kruskal-Wallis or Wilcoxon rank sum test should be used

#But as of now, I don't have any groups to compare.
# Detailed comparison using these metric is provided in
# https://rpubs.com/dillmcfarlan/R_microbiotaSOP




