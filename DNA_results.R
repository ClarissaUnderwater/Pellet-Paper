library(tidyverse) 
library(phyloseq)
library(vegan)
library(DESeq2)
library(dendextend) 
library(viridis)
library(dplyr) 
library(forcats)
library(tidyr)
library(colorspace)
library(pheatmap)
library(Hmisc)
library(reshape2)
library(tibble)


######
#setwd("your/directory")

# Input:
# ASVs_counts_16S.tsv (16S count data).
# ASVs_taxonomy_16S.tsv (16S taxonomy data).
# ASVs_counts_18S.tsv (18S count data).
# combined_18S_taxonomy.tsv (processed 18S taxonomy data, generated in the script).
# ASVs_taxonomy_18S.tsv, ASVs_taxonomy_18S_PR2.tsv, and ASVs_taxonomy_18S_MZG.tsv (used for taxonomy combination).
# ST1_BATS_pellets.xlsx

############################################
#### combining 18S taxonomies (done)
############################################ 

ASV_taxonomy_18S <- (read.table("ASVs_taxonomy_18S.tsv", header=T,row.names=1, check.names=F, sep="\t"))
ASV_taxonomy_18S_PR2 <- (read.table("ASVs_taxonomy_18S_PR2.tsv", header=T,row.names=1, check.names=F, sep="\t"))
ASV_taxonomy_18S_MZG <- (read.table("ASVs_taxonomy_18S_MZG.tsv", header=T,row.names=1, check.names=F, sep="\t"))

ASV_taxonomy_18S <- ASV_taxonomy_18S %>% rownames_to_column("ID")
ASV_taxonomy_18S_PR2 <- ASV_taxonomy_18S_PR2 %>% rownames_to_column("ID")
ASV_taxonomy_18S_MZG <- ASV_taxonomy_18S_MZG %>% rownames_to_column("ID")

# Assuming all tables are already loaded with column names as specified
# Join ASV_taxonomy_18S with ASV_taxonomy_18S_PR2 and ASV_taxonomy_18S_MZG by ID
combined_18S_taxonomy <- ASV_taxonomy_18S %>%
  left_join(ASV_taxonomy_18S_PR2, by = "ID", suffix = c("", "_PR2")) %>%
  left_join(ASV_taxonomy_18S_MZG, by = "ID", suffix = c("", "_MZG")) 

# Replace NAs column by column, preferring ASV_taxonomy_18S, then ASV_taxonomy_18S_PR2, then ASV_taxonomy_18S_MZG
combined_18S_taxonomy <- combined_18S_taxonomy %>%
  mutate(
    domain = coalesce(domain, domain_PR2),
    phylum = coalesce(phylum, phylum_PR2),
    class = coalesce(class, class_PR2),
    order = coalesce(order, order_PR2),
    family = coalesce(family, family_PR2),
    genus = coalesce(genus, genus_PR2, Genus),
    species = coalesce(species, species_PR2, Species)
  ) %>%
  select(ID, domain, phylum, class, order, family, genus, species)  # Keep only original columns

# View result
head(combined_18S_taxonomy)

#Set 'ID' as row names and remove the 'ID' column
tax_18S <- combined_18S_taxonomy %>% column_to_rownames("ID")
tax_18S <- as.matrix( tax_18S)

write.table(combined_18S_taxonomy, "combined_18S_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

rm(ASV_taxonomy_18S, ASV_taxonomy_18S_MZG, ASV_taxonomy_18S_PR2)


#####################################################################
# loading data
#####################################################################
Pellets <- read_excel("ST1_BATS_pellets.xlsx", sheet = "all_pellets"
                      , na = c("NA", "na", "", " "))
Pellets <- Pellets[-1, ]  # Remove the second row containing units


# convert numbers to numeric values
Pellets[] <- lapply(Pellets, function(x) {
  if (is.character(x) || is.factor(x)) {
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    if (all(!is.na(x_num) | is.na(x))) {
      return(x_num)
    } else {
      return(x)  # keep original if not safely numeric
    }
  } else {
    return(x)
  }
})

# Convert to numeric (non-numeric entries become NA)
Pellets$sinking_numeric <- as.numeric(Pellets$`sinking velocity`)
Pellets$sinking_numeric[is.infinite(Pellets$sinking_numeric)] <- NA

sample_info_tab <- subset(Pellets, !is.na(Chao16S))


# read 16S
count_16S <- read.table("ASVs_counts_16S.tsv", header=T, row.names=1, check.names=F, sep="\t")
count_16S <- count_16S[, order(names(count_16S))]
tax_16S <- as.matrix(read.table("ASVs_taxonomy_16S.tsv", header=T,row.names=1, check.names=F, sep="\t"))

# read 18S
count_18S <- read.table("ASVs_counts_18S.tsv", header=T, row.names=1, check.names=F, sep="\t")
count_18S <- count_18S[, order(names(count_18S))]
tax_18S <- as.matrix(read.table("combined_18S_taxonomy.tsv", header=T,row.names=2, check.names=F, sep="\t"))
# see combine taxonomy script at end

# match
count_16S <- count_16S[, as.character(sample_info_tab$DNA_name)]
count_18S <- count_18S[, as.character(sample_info_tab$DNA_name)]

sample_info_tab$Taxon <- factor(sample_info_tab$Taxon, levels = c("Salp big", "Salp", "Snail", "Pteropod","Copepod", "Euphausid", "Shrimp"))
sample_info_tab$info_ID <- with(sample_info_tab,fct_reorder(info_ID, as.numeric(Taxon)))

taxon_colors <- c(
  "Salp big" = "#5C351B",
  "Salp" = "#B79600",
  "Copepod" = "#00BFF2",
  "Pteropod" = "#69C380",
  "Snail" = "#F48100",
  "Euphausid" = "purple3",
  "Shrimp" = "#C10000"
)
#####################################################################
# beta diversity, figure 2
#####################################################################

#16S dendrogram
deseq_counts <- DESeqDataSetFromMatrix(count_16S, colData = sample_info_tab, design = ~type) 
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols
labels(euc_dend) <- as.character(sample_info_tab$info_ID[order.dendrogram(euc_dend)])
plot(euc_dend, ylab="VST Euc. dist.")

#18S dendrogram
deseq_counts_18S <- DESeqDataSetFromMatrix(count_18S, colData = sample_info_tab, design = ~Taxon) 
deseq_counts_vst_18S <- varianceStabilizingTransformation(deseq_counts_18S)
vst_trans_count_tab_18S <- assay(deseq_counts_vst_18S)
euc_dist_18S <- dist(t(vst_trans_count_tab_18S))
euc_clust_18S <- hclust(euc_dist_18S, method="ward.D2")
euc_dend_18S <- as.dendrogram(euc_clust_18S, hang=0.1)
dend_cols_18S <- as.character(sample_info_tab$color[order.dendrogram(euc_dend_18S)])
labels_colors(euc_dend_18S) <- dend_cols_18S
labels(euc_dend_18S) <- as.character(sample_info_tab$info_ID[order.dendrogram(euc_dend_18S)])
plot(euc_dend_18S, ylab="VST Euc. dist.")

# ordination plots
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T) # makes ordination
sample_info_tab_phy <- sample_data(sample_info_tab) # make data table

rownames(sample_info_tab_phy) <- sample_names(vst_count_phy)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy) # combine

vst_count_phy_18S <- otu_table(vst_trans_count_tab_18S, taxa_are_rows=T)
vst_physeq_18S <- phyloseq(vst_count_phy_18S, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
vst_pcoa_18S <- ordinate(vst_physeq_18S, method="MDS", distance="euclidean")
eigen_vals_18S <- vst_pcoa_18S$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

#pdf("ordination_plot_16S.pdf", width=8, height=6)
ordination_plot <-
plot_ordination(vst_physeq, vst_pcoa, col="Taxon") + 
  geom_point(size=1) + 
  labs(col="Taxon") + 
  geom_text(aes(label=sample_info_tab$info_ID, hjust=0.3, vjust=-0.4, size = 1)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values = taxon_colors) +
  theme_bw() + theme(legend.position="none")
print(ordination_plot)
#dev.off() 

#pdf("ordination_plot_18S.pdf", width=8, height=6)
ordination_plot <-
  plot_ordination(vst_physeq_18S, vst_pcoa_18S, color="Taxon") + 
  geom_point(size=1) + labs(col="Taxon") + 
  geom_text(aes(label=sample_info_tab$info_ID, hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals_18S[2]/eigen_vals_18S[1])) + ggtitle("PCoA") + 
  scale_color_manual(values = taxon_colors) +
  theme_bw() + theme(legend.position="none")
print(ordination_plot)
#dev.off()

#pdf("ordination_plot_18S_crop.pdf", width=8, height=6)
ordination_plot <-
  plot_ordination(vst_physeq_18S, vst_pcoa_18S, color="Taxon") + 
  geom_point(size=1) + labs(col="Taxon") + 
  geom_text(aes(label=sample_info_tab$info_ID, hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals_18S[2]/eigen_vals_18S[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=taxon_colors) + 
  ylim(-21,0)+
  xlim(-25,-5)+
  theme_bw() + theme(legend.position="none")
print(ordination_plot)
#dev.off()

#####################################################################
# alpha diversity
#####################################################################
# plot rarefaction (slow)
#rarecurve(t(count_16S), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs", label=F)
#abline(v=(min(rowSums(t(count_16S))))) # and adding a vertical line at the fewest seqs in any sample
#rarecurve(t(count_18S), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs", label=F)
  
############ create a phyloseq object using our un-transformed count table
 count_tab_phy <- otu_table(count_16S, taxa_are_rows=T)
 tax_tab_phy <- tax_table(tax_16S)
 ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
 
 count_tab_phy_18S <- otu_table(count_18S, taxa_are_rows=T)
 tax_tab_phy_18S <- tax_table(tax_18S)
 ASV_physeq_18S <- phyloseq(count_tab_phy_18S, tax_tab_phy_18S, sample_info_tab_phy)
 
 ############ Chao1 is richness (number of ASVs) and Shannon is diversity (evenness of richness)
 
 # call the plot_richness() function on our phyloseq object
 plot_richness(ASV_physeq, x="Station", color="Taxon", measures=c("Chao1", "Shannon")) + 
   theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
 
 # plot against respiration
 plot_richness(ASV_physeq, x="resp", color="Taxon", measures=c("Chao1", "Shannon")) + 
    theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


plot_richness(ASV_physeq_18S, shape="Station", color="Taxon", measures=c("Chao1", "Shannon")) + 
   scale_shape_manual(values=c(20,17,4,13,5)) +
   #scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$color)])) +
   theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

 plot_richness(ASV_physeq, x="Taxon",  color="Station", measures=c("Chao1", "Shannon")) + 
   scale_shape_manual(values=c(20,17,4,13,5)) +
   #scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$color)])) +
   theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
 
 
 # Calculate the alpha diversity metrics for ASV_physeq 
 alpha_diversity_16S <- estimate_richness(ASV_physeq, measures = c("Chao1", "Shannon"))
 alpha_diversity_16S$DNA_name <- sub("^X", "", rownames(alpha_diversity_16S))
 colnames(alpha_diversity_16S) <- c("Chao1_16S", "SE_Chao1_16S","Shannon_16S","DNA_name")
 alpha_diversity_18S <- estimate_richness(ASV_physeq_18S, measures = c("Chao1", "Shannon"))
 alpha_diversity_18S$DNA_name <- sub("^X", "", rownames(alpha_diversity_18S))
 colnames(alpha_diversity_18S) <- c("Chao1_18S", "SE_Chao1_18S","Shannon_18S","DNA_name")
 
 # export
# write.csv(alpha_diversity_16S, file = "alpha_diversity_16S.csv", row.names = TRUE)
 # write.csv(alpha_diversity_18S, file = "alpha_diversity_18S.csv", row.names = TRUE)
 
#####################################################################
# taxonomy with phyloseq
#####################################################################
 # using phyloseq to make a count table that has summed all ASVs that were in the same phylum
 phyla_counts_16S <- otu_table(tax_glom(ASV_physeq, taxrank="phylum")) 
 phyla_counts_18S <- otu_table(tax_glom(ASV_physeq_18S, taxrank="class")) 
 
 # making a vector of phyla names and set as row names
 phyla_tax_vec_16S <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="phylum"))[,"phylum"]) 
 phyla_tax_vec_18S <- as.vector(tax_table(tax_glom(ASV_physeq_18S, taxrank="class"))[,"class"]) 
 rownames(phyla_counts_16S) <- as.vector(phyla_tax_vec_16S)
 rownames(phyla_counts_18S) <- as.vector(phyla_tax_vec_18S)
 
 # we also have to account for NA sequences that weren't assigned any taxonomy even at the phylum level 
 unclassified_tax_counts_16S <- colSums(count_16S) - colSums(phyla_counts_16S)
 unclassified_tax_counts_18S <- colSums(count_18S) - colSums(phyla_counts_18S)
 
 # and we'll add this row to our phylum count table:
 phyla_and_unidentified_counts_16S <- rbind(phyla_counts_16S, "Unclassified"=unclassified_tax_counts_16S)
 phyla_and_unidentified_counts_18S <- rbind(phyla_counts_18S, "Unclassified"=unclassified_tax_counts_18S)

 
 #####################################################################
 # details proteobacteria for 16S
 #####################################################################
 
 # remove the Proteobacteria, so we can next add them back in broken down by class
 temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_16S[!row.names(phyla_and_unidentified_counts_16S) %in% "Proteobacteria", ]
 
 # making count table broken down by class (contains classes beyond the Proteobacteria too at this point)
 class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="class")) 
 
 # making a table that holds the phylum and class level info
 class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="class")) 
 
 phy_tmp_vec <- class_tax_phy_tab[,2]
 class_tmp_vec <- class_tax_phy_tab[,3]
 rows_tmp <- row.names(class_tax_phy_tab)
 class_tax_tab <- data.frame("phylum"=phy_tmp_vec, "class"=class_tmp_vec, row.names = rows_tmp)
 
 # making a vector of just the Proteobacteria classes
 proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$phylum == "Proteobacteria", "class"])
 
 # changing the row names like above so that they correspond to the taxonomy, rather than an ASV identifier
 rownames(class_counts_tab) <- as.vector(class_tax_tab$class) 
 
 # making a table of the counts of the Proteobacterial classes
 proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 
 
 # there are also possibly some some sequences that were resolved to the level of Proteobacteria, but not any further, and therefore would be missing from our class table
 proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_16S[row.names(phyla_and_unidentified_counts_16S) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)
 
 # now combining the tables:
 major_taxa_counts_16S <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)
 major_taxa_counts_18S <- phyla_and_unidentified_counts_18S
 
 # check we didn't miss any other sequences: compare the column sums to see if they are the same
 identical(colSums(major_taxa_counts_16S), colSums(count_16S)) 
 
 #####################################################################
 # filter tables
 #####################################################################

# 16S
# generate a proportions table for summarizing:
 major_taxa_proportions_16S <- apply(major_taxa_counts_16S, 2, function(x) x/sum(x)*100)
 dim(major_taxa_proportions_16S)
 # only keep rows (taxa) that make up greater than 0.1% in any individual sample
 temp_filt_major_taxa_proportions_16S <- data.frame(major_taxa_proportions_16S[apply(major_taxa_proportions_16S, 1, max) > 0.1, ])
 dim(temp_filt_major_taxa_proportions_16S) 
 # add a row called "Other" that keeps track of how much we filtered out (which will also keep our totals at 100%)
 filtered_proportions_16S <- colSums(major_taxa_proportions_16S) - colSums(temp_filt_major_taxa_proportions_16S)
 filt_major_taxa_proportions_16S <- rbind(temp_filt_major_taxa_proportions_16S, "Other"=filtered_proportions_16S)
 # copy table
 filt_major_taxa_proportions_16S_for_plot <- filt_major_taxa_proportions_16S
 # and add a column of the taxa names so that it is within the table, rather than just as row names (this makes working with ggplot easier)
 filt_major_taxa_proportions_16S_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_16S_for_plot)
 # now we'll transform the table into narrow, or long, format (also makes plotting easier)
 filt_major_taxa_proportions_16S_for_plot.g <- pivot_longer(filt_major_taxa_proportions_16S_for_plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame()
 # take a look at the new table and compare it with the old one
 head(filt_major_taxa_proportions_16S_for_plot.g)
 head(filt_major_taxa_proportions_16S_for_plot)

 
 # 18S
 # generate a proportions table for summarizing:
 major_taxa_proportions_18S <- apply(major_taxa_counts_18S, 2, function(x) x/sum(x)*100)
 dim(major_taxa_proportions_18S)
 # only keep rows (taxa) that make up greater than 5% in any individual sample
 temp_filt_major_taxa_proportions_18S <- data.frame(major_taxa_proportions_18S[apply(major_taxa_proportions_18S, 1, max) > 5, ])
 dim(temp_filt_major_taxa_proportions_18S) 
 # add a row called "Other" that keeps track of how much we filtered out (which will also keep our totals at 100%)
 filtered_proportions_18S <- colSums(major_taxa_proportions_18S) - colSums(temp_filt_major_taxa_proportions_18S)
 filt_major_taxa_proportions_18S <- rbind(temp_filt_major_taxa_proportions_18S, "Other"=filtered_proportions_18S)
 # copy table
 filt_major_taxa_proportions_18S_for_plot <- filt_major_taxa_proportions_18S
 # and add a column of the taxa names so that it is within the table, rather than just as row names (this makes working with ggplot easier)
 filt_major_taxa_proportions_18S_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_18S_for_plot)
 # now we'll transform the table into narrow, or long, format (also makes plotting easier)
 filt_major_taxa_proportions_18S_for_plot.g <- pivot_longer(filt_major_taxa_proportions_18S_for_plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame()
 # take a look at the new table and compare it with the old one
 head(filt_major_taxa_proportions_18S_for_plot.g)
 head(filt_major_taxa_proportions_18S_for_plot)

 #####################################################################
 #### preparing and merging metadata
 #####################################################################
 # table with "color" and "characteristics" of each sample 
 sample_info_for_merge<-data.frame("Taxon"=sample_info_tab$Taxon, 
                                   "Station"=sample_info_tab$Station,
                                   "Organism"=sample_info_tab$detailed,
                                   "shape"=sample_info_tab$shape,
                                   "resp"=sample_info_tab$Respiration,
                                   "sinking_v"=sample_info_tab$`sinking velocity`,
                                   "Perimeter_mm"=sample_info_tab$perimeter,
                                   "degr_less"=sample_info_tab$degradation,
                                   "Area_mm"=sample_info_tab$`2D area`,
                                   "info_ID"=sample_info_tab$info_ID,
                                   "animal_ID"=sample_info_tab$animal_ID,
                                   stringsAsFactors=F)
 head(sample_info_for_merge)
 
  # Reorder animal_ID by Taxon
# sample_info_for_merge$info_ID <- reorder(sample_info_for_merge$info_ID,sample_info_for_merge$Taxon)
 sample_info_for_merge$info_ID <- with(sample_info_for_merge, fct_reorder(info_ID, as.numeric(Taxon))
 )
 
 #####################################################################
 #### proportions table
 #####################################################################
 
 # 16S
 # removing X
 filt_major_taxa_proportions_16S_for_plot.g$Sample <- gsub("X","",filt_major_taxa_proportions_16S_for_plot.g$Sample)
 head(filt_major_taxa_proportions_16S_for_plot.g)
 # merging this table with the plotting table we just made
 filt_major_taxa_proportions_16S_for_plot.g2 <- merge(filt_major_taxa_proportions_16S_for_plot.g, sample_info_for_merge)
 head(filt_major_taxa_proportions_16S_for_plot.g2)
 
 
 # 18S
 # removing X
 filt_major_taxa_proportions_18S_for_plot.g$Sample <- gsub("X","",filt_major_taxa_proportions_18S_for_plot.g$Sample)
 head(filt_major_taxa_proportions_18S_for_plot.g)
 # merging this table with the plotting table we just made
 filt_major_taxa_proportions_18S_for_plot.g2 <- merge(filt_major_taxa_proportions_18S_for_plot.g, sample_info_for_merge)
 head(filt_major_taxa_proportions_18S_for_plot.g2)

 
 #####################################################################
 #### stacked bar charts, figure 1
 #####################################################################
 
 ggplot(filt_major_taxa_proportions_16S_for_plot.g2, aes(x=info_ID, y=Proportion, fill=Major_Taxa)) +
   geom_bar(width=0.6, stat="identity", color = "black") +
   theme_bw() +
   # facet_wrap(~hull)+
   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
   labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")
 ggsave("16S_barplot_tax.pdf", width=9, height=12, dpi=300)
 
 ggplot(filt_major_taxa_proportions_18S_for_plot.g2, aes(x=info_ID, y=Proportion, fill=Major_Taxa)) +
   geom_bar(width=0.6, stat="identity", color="black") +
   #scale_fill_manual(values = palette_colors)+
   theme_bw() +
   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
   labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")
 ggsave("bar_18S_5pct.pdf", width=9, height=12, dpi=300)
 
 
 #####################################################################
 #### ASV correlations, figure S2
 #####################################################################
 
 #major_class_proportions_16S contains rownames which are microorganisms, 
 # column names which are samples and percentages in the table
 #major_class_proportions_18S contains rownames which are macroorganisms, 
 # column names which are the same samples and percentages in the table
 # how to calculate a correlation matrix testing the co-occurrence of all possible pairs of micro- and macroorganisms?
 
 
 # Create new data frames and rename the variables
 #df_16S <- as.data.frame(major_taxa_proportions_16S)
 #df_16S <- as.data.frame(filt_major_taxa_proportions_16S)
 #df_16S <- data.frame(df_16S[apply(df_16S, 1, max) > 5, ])
 df_16S <- df_16S[rownames(df_16S) != "Unclassified", ]
 
 #df_18S <- as.data.frame(major_taxa_proportions_18S)
 df_18S <- as.data.frame(filt_major_taxa_proportions_18S)
 df_18S <- data.frame(df_18S[apply(df_18S, 1, max) > 5, ])
 df_18S <- df_18S[rownames(df_18S) != "Unclassified", ]
 
 df_16S_t <- t(df_16S)
 df_18S_t <- t(df_18S)
 
 correlation_matrix <- cor(df_18S_t, df_16S_t, use = "pairwise.complete.obs")
 
 #cor_test_result <- cor.test(df_18S_t, df_16S_t, use = "pairwise.complete.obs")
 #cor_test_result$p.value
 
 p_matrix <- matrix(NA, nrow = ncol(df_18S_t), ncol = ncol(df_16S_t),
                    dimnames = list(colnames(df_18S_t), colnames(df_16S_t)))
 
 r_matrix <- p_matrix  # To store correlation coefficients
 
 for (i in seq_len(ncol(df_18S_t))) {
   for (j in seq_len(ncol(df_16S_t))) {
     x <- df_18S_t[, i]
     y <- df_16S_t[, j]
     # Check for valid input lengths
     if (sum(!is.na(x) & !is.na(y)) >= 3) {  # Minimum 3 points for cor.test
       test <- cor.test(x, y, use = "pairwise.complete.obs")
       p_matrix[i, j] <- test$p.value
       r_matrix[i, j] <- test$estimate
     }
   }
 }
 
 
 rownames(correlation_matrix) <- as.character(rownames(correlation_matrix))
 colnames(correlation_matrix) <- as.character(colnames(correlation_matrix))
 
 rownames(correlation_matrix) <- trimws(rownames(correlation_matrix))
 colnames(correlation_matrix) <- trimws(colnames(correlation_matrix))
 
 sorted_correlation_matrix <- correlation_matrix[order(rownames(correlation_matrix)), 
                                                 order(colnames(correlation_matrix))]
 
 
 # Create and save the heatmap with correlation values displayed
# pdf("heatmap_with_values.pdf", width = 12, height = 6)
 heatmap <-
   pheatmap(
     sorted_correlation_matrix,
     cluster_rows = FALSE, # Hierarchical clustering for rows
     cluster_cols = FALSE, # Hierarchical clustering for columns
     main = "Heatmap of Microorganism-Macroorganism Correlations",
     display_numbers = TRUE, # Display correlation coefficients on heatmap
     number_format = "%.2f", # Format numbers to two decimal places
     fontsize_number = 10,   # Font size for the numbers
     color = colorRampPalette(c("cornflowerblue","white", "brown3"))(50) # Color gradient
   )
 print(heatmap)
# dev.off()
 
 sorted_p_matrix <- p_matrix[order(rownames(p_matrix)), 
                             order(colnames(p_matrix))]
 
 pheatmap(sorted_p_matrix,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          display_numbers = TRUE,
          fontsize = 10,
          main = "Heatmap of p-values")
 
 #####################################################################
## statistical tests for beta diversity
 ##################################################################### 
 #####################################################################
 # 16S
 #############################################################################################################################
 ############ do they group by station?
 #16S dendrogram
 deseq_counts_16S <- DESeqDataSetFromMatrix(count_16S, colData = sample_info_tab, design = ~Station) 
 deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
 vst_trans_count_tab <- assay(deseq_counts_vst)
 euc_dist <- dist(t(vst_trans_count_tab))
 euc_clust <- hclust(euc_dist, method="ward.D2")
 euc_dend <- as.dendrogram(euc_clust, hang=0.1)
 dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
 labels_colors(euc_dend) <- dend_cols
 labels(euc_dend) <- as.character(sample_info_tab$info_ID[order.dendrogram(euc_dend)])
 plot(euc_dend, ylab="VST Euc. dist.")
 
 rel_abund_16S <- apply(count_16S, 2, function(x) x / sum(x))
 bray_dist_16S <- vegdist(t(rel_abund_16S), method = "bray") # Compute Bray–Curtis dissimilarity
 bray_no_salp <- vegdist(t(rel_abund_16S[, df_no_salp$DNA_name]), method = "bray")
 
 adonis2(bray_dist_16S ~ Station, data = sample_info_tab) #Test for compositional differences (PERMANOVA) by station 
 betadisper_16S <- betadisper(bray_dist_16S, sample_info_tab$Station) # Test for dispersion differences (PERMDISP) by station
 anova(betadisper_16S)
 
 # no salp 
 betadisper_no_salp <- betadisper(bray_no_salp, df_no_salp$Station)
 anova(betadisper_no_salp)
 adonis2(bray_no_salp ~ Station, data = df_no_salp) #Test for compositional differences (PERMANOVA) by station 
 
 
 ############ do they group by taxon?
 #16S dendrogram
 deseq_counts <- DESeqDataSetFromMatrix(count_16S, colData = sample_info_tab, design = ~Taxon) 
 deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
 vst_trans_count_tab <- assay(deseq_counts_vst)
 euc_dist <- dist(t(vst_trans_count_tab))
 euc_clust <- hclust(euc_dist, method="ward.D2")
 euc_dend <- as.dendrogram(euc_clust, hang=0.1)
 dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
 labels_colors(euc_dend) <- dend_cols
 labels(euc_dend) <- as.character(sample_info_tab$info_ID[order.dendrogram(euc_dend)])
 plot(euc_dend, ylab="VST Euc. dist.")
 
 rel_abund_16S <- apply(count_16S, 2, function(x) x / sum(x))
 bray_dist_16S <- vegdist(t(rel_abund_16S), method = "bray") # Compute Bray–Curtis dissimilarity
 bray_no_salp <- vegdist(t(rel_abund_16S[, df_no_salp$DNA_name]), method = "bray")
 
 adonis2(bray_dist_16S ~ Taxon, data = sample_info_tab) #Test for compositional differences (PERMANOVA) by station 
 betadisper_16S <- betadisper(bray_dist_16S, sample_info_tab$Taxon) # Test for dispersion differences (PERMDISP) by station
 anova(betadisper_16S)
 
 # no salp 
 bray_no_salp <- vegdist(t(rel_abund_16S[, df_no_salp$DNA_name]), method = "bray")
 betadisper_no_salp <- betadisper(bray_no_salp, df_no_salp$Taxon)
 anova(betadisper_no_salp)
 adonis2(bray_no_salp ~ Taxon, data = df_no_salp) #Test for compositional differences (PERMANOVA) by station 
 
 ############ do they group by organism?
 #16S dendrogram
 deseq_counts_16S <- DESeqDataSetFromMatrix(count_16S, colData = sample_info_tab, design = ~Organism) 
 deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
 vst_trans_count_tab <- assay(deseq_counts_vst)
 euc_dist <- dist(t(vst_trans_count_tab))
 euc_clust <- hclust(euc_dist, method="ward.D2")
 euc_dend <- as.dendrogram(euc_clust, hang=0.1)
 dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
 labels_colors(euc_dend) <- dend_cols
 labels(euc_dend) <- as.character(sample_info_tab$info_ID[order.dendrogram(euc_dend)])
 plot(euc_dend, ylab="VST Euc. dist.")
 
 rel_abund_16S <- apply(count_16S, 2, function(x) x / sum(x))
 bray_dist_16S <- vegdist(t(rel_abund_16S), method = "bray") # Compute Bray–Curtis dissimilarity
 bray_no_salp <- vegdist(t(rel_abund_16S[, df_no_salp$DNA_name]), method = "bray")
 
 # include only organisms with at least 2 pellets and no salp
 table(df_no_salp$Organism)
 organisms_with_reps <- names(which(table(df_no_salp$Organism) >= 2))
 df_replicated <- subset(df_no_salp, Organism %in% organisms_with_reps)
 shared_ids <- intersect(colnames(rel_abund_16S), df_replicated$DNA_name)
 rel_abund_16S_repl <- rel_abund_16S[, shared_ids]
 df_replicated <- df_replicated[df_replicated$DNA_name %in% shared_ids, ]
 bray_repl <- vegdist(t(rel_abund_16S_repl), method = "bray")
 
 adonis2(bray_repl ~ Organism, data = df_replicated)
 betadisper_repl <- betadisper(bray_repl, df_replicated$Organism)
 anova(betadisper_repl)
 
 ############ do they group by organism only euphausiid?
 # include only organisms with at least 2 pellets
 table(df_no_salp$Organism)
 df_replicated <- subset(df_no_salp, tx %in% "EP")
 shared_ids <- intersect(colnames(rel_abund_16S), df_replicated$DNA_name)
 rel_abund_16S_repl <- rel_abund_16S[, shared_ids]
 df_replicated <- df_replicated[df_replicated$DNA_name %in% shared_ids, ]
 bray_repl <- vegdist(t(rel_abund_16S_repl), method = "bray")
 
 adonis2(bray_repl ~ Organism, data = df_replicated)
 betadisper_repl <- betadisper(bray_repl, df_replicated$Organism)
 anova(betadisper_repl)
 
 
 ### do they group by station or taxon or organism?
 rel_abund_16S <- apply(count_16S, 2, function(x) x / sum(x))
 bray_dist_16S <- vegdist(t(rel_abund_16S), method = "bray") # Compute Bray–Curtis dissimilarity
 
 adonis2(bray_dist_16S ~ Taxon, data = sample_info_tab) #Test for compositional differences (PERMANOVA) by station 
 betadisper_16S <- betadisper(bray_dist_16S, sample_info_tab$Taxon) # Test for dispersion differences (PERMDISP) by station
 
 #adonis2(bray_dist_16S ~ Station, data = sample_info_tab) #Test for compositional differences (PERMANOVA) by station 
 #betadisper_16S <- betadisper(bray_dist_16S, sample_info_tab$Station) # Test for dispersion differences (PERMDISP) by station
 anova(betadisper_16S)
 
 # is it driven by the the salp?
 df_no_salp <- subset(sample_info_tab, tx != "SP") 
 bray_no_salp <- vegdist(t(rel_abund_16S[, df_no_salp$DNA_name]), method = "bray")
 betadisper_no_salp <- betadisper(bray_no_salp, df_no_salp$Taxon)
 anova(betadisper_no_salp)
 
 # is it driven by the the small salp?
 df_no_single <- subset(sample_info_tab, Taxon != "Salp") 
 bray_no_salp <- vegdist(t(rel_abund_16S[, df_no_single$DNA_name]), method = "bray")
 betadisper_no_salp <- betadisper(bray_no_salp, df_no_single$Taxon)
 anova(betadisper_no_salp)
 
 ## which  microbes are different by taxa?
 count_16S_no_salp <- count_16S[, df_no_salp$DNA_name]
 dds <- DESeqDataSetFromMatrix(countData = count_16S_no_salp,
                               colData = df_no_salp,
                               design = ~ Taxon)
 dds <- DESeq(dds)
 res <- results(dds)
 res_sig <- res[which(res$padj < 0.05), ]
 
 tax_df <- as.data.frame(tax_16S) |> rownames_to_column("ASV")
 groups <- levels(colData(dds)$Taxon)
 pairs  <- t(combn(groups, 2))             # rows are [B, A]
 
 one_contrast <- function(B, A){
   results(dds, contrast = c("Taxon", B, A)) |>
     as.data.frame() |>
     rownames_to_column("ASV") |>
     mutate(contrast  = paste(B, "vs", A),
            direction = case_when(
              log2FoldChange >  0 ~ paste("Higher in", B),
              log2FoldChange <  0 ~ paste("Higher in", A),
              TRUE                ~ "No change"
            ))
 }
 
 res_all <- do.call(rbind, apply(pairs, 1, \(x) one_contrast(B = x[1], A = x[2])))
 
 res_sig <- res_all |>
   left_join(tax_df, by = "ASV") |>
   filter(!is.na(padj), padj < 0.05) |>
   arrange(contrast, padj, desc(abs(log2FoldChange)))
 
 # save
 #write.table(res_sig, "DESeq2_16S_byTaxon_sig_ASVs.tsv", sep="\t", quote=FALSE, row.names=FALSE)
 
 
 # df_no_salp$Taxon or $Organism identifies the animal
 dominant_taxa <- rel_abund_16S_no_salp %>%
   as.data.frame() %>%
   rownames_to_column("ASV") %>%
   pivot_longer(-ASV, names_to = "Sample", values_to = "RelAbund") %>%
   left_join(df_no_salp[, c("DNA_name","Taxon")],
             by = c("Sample"="DNA_name")) %>%
   group_by(Taxon, ASV) %>%
   summarise(mean_abundance = mean(RelAbund), .groups="drop") %>%
   arrange(Taxon, desc(mean_abundance)) %>%
   group_by(Taxon) %>%
   slice_max(mean_abundance, n = 10)   # top 10 per animal type
 
 # Subset the raw counts to no-salp samples
 count_16S_no_salp <- count_16S[, df_no_salp$DNA_name]
 count_18S_no_salp <- count_18S[, df_no_salp$DNA_name]
 
 # Compute relative abundances (each column sums to 1)
 rel_abund_16S_no_salp <- apply(count_16S_no_salp, 2, function(x) x / sum(x))
 rel_abund_18S_no_salp <- apply(count_18S_no_salp, 2, function(x) x / sum(x))
 
 # Check the result
 dim(rel_abund_16S_no_salp)
 colSums(rel_abund_16S_no_salp)[1:5]   # should all be 1
 
 tax_df <- as.data.frame(tax_16S) |> rownames_to_column("ASV")
 dominant_taxa <- left_join(dominant_taxa, tax_df, by = "ASV")
 
 library(indicspecies)
 comm <- t(rel_abund_16S)
 indval <- multipatt(comm, sample_info_tab$Taxon, func="r.g", duleg=TRUE, control=how(nperm=999))
 summary(indval)
 
 
 comm <- t(rel_abund_18S)
 indval <- multipatt(comm, sample_info_tab$Taxon, func="r.g", duleg=TRUE, control=how(nperm=999))
 summary(indval)
 #############################################################################################################################
 # test for links between 16S and 18S
 #############################################################################################################################
 
 # raw or relative abundance matrices
 # rows = ASVs, columns = samples (now including salps)
 shared_samples <- intersect(colnames(rel_abund_16S), colnames(rel_abund_18S))
 
 m16 <- rel_abund_16S[, shared_samples]
 m18 <- rel_abund_18S[, shared_samples]
 
 # optional: filter rare ASVs to reduce noise
 keep16 <- rowSums(m16 > 0) >= 3     # present in ≥3 samples
 keep18 <- rowSums(m18 > 0) >= 3
 m16 <- m16[keep16, ]
 m18 <- m18[keep18, ]
 
 
 # 1) Align samples and drop zero-variance features
 shared <- intersect(colnames(rel_abund_16S), colnames(rel_abund_18S))
 m16 <- rel_abund_16S[, shared, drop=FALSE]
 m18 <- rel_abund_18S[, shared, drop=FALSE]
 
 # keep features present in >= 3 samples and with variance > 0
 keep16 <- rowSums(m16 > 0) >= 3 & apply(m16, 1, var) > 0
 keep18 <- rowSums(m18 > 0) >= 3 & apply(m18, 1, var) > 0
 m16 <- m16[keep16, , drop=FALSE]
 m18 <- m18[keep18, , drop=FALSE]
 
 # 2) Build a single matrix with FEATURES as COLUMNS (required by rcorr)
 m_all <- rbind(m16, m18)     # features = rows, samples = cols
 X     <- t(m_all)            # now: samples = rows, features = cols
 
 n16 <- nrow(m16)             # number of 16S features
 n18 <- nrow(m18)             # number of 18S features
 
 # 3) Correlate across features (columns)
 cor_matrix <- rcorr(as.matrix(X), type = "spearman")
 R <- cor_matrix$r
 P <- cor_matrix$P
 
 # 4) Extract ONLY 16S–18S block (first n16 are 16S, next n18 are 18S)
 # (because we stacked m16 above m18 in rbind)
 idx16 <- 1:n16
 idx18 <- (n16 + 1):(n16 + n18)
 
 R_16S_18S <- R[idx16, idx18, drop=FALSE]
 P_16S_18S <- P[idx16, idx18, drop=FALSE]
 
 # 5) Turn into a tidy table and filter strong links
 cor_df <- melt(R_16S_18S, varnames = c("ASV_16S", "ASV_18S"), value.name = "rho") %>%
   mutate(ASV_16S = rownames(m16)[ASV_16S],
          ASV_18S = rownames(m18)[ASV_18S])
 
 p_df <- melt(P_16S_18S, varnames = c("ASV_16S", "ASV_18S"), value.name = "p")
 
 links <- left_join(cor_df, p_df, by = c("ASV_16S","ASV_18S")) %>%
   filter(!is.na(rho), !is.na(p),
          abs(rho) >= 0.6, p < 0.05) %>%
   arrange(desc(abs(rho)))
 
 # 6) (Optional) annotate taxonomy
 tax16_df <- as.data.frame(tax_16S) %>% tibble::rownames_to_column("ASV_16S")
 tax18_df <- as.data.frame(tax_18S) %>% tibble::rownames_to_column("ASV_18S")
 links_annot <- links %>%
   left_join(tax16_df, by="ASV_16S") %>%
   left_join(tax18_df, by="ASV_18S", suffix=c("_16S","_18S"))
 
 
 # add FDR to your links table
 links_annot$padj <- p.adjust(links_annot$p, method = "BH")
 
 # keep only robust links (edit thresholds as you like)
 links_filt <- subset(links_annot, !is.na(rho) & !is.na(padj) &
                        abs(rho) >= 0.6 & padj < 0.05)
 
 summary_phylum_class <- links_filt |>
   dplyr::group_by(phylum_16S, class_18S) |>
   dplyr::summarise(
     n_links       = dplyr::n(),
     mean_rho      = mean(rho, na.rm = TRUE),
     strongest_rho = max(abs(rho), na.rm = TRUE),
     direction     = ifelse(mean(rho, na.rm = TRUE) > 0, "positive", "negative"),
     .groups = "drop"
   ) |>
   dplyr::arrange(dplyr::desc(n_links), dplyr::desc(abs(mean_rho)))
 
 summary_genus_order <- links_filt |>
   dplyr::group_by(genus_16S, order_18S) |>
   dplyr::summarise(
     n_links       = dplyr::n(),
     mean_rho      = mean(rho, na.rm = TRUE),
     strongest_rho = max(abs(rho), na.rm = TRUE),
     direction     = ifelse(mean(rho, na.rm = TRUE) > 0, "positive", "negative"),
     .groups = "drop"
   ) |>
   dplyr::arrange(dplyr::desc(n_links), dplyr::desc(abs(mean_rho)))
 
 top_pairs <- links_filt |>
   dplyr::select(ASV_16S, genus_16S, family_16S, phylum_16S,
                 ASV_18S, genus_18S, family_18S, class_18S, order_18S,
                 rho, padj) |>
   dplyr::arrange(dplyr::desc(abs(rho)), padj)
 
 ggplot(summary_phylum_class,
        aes(x = phylum_16S, y = class_18S, size = n_links, fill = mean_rho)) +
   geom_point(shape = 21, color = "black", alpha = 0.9) +
   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   labs(x = "Bacterial phylum (16S)", y = "Eukaryotic class (18S)",
        size = "# ASV links", fill = "Mean ρ")
 
 write.table(summary_phylum_class, file = "summary_phylum_class.csv", sep = ",", row.names = FALSE, col.names = TRUE)
 
 unique(summary_phylum_class$class_18S)
 unique(summary_phylum_class$phylum_16S )
 #############################################################################################################################
 # 18S
 ############################################################
 #18S dendrogram
 
 ############ do they group by station?
 deseq_counts_18S <- DESeqDataSetFromMatrix(count_18S, colData = sample_info_tab, design = ~Station) 
 deseq_counts_vst_18S <- varianceStabilizingTransformation(deseq_counts_18S)
 vst_trans_count_tab_18S <- assay(deseq_counts_vst_18S)
 euc_dist_18S <- dist(t(vst_trans_count_tab_18S))
 euc_clust_18S <- hclust(euc_dist_18S, method="ward.D2")
 euc_dend_18S <- as.dendrogram(euc_clust_18S, hang=0.1)
 dend_cols_18S <- as.character(sample_info_tab$color[order.dendrogram(euc_dend_18S)])
 labels_colors(euc_dend_18S) <- dend_cols_18S
 labels(euc_dend_18S) <- as.character(sample_info_tab$info_ID[order.dendrogram(euc_dend_18S)])
 
 plot(euc_dend_18S, ylab="VST Euc. dist.")
 rel_abund_18S <- apply(count_18S, 2, function(x) x / sum(x))
 bray_dist_18S <- vegdist(t(rel_abund_18S), method = "bray") # Compute Bray–Curtis dissimilarity
 bray_no_salp <- vegdist(t(rel_abund_18S[, df_no_salp$DNA_name]), method = "bray")
 
 adonis2(bray_dist_18S ~ Station, data = sample_info_tab) #Test for compositional differences (PERMANOVA) by station 
 betadisper_18S <- betadisper(bray_dist_18S, sample_info_tab$Station) # Test for dispersion differences (PERMDISP) by station
 anova(betadisper_18S)
 
 # no salp 
 betadisper_no_salp <- betadisper(bray_no_salp, df_no_salp$Station)
 anova(betadisper_no_salp)
 adonis2(bray_no_salp ~ Station, data = df_no_salp) #Test for compositional differences (PERMANOVA) by station 
 
 
 ############ do they group by taxon?
 deseq_counts_18S <- DESeqDataSetFromMatrix(count_18S, colData = sample_info_tab, design = ~Taxon) # test either by station or type, rerun accordingly until anova 
 deseq_counts_vst_18S <- varianceStabilizingTransformation(deseq_counts_18S)
 vst_trans_count_tab_18S <- assay(deseq_counts_vst_18S)
 euc_dist_18S <- dist(t(vst_trans_count_tab_18S))
 euc_clust_18S <- hclust(euc_dist_18S, method="ward.D2")
 euc_dend_18S <- as.dendrogram(euc_clust_18S, hang=0.1)
 dend_cols_18S <- as.character(sample_info_tab$color[order.dendrogram(euc_dend_18S)])
 labels_colors(euc_dend_18S) <- dend_cols_18S
 labels(euc_dend_18S) <- as.character(sample_info_tab$info_ID[order.dendrogram(euc_dend_18S)])
 
 plot(euc_dend_18S, ylab="VST Euc. dist.")
 rel_abund_18S <- apply(count_18S, 2, function(x) x / sum(x))
 bray_dist_18S <- vegdist(t(rel_abund_18S), method = "bray") # Compute Bray–Curtis dissimilarity
 bray_no_salp <- vegdist(t(rel_abund_18S[, df_no_salp$DNA_name]), method = "bray")
 
 adonis2(bray_dist_18S ~ Taxon, data = sample_info_tab) #Test for compositional differences (PERMANOVA) by station 
 betadisper_18S <- betadisper(bray_dist_18S, sample_info_tab$Taxon) # Test for dispersion differences (PERMDISP) by station
 anova(betadisper_18S)
 
 # no salp 
 bray_no_salp <- vegdist(t(rel_abund_18S[, df_no_salp$DNA_name]), method = "bray")
 betadisper_no_salp <- betadisper(bray_no_salp, df_no_salp$Taxon)
 anova(betadisper_no_salp)
 adonis2(bray_no_salp ~ Taxon, data = df_no_salp) #Test for compositional differences (PERMANOVA) by station 
 
 
 ############ do they group by organism?
 deseq_counts_18S <- DESeqDataSetFromMatrix(count_18S, colData = sample_info_tab, design = ~Organism) # test either by station or type, rerun accordingly until anova 
 deseq_counts_vst_18S <- varianceStabilizingTransformation(deseq_counts_18S)
 vst_trans_count_tab_18S <- assay(deseq_counts_vst_18S)
 euc_dist_18S <- dist(t(vst_trans_count_tab_18S))
 euc_clust_18S <- hclust(euc_dist_18S, method="ward.D2")
 euc_dend_18S <- as.dendrogram(euc_clust_18S, hang=0.1)
 dend_cols_18S <- as.character(sample_info_tab$color[order.dendrogram(euc_dend_18S)])
 labels_colors(euc_dend_18S) <- dend_cols_18S
 labels(euc_dend_18S) <- as.character(sample_info_tab$info_ID[order.dendrogram(euc_dend_18S)])
 
 plot(euc_dend_18S, ylab="VST Euc. dist.")
 rel_abund_18S <- apply(count_18S, 2, function(x) x / sum(x))
 bray_dist_18S <- vegdist(t(rel_abund_18S), method = "bray") # Compute Bray–Curtis dissimilarity
 bray_no_salp <- vegdist(t(rel_abund_18S[, df_no_salp$DNA_name]), method = "bray")
 
 # include only organisms with at least 2 pellets and no salp
 table(df_no_salp$Organism)
 organisms_with_reps <- names(which(table(df_no_salp$Organism) >= 2))
 df_replicated <- subset(df_no_salp, Organism %in% organisms_with_reps)
 shared_ids <- intersect(colnames(rel_abund_18S), df_replicated$DNA_name)
 rel_abund_18S_repl <- rel_abund_18S[, shared_ids]
 df_replicated <- df_replicated[df_replicated$DNA_name %in% shared_ids, ]
 bray_repl <- vegdist(t(rel_abund_18S_repl), method = "bray")
 
 adonis2(bray_repl ~ Organism, data = df_replicated)
 betadisper_repl <- betadisper(bray_repl, df_replicated$Organism)
 anova(betadisper_repl)
 
 #
 
 