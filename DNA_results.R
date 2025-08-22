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
 
