# Pi2GAIA (Picrust2 Gene Abundance Informatics Analysis), a workflow for easier visualisation and analyses of PICRUSt2 Output.

Downstream analyses workflow for environmental earth (soil, water etc.) PICRUSt2 KEGG Ortholog output data.

Work in progress.

# Build the script environment by getting the required/recommended tools
############################################################################################################################

# Load in data changing packages
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(tidyverse)
library(data.table)
library(ggforce)
library(ggraph)
library(biomformat)
library(Biostrings)

# Load in analyses/visualising pacakges
library(microbiomeMarker)
library(vegan)
library(phyloseq)
library(ggtree)
library(ape)
library(ggplot2)
library(pathview)
library(pheatmap)
library(plotly)
library(BiocGenerics)


# Make the phyloseq object with feature, sample, tree and refseq data
############################################################################################################################

# Set the working directory
setwd("C:/Users/timwa/Downloads")

# Remove prefixes from KO identifiers in the mapping table
ko2kegg <- read_tsv("C:/Users/timwa/OneDrive/Documenten/ko2kegg_mapping_file.tsv")
ko2kegg$KO_number <- sub("ko:", "", ko2kegg$ko_number)
ko2kegg$pathway <- sub("path:", "", ko2kegg$pathway)

#Load in the metadata
sam_tab <- ("C:/Users/timwa/Downloads/data_farmer_2022.txt")
metadata <- read_tsv("C:/Users/timwa/Downloads/data_farmer_2022.txt")

#load in the required files for 16S data
feature_tab <- ("C:/Users/timwa/Downloads/pred_metagenome_unstrat.tsv")
newick_tree <- "C:/Users/timwa/Downloads/placed_seqs_16S.tre"
refseq_fasta <- "C:/Users/timwa/Downloads/silva_species_assignment_v138.1.fa/silva_species_assignment_v138.1.fa"

# Read the feature table
feature_tab_df <- read_delim(feature_tab)
write_tsv(feature_tab_df, file = "C:/Users/timwa/Downloads/pred_metagenome_16S_unstrat_KO_filtered.tsv") #(16S, ITS, 18S)

# Create a phyloseq object with the filtered feature_table and sample_data
ko_phyloseq <- import_picrust2("pred_metagenome_16S_unstrat_KO_filtered.tsv", sam_tab, trait = "KO") #(16S, ITS, 18S)

# Retrieve KO labels
ko_labels <- rownames(otu_table(ko_phyloseq))

### Normalizing of data was already done in PICRUSt2, however transforming maybe needed to make interpretation possible
### Note: Transform data by using the method that suits the structure and nature of your data (consult literature)

# Log2 transform the abundance data per sample
abundance_matrix <- apply(otu_table(ko_phyloseq), 2, rank)
abundance_matrix <- log2(abundance_matrix + 1)

# Reassign KO labels after transformation
rownames(abundance_matrix) <- ko_labels

# Create a new phyloseq object with the normalized abundance data
ko_phyloseq <- phyloseq(otu_table(as.data.frame(abundance_matrix), taxa_are_rows = TRUE),
                        tax_table(tax_table(ko_phyloseq)),
                        sample_data(metadata))

# Read the Newick tree file and the RefSeq FASTA file
tree <- read.tree(newick_tree)
fasta_sequences <- readDNAStringSet(refseq_fasta)
# Assign the Newick tree and RefSeq FASTA file to the phyloseq object
ko_phyloseq@phy_tree <- tree
ko_phyloseq@refseq <- fasta_sequences


# Perform PCoA and PERMANOVA analyses to study microbiome composition (3D version at end of script for 3 principals)
############################################################################################################################

# Perform PCoA on abundance data using bray distance matrix
ko_data_pcoa <- ordinate(ko_phyloseq, method = "PCoA", distance = "bray", na.rm = TRUE)
#Extract PCoA scores
pcoa_scores <- data.frame(ko_data_pcoa$vectors)

# Add metadata columns to pcoa_scores data frame
pcoa_scores$farm_type <- metadata$farm_type
pcoa_scores$comp_type <- paste(metadata$farm_type, metadata$comp_type, sep = "_")
pcoa_scores$soil_type <- metadata$soil_type
pcoa_scores$sample <- rownames(pcoa_scores)

# Create a scatter plot with points colored by soil_type and shaped by comp_type
pcoa_plot <- ggplot(pcoa_scores, aes(x = Axis.1, y = Axis.2, color = soil_type)) +
  geom_point(alpha = 0.7, size = 3, aes(shape = comp_type)) +
  scale_shape_manual(values = c(19, 17, 1, 2)) + # Circle filled, triangle filled, circle hollow, triangle hollow
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", color = "Soil Type", shape = "Company Practices",
       title = "PCoA Plot of 16S relative KO abundance based on soil type and company practices") +
  theme_minimal() +
  geom_text(aes(label = sample), hjust = -0.1, vjust = -0.5, size = 3) +
  stat_ellipse(aes(group = soil_type), level = 0.95, alpha = 0.5) +
  coord_fixed() +
  theme(legend.position = "right",
        legend.justification = "top",
        legend.box.margin = margin(0, 0, 10, 0))
print(pcoa_plot)

# Create 3D scatter plot with color and symbol mappings
plot_ly(data = pcoa_scores,
        x = ~Axis.1, y = ~Axis.2, z = ~Axis.3,
        color = ~soil_type, symbol = ~comp_type,
        type = "scatter3d", mode = "markers",
        marker = list(size = 5),
        text = ~paste("Sample: ", sample),
        hoverinfo = "text",
        colors = "Dark2",
        symbols = c("circle", "square", "diamond", "triangle-up-open"),
        coloraxis = list(colorbar = list(title = "Soil Type")),
        symbolaxis = list(showlegend = TRUE)) %>%
  layout(scene = list(xaxis = list(title = "PCoA Axis 1"),
                      yaxis = list(title = "PCoA Axis 2"),
                      zaxis = list(title = "PCoA Axis 3"),
                      aspectmode = "cube"),
         title = "PCoA Plot of 16S relative KO abundance (3D)")

# Calculate Bray-Curtis distance matrix directly from the transformed abundance data
dist_matrix <- t(otu_table(ko_phyloseq))
dist_matrix <- vegdist(otu_table(dist_matrix), method = "bray")

# Perform PERMANOVA on distance matrix according to grouping factor
permanova_result <- adonis2(dist_matrix ~ comp_type, data = metadata, permutations = 10000)
print(permanova_result)

#perform regression analysis on the axes and metadata factors to look for what the coordinates in the plot could mean
lm_axis1 <- lm(Axis.1 ~ metadata$pH, data = pcoa_scores)
summary(lm_axis1)

lm_axis2 <- lm(Axis.2 ~ metadata$pH, data = pcoa_scores)
summary(lm_axis2)

# Create a weighted combination of Axis 1 and Axis 2
pcoa_scores$Combined_Axes <- pcoa_scores$Axis.1 + 0.5 * pcoa_scores$Axis.2

# Perform regression analysis using the combined axis
lm_combined <- lm(formula = metadata$pH ~ pcoa_scores$Combined_Axes)
summary(lm_combined)


# Perform pathview analyses for pathway(s) of interest to see which genes are present or absent
############################################################################################################################

# Perform pathway analysis using pathview
pathway_id <- "01053" #enter here the pathway of interest (nutrient metabolism, hormones, 2nd metabolites, stress receptors)
gene_data <- as.matrix(otu_table(ko_phyloseq))
gene_data <- as.matrix(gene_data)
gene_data <- abs(gene_data)
gene_data <- format(gene_data, scientific = FALSE, as.numeric = TRUE)
gene_data <- cbind(ko_number = rownames(gene_data), gene_data)

pathview_result <- pathview(
  gene.data = gene_data,
  pathway.id = pathway_id,
  gene.idtype = "KO",
  species = "ko",
  gene.id = TRUE,
  ko2path = TRUE,
  kegg.native = TRUE)


# Extract the KO numbers/genes if interest from pathview data and/or directly from the abundance data (01053 = siderophores)
############################################################################################################################

# get abundance data and convert it to a data frame
gene_data <- otu_table(ko_phyloseq)
gene_data <- data.frame(gene_data)
# Get siderophores ko numbers and convert it to a data frame
pf_ko_numbers <- pathview_result$plot.data.gene$kegg.names
pf_ko_numbers_df <- data.frame(ko_number = pf_ko_numbers)

# Add ko_numbers as first column
gene_data <- cbind(ko_number = rownames(gene_data), gene_data)
#extract the ko number data with the pathview ko numbers
pf_ko_numbers <- merge(gene_data, pf_ko_numbers_df, by = "ko_number")
# Restore row names from the 'ko_number' column
pf_ko_numbers <- unique(pf_ko_numbers$ko_number)
# Filter gene_data_siderophores to extract abundance data for siderophores_ko_numbers
abundance_matrix <- gene_data[gene_data$ko_number %in% pf_ko_numbers, ]
# Remove the first column (KO numbers) from the filtered data
abundance_matrix <- select(abundance_matrix, -ko_number)

# Replace KO numbers with siderophores genes and descriptions
siderophores_ko_labels <- c("entA: 2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase",
                            "entB: Bifunctional isochorismate lyase / aryl carrier protein",
                            "entC: Isochorismate synthase",
                            "entE: 2,3-dihydroxybenzoate—[aryl-carrier protein] ligase",
                            "entF: L-serine—[L-seryl-carrier protein] ligase",
                            "dhbF: Glycine—[glycyl-carrier protein] ligase",
                            "mbtI: Salicylate synthetase", "mbtA: Mycobactin salicyl-AMP ligase",
                            "mbtB: Mycobactin phenyloxazoline synthetase",
                            "mbtE: Mycobactin peptide synthetase", "mbtC: Mycobactin polyketide synthetase",
                            "mbtD: Mycobactin polyketide synthetase", "mbtF: Mycobactin peptide synthetase",
                            "mbtG: Mycobactin lysine-N-oxygenase", "mbtH: 2,3-dihydroxybenzoate-AMP ligase",
                            "pchF: L-cysteine—[L-cysteinyl-carrier protein] ligase",
                            "asbF: 3-dehydroshikimate dehydratase",
                            "mxcG: Nonribosomal peptide synthetase MxcG", "mxcL: Aminotransferase MxcL")
# Replace row names in abundance_matrix with siderophores_ko_labels
rownames(abundance_matrix) <- siderophores_ko_labels


# Create input for heatmap of nutrient cycles
############################################################################################################################

# Filter the data set based on CHNOPS cycles KO numbers
ko_data_heatmap <- otu_table(ko_phyloseq)
ko_data_heatmap <- as.matrix(ko_data_heatmap)
ko_data_heatmap <- as.data.frame(ko_data_heatmap)
ko_data_heatmap <- cbind(ko_number = rownames(ko_data_heatmap), ko_data_heatmap)

# Define the reordered KO numbers, gene symbols, and gene names
chnops_ko_numbers <- c(
  "K00117",
  "K01077",
  "K01113",
  "K00937",
  "K01093",
  "K02584",
  "K02588",
  "K02591",
  "K02586",
  "K02592",
  "K02585",
  "K00362",
  "K00363",
  "K00367",
  "K00368",
  "K00370",
  "K00371",
  "K00372",
  "K00376",
  "K02305",
  "K02567",
  "K02568",
  "K03385",
  "K11180",
  "K11181",
  "K12339",
  "K17225",
  "K17224",
  "K17218",
  "K02048",
  "K00385",
  "K00392",
  "K17222",
  "K00925",
  "K00625",
  "K01895",
  "K01905",
  "K04073",
  "K00244",
  "K01428",
  "K01429",
  "K01430")

chnops_ko_labels <- c(
  "gcd: Quinoprotein glucose dehydrogenase",
  "phoA: alkaline phosphatase",
  "phoD: alkaline phosphatase D",
  "ppk: Polyphosphate kinase",
  "appA: 4-phytase / acid phosphatase",
  "nifA: Nif-specific regulatory protein",
  "nifH: Nitrogenase iron protein",
  "nifK: Nitrogenase molybdenum-iron protein beta chain",
  "nifD: Nitrogenase molybdenum-iron protein alpha chain",
  "nifN: Nitrogenase molybdenum-iron protein",
  "nifB: Nitrogen fixation protein",
  "nirB: nitrite reductase (NADH) large subunit",
  "nirD: nitrite reductase (NADH) small subunit",
  "narB: ferredoxin-nitrate reductase",
  "nirK: nitrite reductase (NO-forming)",
  "narG/narZ/nxrA: nitrate reductase / nitrite oxidoreductase alpha subunit",
  "narH/narY/nxrB: nitrate reductase / nitrite oxidoreductase beta subunit",
  "nasC/nasA: assimilatory nitrate reductase catalytic subunit",
  "nosZ: nitrous-oxide reductase",
  "norC: nitric oxide reductase subunit C",
  "napA: nitrate reductase (cytochrome)",
  "napB: nitrate reductase (cytochrome) electron transfer subunit",
  "nrfA: nitrite reductase (cytochrome c-552)",
  "dsrA: dissimilatory sulfite reductase alpha subunit",
  "dsrB: dissimilatory sulfite reductase beta subunit",
  "cysM: S-sulfo-L-cysteine synthase (O-acetyl-L-serine-dependent)",
  "soxC: sulfane dehydrogenase subunit SoxC",
  "soxB: S-sulfosulfanyl-L-cysteine sulfohydrolase",
  "sqr: Sulfide:quinone oxidoreductase",
  "asrC: Anaerobic sulfite reductase subunit C",
  "SIR: Sulfite reductase (ferredoxin)",
  "soxA: L-cysteine S-thiosulfotransferase",
  "soxX: L-cysteine S-thiosulfotransferase",
  "ackA: Acetate kinase",
  "pta: Phosphate acetyltransferase",
  "ackdA: Acetate---CoA ligase (ADP-forming) subunit alpha",
  "adhE: Acetaldehyde/alcohol dehydrogenase",
  "ack: Acetate kinase",
  "mhpF: acetaldehyde dehydrogenase",
  "ureC: urease subunit alpha",
  "ureB: urease subunit beta",
  "ureA: urease subunit gamma")

#Extract the abundance data from the phyloseq object relating to nutrient cycling
chnops_ko_data <- ko_data_heatmap[ko_data_heatmap$ko_number %in% chnops_ko_numbers, ]
chnops_ko_data <- as.matrix(chnops_ko_data)
# Extract the abundance matrix
abundance_matrix <- as.matrix(chnops_ko_data[,-1])
# Add ko numbers as row names
rownames(abundance_matrix) <- chnops_ko_data[,1]
# Extract KO numbers and Set the original KO numbers as row names
ko_numbers <- chnops_ko_numbers
rownames(abundance_matrix) <- chnops_ko_labels


# Perform Kruskal Wallis DAA analysis and create heatmap ordered based in grouping factor with selected abundance data
############################################################################################################################

# Convert all elements in the matrix to numeric
saved_row_names <- rownames(abundance_matrix)
abundance_matrix <- apply(abundance_matrix, 2, as.numeric)
rownames(abundance_matrix) <- saved_row_names

# Extract the grouping factor for Kruskal-Wallis test
grouping_factor <- metadata$soil_type
# Perform Kruskal-Wallis test for each row (microbial feature) in the abundance matrix
kruskal_results <- apply(abundance_matrix, 1, function(row) {
  kruskal.test(row ~ grouping_factor)})
print(kruskal_results)
# account for multiple testing with bonferroni
kruskal_results_adjusted <- apply(abundance_matrix, 1, function(row) {
  kruskal_test <- kruskal.test(row ~ grouping_factor)
  p_value <- kruskal_test$p.value
  adjusted_p_value <- p.adjust(p_value, method = "bonferroni")
  return(adjusted_p_value)})
print(kruskal_results_adjusted)

# Create the annotation data frame with the correct sample numbers
annot_data <- data.frame(sample = metadata$sample, soil_type = factor(metadata$soil_type, ordered = TRUE))
abundance_matrix <- t(abundance_matrix)
# Ensure soil_type in annotation data is ordered factor
annot_data$soil_type <- factor(annot_data$soil_type, levels = c("klei", "zand", "zavel"), ordered = TRUE)
# Order the annotation data by soil_type
abundance_matrix <- t(abundance_matrix)
annot_data <- annot_data[order(annot_data$soil_type), ]
# Reorder the abundance_matrix columns based on the soil_type in annotation data
abundance_matrix <- abundance_matrix[, annot_data$sample]
# Create annotation dataframe with row names corresponding to column names of abundance matrix
annotation <- data.frame(soil_type = annot_data$soil_type[match(colnames(abundance_matrix), annot_data$sample)])

# Define annotation colors
annotation_colors <- list(soil_type = c("klei" = "purple", "zand" = "blue", "zavel" = "green"))
# Define a custom color palette with adjusted intensity of red
my_palette <- colorRampPalette(c("#FFFFCC", "#FFD700", "#D6604D"))(20)
# Ensure annotation row names match abundance_matrix column names
rownames(annotation) <- colnames(abundance_matrix)

# Plot the heatmap with correct ordering based on soil_type and annotation colors
pheatmap_plot <- pheatmap(abundance_matrix,
                          cluster_rows = FALSE,
                          cluster_cols = FALSE,
                          display_numbers = FALSE,
                          show_rownames = TRUE,
                          display_colnames = TRUE,
                          annotation_col = annotation,
                          annotation_colors = annotation_colors,
                          color = my_palette,
                          main = "16S relative KO abundance of siderophores synthesis on soil type heatmap")


# Create sankey input data with the stratified output and BIOM file of chosen dataset (16S/ITS/18S)
############################################################################################################################

# Read the BIOM file containing taxonomy information
biom_file <- "C:/Users/timwa/Downloads/asv16stax.biom"
biom <- import_biom(biom_file)

# Read the strat picrust2 output containing asv and ko numbers
stratified_data <- read_tsv("C:/Users/timwa/OneDrive/Documenten/Rstudio/pred_metagenome_contrib_KO_16S.tsv")

# Extract ASV numbers and associated KO numbers
ko_to_asv_mapping <- stratified_data %>%
  group_by(function.) %>%
  summarise(ASVs = paste(taxon, collapse = "|"))

# Separate concatenated ASVs
ko_to_asv_mapping <- ko_to_asv_mapping %>%
  separate_rows(ASVs, sep = "\\|") %>%
  mutate(ASVs = trimws(ASVs))

# Merge KO-to-ASV mapping with taxonomy information
ko_to_asv_taxonomy <- merge(ko_to_asv_mapping, tax_table(biom), by.x = "ASVs", by.y = "row.names", all.x = TRUE)

# Define KO numbers siderophores
ko_numbers_sid <- list(
  "K00216", "K01252", "K01851", "K02363", "K02364", "K04780", "K04781", "K04784", "K04787", "K04788", "K04789", "K04790", "K04791", "K04792", "K04793", "K12240", "K15652", "K15653", "K15681")
# Define KO labels siderophores
ko_labels_sid <- c("entA", "entB", "entC", "entE", "entF", "dhbF", "mbtI", "mbtA", "mbtB", "mbtE", "mbtC", "mbtD", "mbtF", "mbtG", "mbtH", "pchF", "asbF", "mxcG", "mxcL")

# Define KO numbers nutrient cycles
ko_numbers <- list(
  P = c("K00117", "K01077", "K01113", "K00937", "K01093"),
  N = c("K02584", "K02588", "K02591", "K02586", "K02592", "K02585", "K00362", "K00363", "K00367", "K00368", "K00370", "K00371", "K00372", "K00376", "K02305", "K02567", "K02568", "K03385"),
  S = c("K11180", "K11181", "K12339", "K17225", "K17224", "K17218", "K02048", "K00385", "K00392", "K17222"),
  A = c("K00925", "K00625", "K01895", "K01905", "K04072", "K04073", "K00244"),
  U = c("K01428", "K01429", "K01430"))

# Define KO labels corresponding to nutrient KO numbers with symbols only
ko_labels <- list(
  P = c("gcd", "phoA", "phoD", "ppk", "appA"),
  N = c("nifA", "nifH", "nifK", "nifD", "nifN", "nifB", "nirB", "nirD",
        "narB", "nirK", "narG", "narH", "nasC", "nosZ", "norC", "napA", "napB", "nrfA"),
  S = c("dsrA", "dsrB", "cysM", "soxC", "soxB", "sqr", "asrC", "SIR", "soxA", "soxX"),
  A = c("ackA", "pta", "ackdA", "adhE", "adhE", "ack", "mhpF"),
  U = c("ureC", "ureB", "ureA"))

# Make mapping df for nutrient cycling sankey plot
ko_group_map <- data.frame(
  KO = unlist(ko_numbers),
  Group = rep(names(ko_numbers), lengths(ko_numbers)))


# Sankey plot for nutrient cycling groups
############################################################################################################################

# Filter and aggregate ASV contributions by nutrient groups
taxonomy_contributions <- ko_to_asv_taxonomy %>%
  inner_join(ko_group_map, by = c("function." = "KO")) %>%
  group_by(Rank2, Group) %>%
  summarise(Count = n(), .groups = "drop")

# Calculate total contributions per nutrient group
total_contributions <- taxonomy_contributions %>%
  group_by(Group) %>%
  summarise(total_count = sum(Count), .groups = "drop")

# Calculate percentage contribution for each taxon to each nutrient group
taxonomy_contributions <- taxonomy_contributions %>%
  left_join(total_contributions, by = "Group") %>%
  mutate(percent = (Count / total_count) * 100)

# Identify the top contributing taxa for each unique KO number (function.)
top_contributors <- taxonomy_contributions %>%
  group_by(Group, Rank2) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Group) %>%
  arrange(desc(Count)) %>%
  ungroup()

# Filter and extract taxonomy labels for ASVs in top_contributors
taxonomy_labels <- ko_to_asv_taxonomy %>%
  filter(Rank2 %in% unique(top_contributors$Rank2)) %>%
  distinct(Rank2) %>%
  mutate(label = paste(Rank2))

# Filter taxonomy to include only considerable contributions with threshold (1-5>%)
threshold_percent <- 5
significant_contributions <- taxonomy_contributions %>%
  filter(percent >= threshold_percent)

# Create nodes dataframe for left and right sides (taxa and nutrient groups) using significant contributions
left_nodes <- data.frame(name = unique(significant_contributions$Rank2), side = "left")
right_nodes <- data.frame(name = unique(significant_contributions$Group), side = "right")

# Combine nodes dataframe for left and right sides
nodes <- bind_rows(left_nodes, right_nodes) %>%
  arrange(side) %>%
  mutate(id = row_number() - 1)

# Create links dataframe using significant contributions
links <- significant_contributions %>%
  left_join(nodes, by = c("Rank2" = "name")) %>%
  left_join(nodes, by = c("Group" = "name")) %>%
  mutate(
    source = match(Rank2, nodes$name) - 1,
    target = match(Group, nodes$name) + max(source) + 1,
    value = Count,
    label = paste0(round(percent, 2), "%"))

# Define color palette for nodes based on Rank2 taxonomy category
taxon_color_palette <- rainbow(length(unique(nodes$name)))

# Create Sankey plot using plotly with colored links based on Rank2 taxonomy category
sankey_plotly <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = list(
    pad = 15,
    thickness = 20,
    line = list(color = "black", width = 0.5),
    label = nodes$name,
    color = nodes$side),
  link = list(
    source = links$source,
    target = links$target,
    value = links$value,
    label = paste0(links$Group, "<br>", links$label),
    color = taxon_color_palette[as.numeric(factor(links$Rank2))]))

# Customize layout and display the plot
sankey_plotly <- sankey_plotly %>%
  layout(
    title = "16S Taxa Contribution to Nutrient Groups (CHNOPS)",
    font = list(size = 12))
# Display the plot
sankey_plotly


# sankey plot for biosynthesis siderophores
############################################################################################################################

# Filter and aggregate ASV contributions for nutrient cycling KO numbers
taxonomy_contributions <- ko_to_asv_taxonomy %>%
  filter(function. %in% unlist(ko_numbers_sid)) %>%
  group_by(Rank2, function.) %>%
  summarise(Count = n(), .groups = "drop")

# Calculate total contributions per nutrient group
total_contributions <- taxonomy_contributions %>%
  group_by(function.) %>%
  summarise(total_count = sum(Count), .groups = "drop")

# Calculate percentage contribution for each taxon to each KO
taxonomy_contributions <- taxonomy_contributions %>%
  left_join(total_contributions, by = "function.") %>%
  mutate(percent = (Count / total_count) * 100) %>%
  filter(percent > 1)

# Identify the top contributing taxa and their associated KO labels
top_contributors <- taxonomy_contributions %>%
  group_by(Rank2) %>%
  top_n(n = 1, wt = Count) %>%
  ungroup()

# Filter and extract taxonomy labels for ASVs in top_contributors
taxonomy_labels <- ko_to_asv_taxonomy %>%
  filter(Rank2 %in% unique(top_contributors$Rank2)) %>%
  distinct(Rank2) %>%
  mutate(label = paste(Rank2))

# Filter taxonomy to include only considerable contributions with threshold (1-5>%)
threshold_percent <- 5
taxonomy_contributions <- taxonomy_contributions %>%
  filter(percent >= threshold_percent)

# Merge taxonomy labels with top_contributors to create left_nodes
left_nodes <- merge(top_contributors, taxonomy_labels, by = "Rank2", all.x = TRUE) %>%
  select(name = Rank2, label = label)

# Create nodes dataframe for taxa (left side)
left_nodes <- left_nodes %>%
  mutate(side = "left")

# Create right_nodes dataframe for nutrient groups
right_nodes <- bind_rows(
  Map(data.frame, name = unlist(ko_numbers_sid), label = unlist(ko_labels_sid)),
  .id = "group") %>%
  mutate(
    group = gsub("\\d+", "", group),
    side = "right")

# Combine nodes dataframe for left and right sides
nodes <- bind_rows(
  left_nodes,
  right_nodes)
nodes$id <- seq_len(nrow(nodes)) - 1

# Perform left join using dplyr for links
links <- taxonomy_contributions %>%
  left_join(nodes, by = c("Rank2" = "name"), suffix = c(".contrib", ".nodes"), relationship = "many-to-many") %>%
  mutate(
    source = match(Rank2, nodes$name) - 1,
    target = match(function., nodes$name) + max(source) + 1,
    value = Count,
    label = paste0(round(percent, 2), "% ")) %>%
  distinct(source, target, .keep_all = TRUE)

# Define color palette for nodes based on Rank2 taxonomy category
taxon_color_palette <- rainbow(length(unique(nodes$name)))

# Create Sankey diagram using plotly with colored links
sankey_plotly <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = list(
    pad = 15,
    thickness = 20,
    line = list(color = "black", width = 0.5),
    label = nodes$label,
    color = nodes$name),
  link = list(
    source = links$source,
    target = links$target,
    value = links$value,
    label = links$label,
    color = taxon_color_palette[as.numeric(factor(links$Rank2))]))

# Customize layout and display the plot
sankey_plotly <- sankey_plotly %>%
  layout(
    title = "16S taxa contribution in biosynthesis of siderophores",
    font = list(size = 12))
# Display the plot
sankey_plotly










#Extra supplementary analyses used for testing and looking
############################################################################################################################

# Filtering step used to remove duplicate values in ITS and 18S which caused the phyloseq object to fail
# Load in the ITS and 18S abundance data
feature_tab_ITS <- ("C:/Users/timwa/OneDrive/Documenten/Rstudio/pred_metagenome_unstrat_ITS_KO.tsv")
feature_tab_df <- read_delim(feature_tab_ITS)
feature_tab_18S <- ("C:/Users/timwa/OneDrive/Documenten/Rstudio/pred_metagenome_unstrat_18S_KO.tsv")
feature_tab_df <- read_delim(feature_tab_18S)

# Filter out duplicates based on KO column
feature_tab_df <- feature_tab_df %>%
  distinct(function., .keep_all = TRUE)
feature_tab_df <- na.omit(feature_tab_df)
write_tsv(feature_tab_df, file = "C:/Users/timwa/Downloads/pred_metagenome_18S_unstrat_KO_filtered.tsv")
# Continue with creating the phyloseq object with import_picrust2()

#test if data is suitable with quantro
library(quantro)
qtest <- quantro(object = p, groupFactor = pd$CellType)
summary(qtest)
anova(qtest)
quantroStat(qtest)



# Perform NMDS on abundance data using Bray-Curtis distance matrix
ko_data_nmds <- ordinate(ko_phyloseq, method = "NMDS", distance = "bray")

# Extract NMDS scores
nmds_scores <- data.frame(ko_data_nmds$points)

# Add metadata columns to nmds_scores data frame
nmds_scores$farm_type <- metadata$farm_type
nmds_scores$soil_type <- metadata$soil_type
nmds_scores$comp_type <- metadata$comp_type
nmds_scores$sample <- rownames(nmds_scores)

# Create a scatter plot with points colored by soil_type and shaped by comp_type
nmds_plot <- ggplot(nmds_scores, aes(x = MDS1, y = MDS2, color = soil_type, shape = comp_type)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Soil Type", shape = "Component Type",
       title = "NMDS Plot with 2 Axes") +
  theme_minimal() +
  geom_text(aes(label = sample), hjust = -0.1, vjust = -0.5, size = 3) +
  stat_ellipse(aes(group = soil_type), level = 0.95, alpha = 0.5) +
  coord_fixed()
print(nmds_plot)


# Perform PCA on abundance data
ko_data_pca <- otu_table(ko_phyloseq)
ko_data_pca <- data.frame(ko_data_pca)
ko_data_pca <- as.data.frame(t(ko_data_pca))

# Perform PCA on the summarized data and extract PCA scores
pca_result <- prcomp(ko_data_pca, scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)

# Add metadata columns to pca_scores data frame
pca_scores$farm_type <- metadata$farm_type
pca_scores$soil_type <- metadata$soil_type
pca_scores$comp_type <- metadata$comp_type
pca_scores$sample <- metadata$sample

# Create a scatter plot with points colored by soil_type and shaped by comp_type
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = soil_type, shape = comp_type)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(x = "PC1", y = "PC2", color = "Soil Type", shape = "Component Type",
       title = "PCA Plot of ITS ko Abundance with 3 Principal Components") +
  theme_minimal() +
  geom_text(aes(label = sample), hjust = -0.1, vjust = -0.5, size = 3) +
  stat_ellipse(aes(group = soil_type), level = 0.95, alpha = 0.5) +
  coord_fixed()
print(pca_plot)


# Load necessary library
library(kronaTools)

# Convert PICRUSt2 output to Krona format
output_file <- "krona_plot.html"
ktImportText(picrust2_data, output_file)
