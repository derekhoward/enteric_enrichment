# rename all files by placing them in parent folders
# then each file needs to be renamed to standard for reading in 10X data
# REQ: barcodes.tsv, features.tsv, matrix.mtx
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(Seurat)
library(here)

###############################################################################################################################
# load metadata
###############################################################################################################################
meta <- read_delim(here('./data/all.meta.txt'), delim='\t')
meta <- meta[-1,] # drop first row (type info)

gene_universe <- read_delim(here('./data/gene_symbol_annotations/AllenInstitute_custom_Agilent_Array.txt'), delim = '\t')
# some genes are stuck together with sep = ';' 
gene_universe %<>% separate_rows(Gene_symbol, sep=';')
gene_universe <- unique(gene_universe$Gene_symbol)
length(gene_universe)

###############################################################################################################################
# human neurons data
###############################################################################################################################
hum_neur <- Read10X(data.dir = here('./data/hli.neur'), gene.column=1)
hneur_df <- as.data.frame(hum_neur)
hneur_df$gene_symbol <- rownames(hum_neur)
hneur_long <- hneur_df %>% pivot_longer(-gene_symbol, values_to='count', names_to='cell_id')

# get celltype_ids to filter appropriate metadata
hneur_ids <- colnames(hum_neur)
hneur_meta <- meta %>%
  filter(NAME %in% hneur_ids) %>%
  filter(Annotation != 'Neuron') # this removes duplicated metadata from dataset with all cells

# join in metadata and get mean expression for gene per neuron type
hneur_means <- hneur_long %>% 
  left_join(hneur_meta %>% select(NAME, Annotation), by=c('cell_id' = 'NAME')) %>% 
  group_by(Annotation, gene_symbol) %>% 
  summarise(meanExp = mean(count)) %>% 
  filter(gene_symbol %in% gene_universe)

hneur_means %<>% rename(cell_type = Annotation, expression = meanExp)


###############################################################################################################################
# human glia data
###############################################################################################################################
hum_glia <- Read10X(data.dir = here('./data/hli.glia'), gene.column=1)
hglia_df <- as.data.frame(hum_glia)
hglia_df$gene_symbol <- rownames(hum_glia)
hglia_long <- hglia_df %>% pivot_longer(-gene_symbol, values_to='count', names_to='cell_id')

# extract cell type annotations from tsne data
hglia_tsne <- read_delim('./data/hli.glia/tsne2.txt', delim = '\t')
hglia_tsne <- hglia_tsne[-1,]
hglia_tsne[,2:3] <- lapply(hglia_tsne[,2:3], as.numeric)


hglia_means <- hglia_long %>% 
  left_join(hglia_tsne %>% select(NAME, LABEL), by=c('cell_id' = 'NAME')) %>% 
  group_by(LABEL, gene_symbol) %>% 
  summarise(meanExp = mean(count)) %>% 
  filter(gene_symbol %in% gene_universe)

hglia_means %<>% rename(cell_type = LABEL, expression = meanExp)

###############################################################################################################################
# now that celltype means has been calculated for glia and neurons above
# we should do the log(1+exp) and concatenate the results for the different celltypes together

# then do the zscore for genes across all celltypes
# then use the final ranks for shiny app
###############################################################################################################################
dim(hneur_means)
dim(hglia_means)
length(unique(hneur_means$gene_symbol))
length(unique(hglia_means$gene_symbol))

overlapping_genes <- intersect(hneur_means$gene_symbol, hglia_means$gene_symbol)
length(overlapping_genes)

hneur_means %<>% 
  mutate(log1Exp = log(1+expression))

hglia_means %<>% 
  mutate(log1Exp = log(1+expression))

ent_cells <- bind_rows(hneur_means, hglia_means)
ent_cells %<>% filter(gene_symbol %in% overlapping_genes)

# write raw expression to file
ent_cells %>%
  select(-expression) %>% 
  pivot_wider(names_from='cell_type', values_from='log1Exp') %>% 
  write_csv(here('./data/processed/ent_neurons_glia_logExp.csv'))


normalized_expression <- ent_cells %>% 
  group_by(gene_symbol) %>% 
  mutate(log1ExpZ = (log1Exp - mean(log1Exp)) / sd(log1Exp)) %>% 
  select(-expression, -log1Exp)

# write expression zscores to file
normalized_expression %>% 
  pivot_wider(names_from='cell_type', values_from='log1ExpZ') %>% 
  write_csv(here('./data/processed/ent_neurons_glia_zscores.csv'))
  
normalized_expression %<>%
  group_by(cell_type) %>% 
  mutate(log1ExpZRank = rank(log1ExpZ)) %>%
  select(-log1ExpZ)
final_ranks <- normalized_expression %>% pivot_wider(names_from='cell_type', values_from='log1ExpZRank')

dir.create(here('./data/processed'))
write_csv(final_ranks, here('./data/processed/ent_neurons_glia_ranks.csv'))

#############################################################################################
# Processing dataset with all human data
#############################################################################################
human <- Read10X(data.dir = './data/hli', gene.column=1)

genes_assayed <- rownames(human)
class(genes_assayed)
overlap <- intersect(genes_assayed, gene_universe)

length(overlap)
# select only overlap between gene_universe and genes_assayed 
human[overlap,]
#############################################################################################
# exploring working on chunks of data
meta %>% group_by(Dataset) %>% tally()

# select just human data for all cells
human_celltypes <- meta %>%
  filter(Dataset == 'Human colon all cells (10X)') %>% 
  select(Annotation) %>% unique() %>% .$Annotation

human_meta <- meta %>%
  filter(Dataset == 'Human colon all cells (10X)')

human_meta %>% 
  group_by(Annotation) %>% 
  tally() %>% 
  arrange(desc(n))

#############################################################################################
# Process all celltypes in chunks
#############################################################################################
celltype_means <- list()
for (celltype in human_celltypes) {
  ids <- human_meta %>% filter(Annotation == celltype) %>% .$NAME
  data_subset <- human[overlap, ids]
  subset_mean <- Matrix::rowMeans(data_subset)
  
  celltype_means[[celltype]] <- subset_mean
}

celltype_means <- as.data.frame(celltype_means)
celltype_means$gene_symbol <- rownames(celltype_means) 
celltype_means <- as_tibble(celltype_means) %>% select(gene_symbol, everything())

#############################################################################################
celltype_means_long <- celltype_means %>% 
  pivot_longer(-gene_symbol, values_to='meanExp', names_to='cell_type')

celltype_means_long %<>% mutate(log1Exp = log(meanExp+1))

celltype_means_long %>% 
  select(-meanExp) %>% 
  pivot_wider(names_from='cell_type', values_from='log1Exp') %>% 
  write_csv(here('./data/processed/all_ent_logExp.csv'))

normalized_expression <- celltype_means_long %>% 
  group_by(gene_symbol) %>% 
  mutate(log1ExpZ = (log1Exp - mean(log1Exp)) / sd(log1Exp)) %>% 
  select(-meanExp, -log1Exp)

final_ranks <- normalized_expression %>%
  group_by(cell_type) %>% 
  mutate(log1ExpZRank = rank(log1ExpZ)) %>%
  select(-log1ExpZ) %>% 
  pivot_wider(names_from='cell_type', values_from='log1ExpZRank')

dir.create(here('./data/processed'))
write_csv(final_ranks, here('./data/processed/all_ent_ranks.csv'))
