# rename all files by placing them in parent folders
# then each file needs to be renamed to standard for reading in 10X data
# REQ: barcodes.tsv, features.tsv, matrix.mtx
library(readr)
library(dplyr)
library(magrittr)
library(Seurat)
library(here)

#############################################################################################
# meta 
#############################################################################################
gene_universe <- read_delim(here('./data/gene_symbol_annotations/AllenInstitute_custom_Agilent_Array.txt'), delim = '\t')
# some genes are stuck together with sep = ';' 
library(tidyr)
gene_universe %<>% separate_rows(Gene_symbol, sep=';')
gene_universe <- unique(gene_universe$Gene_symbol)
length(gene_universe)

meta <- read_delim('./data/all.meta.txt', delim='\t')
meta <- meta[-1,] # drop first row (type info)
#############################################################################################
# all human data
human <- Read10X(data.dir = './data/hli', gene.column=1)


genes_assayed <- rownames(human)
class(genes_assayed)
overlap <- intersect(genes_assayed, gene_universe)
length(overlap)

human[overlap,]

#############################################################################################
# Figure out sampling of celltypes to work with a smaller dataset
#############################################################################################
#filtered_meta <- meta %>% filter(Dataset %in% 'Human colon all cells (10X)')
#filtered_meta %>% select(Annotation, Dataset) %>% distinct()
#sampled_meta <- filtered_meta %>% sample_n(1500)
#sort(unique(sampled_meta$Annotation))

#human_celltypes <- unique(filtered_meta$Annotation)
#sort(human_celltypes)
# so the sampled_meta has at least one of each of the celltypes!

#sampled_cell_ids <- sampled_meta$NAME
#############################################################################################
# select sample of data and restrict genes to universe of interest
#############################################################################################
#human_sample <- human[overlap,sampled_cell_ids]
human_sample <- human[overlap,]
dim(human_sample)

hdf <- as.data.frame(human_sample)
hdf$gene_symbol <- rownames(human_sample)
hdf_long <- hdf %>% pivot_longer(-gene_symbol, values_to='count', names_to='cell_id')

# extract cell type annotations from tsne data
h_tsne <- read_delim('./data/hli/tsne2.txt', delim = '\t')
h_tsne <- h_tsne[-1,]
h_tsne[,2:3] <- lapply(h_tsne[,2:3], as.numeric)

h_means <- hdf_long %>% 
  left_join(h_tsne %>% select(NAME, LABEL), by=c('cell_id' = 'NAME')) %>% 
  group_by(LABEL, gene_symbol) %>% 
  summarise(meanExp = mean(count))

h_means %<>% rename(cell_type = LABEL, expression = meanExp)


###############################################################################################################################
# now that celltype means has been calculated for glia and neurons above
# we should do the log(1+exp) and concatenate the results for the different celltypes together

# then do the zscore for genes across all celltypes
# then use the final ranks for shiny app
###############################################################################################################################
dim(h_means)
length(unique(h_means$gene_symbol))

h_means %<>% 
  mutate(log1Exp = log(1+expression))

# write log(expression+1) to file
h_means %>% 
  select(-expression) %>% 
  pivot_wider(names_from='cell_type', values_from='log1Exp') %>% 
  write_csv(here('./data/processed/all_ent_cells_logExp.csv'))

normalized_expression <- h_means %>% 
  group_by(gene_symbol) %>% 
  mutate(log1ExpZ = (log1Exp - mean(log1Exp)) / sd(log1Exp)) %>% 
  select(-expression, -log1Exp)

normalized_expression %>% 
  filter(is.na(log1ExpZ)) %>%
  #filter(!is.na(log1ExpZ)) %>% 
  dim()

# write expression zscores to file
normalized_expression %>% 
  pivot_wider(names_from='cell_type', values_from='log1ExpZ') %>% 
  write_csv(here('./data/processed/all_ent_zscores.csv'))

normalized_expression %<>%
  group_by(cell_type) %>% 
  mutate(log1ExpZRank = rank(log1ExpZ)) %>%
  select(-log1ExpZ)
final_ranks <- normalized_expression %>% pivot_wider(names_from='cell_type', values_from='log1ExpZRank')
write_csv(final_ranks, here('./data/processed/all_ent_ranks.csv'))


