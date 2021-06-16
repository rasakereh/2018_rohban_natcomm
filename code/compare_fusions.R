#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'compare_fusions
Usage:
compare_fusion -t <perturbation_type> -p <plate_list_path> [-d <profiles_dir> -e <metadata_file> -f <feat_list_file>]
Options:
-h --help                                         Show this screen.
-t <perturbation_type> --pert=<perturbation_type> Either chemical or genetic
-d <profiles_dir> --profiles=<profiles_dir>       The parent directory which contains profiles
-e <metadata_file> --metadata=<metadata_file>     Path to a csv file containing the association between the Metadata_broad_sample and Metadata_moa. This could be skipped if it is present in the profiles.
-f <feat_list_file> --feats=<feat_list_file>      Path to a text file containing the list of features to be used for the mean profiles.
-p <plate_list_path> --plates=<plate_list_path>   Path to the plate list text file.
' -> doc

opts <- docopt::docopt(doc)

p <- as.character(opts[["pert"]])
meta.file <- opts[["metadata"]]
feat.list <- opts[["feats"]]
plates <- opts[["plates"]]
profile_dir <- opts[["profiles"]]

library(dplyr)
library(ggplot2)

source("fusion.R")

genetic <- (p == "genetic")
not.same.batch <- !genetic

if (genetic) {
  not.same.batch <- F
}

plate.list <- readr::read_csv(plates, col_names = F) %>% as.matrix() %>% as.vector()
if (!is.null(meta.file)) {
  metadata.df <- readr::read_csv(meta.file)  
} else {
  metadata.df <- NULL
}
if (!is.null(feat.list)) {
  feat.list <- readr::read_csv(feat.list, col_names = F)  
  feat.list <- unname(unlist(feat.list))
}

profile.types <- c('median', 'mad', 'cov', 'location')

print('Loading data...')
# whole.data <- lapply(profile.types, function(profile.type){
#   print(paste("loading", profile.type))
#   read.and.summarize(profile_dir, plate.list, feat.list, profile.type, metadata.df)
# })

# all.feats <- lapply(whole.data, function(dataset) {dataset$feats})
# feats <- all.feats[[1]]
# meta_feats <- whole.data[[1]]$data %>% colnames %>% setdiff(feats)
# metadata.cols <- lapply(whole.data, function(dataset) dataset$data[,meta_feats])
# sample.names <- metadata.cols[[1]] %>% dplyr::select(Metadata_broad_sample)

# print('Imputing missing data...')
# whole.data <- lapply(seq_along(whole.data), function(wholedata, name, index){
#   dataset <- wholedata[[index]]$data[,all.feats[[index]]]
#   for(i in 1:ncol(dataset)){
#     dataset[is.na(dataset[,i]), i] <- mean(dataset[,i], na.rm = TRUE)
#   }
#   dataset
# }, wholedata=whole.data, name=names(whole.data))

# saveRDS(whole.data, 'wholedata.rds')

whole.data <- readRDS('wholedata.rds')

fusion.methods <- c("pseudo-PFA", "MFA", "jNMF", "SNF", "rgcca")

print('Fusing median, mad, cov., and loc.')
affinities.loc <- lapply(fusion.methods, function(fusion.method){
  print(paste("Fusing datasets using", fusion.method))
  affinity.matrix <- fuse.matrices(whole.data, fusion.method)
  rownames(affinity.matrix) <- sample.names
  colnames(affinity.matrix) <- sample.names
})
names(affinities.loc) <- paste('median+mad+cov.+loc. using', fusion.methods)

print('Fusing median, mad and cov.')
location.index <- which(profile.types == 'location')
affinities.no.loc <- lapply(fusion.methods[-location.index], function(fusion.method){
  print(paste("Fusing datasets using", fusion.method))
  affinity.matrix <- fuse.matrices(whole.data, method.name)
  rownames(affinity.matrix) <- sample.names
  colnames(affinity.matrix) <- sample.names
})
names(affinities.no.loc) <- paste('median+mad+cov. using', fusion.methods)

affinities <- c(affinities.loc, affinities.no.loc)

metadata <- metadata.cols[[1]] %>%
  dplyr::select(Metadata_broad_sample, Metadata_moa, Metadata_Plate_Map_Name)

metadata <- metadata %>% mutate(Metadata_moa = str_to_lower(Metadata_moa))

#cmpd_classification <- Vectorize(cmpd_classification, "k0")
cmpd_knn_classification <- Vectorize(cmpd_knn_classification, "k0")

print('Passing cmpd_knn_classification phase ...')
d.affinities <- lapply(affinities, function(cor.mat){
  cmpd_knn_classification(cor.mat, metadata, k, not.same.batch = not.same.batch) 
})

l.affinities <- lapply(d.affinities, function(d.affinity){
  lapply(d.affinity[3, ], function(x) sum(x))
})

print('Generating classification comparison plot ...')
dataset.methods <- names(affinities)
D <- data.frame(method = dataset.methods[[1]], k = k, tp = (unlist(l.affinities[[1]])))
for(dataset.method in dataset.methods[-1])
{
  D <- rbind(D, 
        data.frame(method = dataset.method, k = k, tp = (unlist(l.affinities[[dataset.method]]))))
}
D <- D %>% mutate(method = factor(method, levels = dataset.methods))

g <- ggplot(D, aes(x = k, y = tp, color = method, order = as.character(method))) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = k, minor_breaks = k) +
  ylab("No. of treatments") + 
  ggtitle("No. of treatments with a \n relevant match in their k-NNs") +
  theme_bw() +
  theme(axis.text = element_text(size=20), text = element_text(size=15)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.title=element_blank())
print(g) 
ggsave("classification_comparison.png", g, width = 7, height = 5)

top.prec <- c(seq(from = 0.98, to = 0.997, by = 0.002))
enrichment_top_conn <- Vectorize(enrichment_top_conn, vectorize.args = "top.perc")

print('Calculating enrichment scores...')
sm.affinities <- lapply(affinities, function(cor.mat){
  perpare_sm(sm = cor.mat, metadata = metadata)
})

affinity.results <- lapply(sm.affinities, function(sm.cor.mat){
  enrichment_top_conn(sm = sm.median.mad, metadata = metadata, top.perc = top.prec, not.same.batch = not.same.batch)[3, ] %>% unlist %>% unname()
})

print('Generating enrichment plot ...')
dataset.methods <- names(affinities)
D <- data.frame(top.prec = top.prec * 100, odds.ratio = affinity.results[[1]], method = dataset.methods[[1]])
for(dataset.method in dataset.methods[-1])
{
  D <- rbind(D, data.frame(top.prec = top.prec * 100, odds.ratio = affinity.results[[dataset.method]], method = dataset.method))
}
D <- D %>% mutate(method = factor(method, levels = dataset.methods))
D <- D %>% mutate(top.prec = 100 - top.prec)

g <- ggplot(D, aes(x = top.prec, y = odds.ratio, color = method, order = method)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_continuous(breaks = 100 - rev(top.prec[seq(from = 1, to = length(top.prec), by = 2)] * 100), minor_breaks = 100 - rev(top.prec * 100)) +
  ylab("Folds of enrichment") + 
  xlab("p") +
  ggtitle("Folds of enrichment for top p% connections \n to have same MOAs/Pathways") +
  theme_bw() +
  theme(axis.text = element_text(size=20), text = element_text(size=15)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.title=element_blank())
print(g) 
ggsave("global_comparison.png", g, width = 7, height = 5)
save.image("workspace.RData")
