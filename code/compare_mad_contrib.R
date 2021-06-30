#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'compare_mean_cov
Usage:
compare_mean_cov -p <perturbation_type>
Options:
-h --help                                         Show this screen.
-p <perturbation_type> --pert=<perturbation_type> Either chemical or genetic' -> doc

opts <- docopt::docopt(doc)

p <- as.character(opts[["pert"]])

library(dplyr)
library(ggplot2)

source("moa_evaluations.R")

enrichment.based.classification <- FALSE
k.snf <- 7     # neighborhood size in SNF
t <- 10
k <- 1:10      # k top hits are used for classification
genetic <- (p == "genetic")
not.same.batch <- !genetic
snf.med.mad <- T

if (genetic) {
  not.same.batch <- F
}

cr.melt.mean <- readRDS("cr_median.rds")
cr.melt.cov <- readRDS("cr_cov.rds")
cr.melt.mad <- readRDS("cr_mad.rds")  

sigma.mean <- 0.5 
sigma.cov <- 0.5 
sigma.mad <- 0.5 
sigma.loc <- 0.5 

cr.mean <- cr.melt.mean %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2")  

cr.mad <- cr.melt.mad %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2") 

cr.cov <- cr.melt.cov %>%
  select(Var1, Var2, value) %>%
  group_by(Var1, Var2) %>%
  summarise(value = max(value)) %>%
  reshape2::acast("Var1 ~ Var2")

cr.loc <- matrix(runif(nrow(cr.mad)*ncol(cr.mad)))
rownames(cr.loc) <- rownames(cr.mad)
colnames(cr.loc) <- colnames(cr.mad)

cr.mean <- sim_normalize(cr.mean) 
cr.cov <- sim_normalize(cr.cov) 
cr.loc <- sim_normalize(cr.loc) 
cr.mad <- sim_normalize(cr.mad) 

d <- apply(cr.mean, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.mean) -1 )))
cr.mean <- cr.mean[d, d]

d <- apply(cr.cov, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.cov) -1 )))
cr.cov <- cr.cov[d, d]

d <- apply(cr.mad, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.mad) -1 )))
cr.mad <- cr.mad[d, d]

d <- apply(cr.loc, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.loc) -1 )))
cr.loc <- cr.loc[d, d]

cm.rn <- setdiff(intersect(intersect(rownames(cr.mad), rownames(cr.mean)), intersect(rownames(cr.cov), rownames(cr.loc))), NA)
cr.mean <- cr.mean[cm.rn, cm.rn]
cr.cov <- cr.cov[cm.rn, cm.rn]
cr.mad <- cr.mad[cm.rn, cm.rn]
cr.loc <- cr.loc[cm.rn, cm.rn]
  
cr.mean[is.na(cr.mean)] <- 0
cr.loc[is.na(cr.loc)] <- 0
cr.mad[is.na(cr.mad)] <- 0
cr.cov[is.na(cr.cov)] <- 0

# median, median+cov, median+mad+cov
# median+garbage, median+cov+garbage, median+mad+cov+garbage

af.snf <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
cr.median <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.3 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = sigma.cov)
af.snf <- SNFtool::SNF(list(af.1, af.3), K = k.snf, t = round(3/2 * t), auto_stop=T)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.cov <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mad, K = k.snf, sigma = sigma.mad)
af.3 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = sigma.cov)
af.snf <- SNFtool::SNF(list(af.1, af.2, af.3), K = k.snf, t = round(3/2 * t), auto_stop=T)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.mad.cov <- af.snf

####################################### garbages

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.3 <- SNFtool::affinityMatrix(Diff = 1 - cr.loc, K = k.snf, sigma = sigma.loc)
af.snf <- SNFtool::SNF(list(af.1, af.3), K = k.snf, t = round(3/2 * t), auto_stop=T)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.gbg <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.loc, K = k.snf, sigma = sigma.loc)
af.3 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = sigma.cov)
af.snf <- SNFtool::SNF(list(af.1, af.2, af.3), K = k.snf, t = round(3/2 * t), auto_stop=T)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.cov.gbg <- af.snf

af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = sigma.mean)
af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mad, K = k.snf, sigma = sigma.mad)
af.3 <- SNFtool::affinityMatrix(Diff = 1 - cr.cov, K = k.snf, sigma = sigma.cov)
af.4 <- SNFtool::affinityMatrix(Diff = 1 - cr.loc, K = k.snf, sigma = sigma.loc)
af.snf <- SNFtool::SNF(list(af.1, af.2, af.3, af.4), K = k.snf, t = round(3/2 * t), auto_stop=T)
rownames(af.snf) <- rownames(af.1)
colnames(af.snf) <- colnames(af.1)
cr.median.mad.cov.gbg <- af.snf

metadata <- cr.melt.mean %>%
  select(Var1, Metadata_moa.x, Metadata_Plate_Map_Name.x) %>%
  unique() %>%
  mutate(Metadata_broad_sample = Var1, Metadata_moa = Metadata_moa.x, Metadata_Plate_Map_Name = Metadata_Plate_Map_Name.x) %>%
  select(-Var1, -Metadata_moa.x, -Metadata_Plate_Map_Name.x)

metadata <- metadata %>% mutate(Metadata_moa = str_to_lower(Metadata_moa))

#cmpd_classification <- Vectorize(cmpd_classification, "k0")
cmpd_knn_classification <- Vectorize(cmpd_knn_classification, "k0")

if (enrichment.based.classification) {
  d.mean <- cmpd_classification(sm = cr.mean, metadata = metadata, k0 = k)
  d.mix <- cmpd_classification(sm = cr.mix, metadata = metadata, k0 = k)
  
  if (length(k) == 1) {
    d.mean <- (d.mean %>% t) %>% apply(., 2, function(x) {y <- data.frame(as.vector(x)); names(y) <- names(x); y}) 
    d.mix <- (d.mix %>% t) %>% apply(., 2, function(x) {y <- data.frame(as.vector(x)); names(y) <- names(x); y}) 
    
    d.mean <- d.mean %>% as.data.frame()
    d.mix <- d.mix %>% as.data.frame()
    
    d.mean.sel <- d.mean %>%
      filter(pass)
    
    d.mix.sel <- d.mix %>%
      filter(pass)
    
    d.mean.sel %>% NROW() %>% print
    d.mix.sel %>% NROW() %>% print
    
    diff.cmpds <- setdiff(d.mix.sel$Var1, d.mean.sel$Var1)
    d.mix %>%
      filter(Var1 %in% diff.cmpds) %>%
      select(-Metadata_moa.x, -pass) %>%
      left_join(d.mean, by = "Var1") %>%
      select(-pass) %>%
      rename(p.val.mix = p.val.x, p.val.mean = p.val.y, Compound = Var1) %>%
      arrange(p.val.mix) %>%
      knitr::kable() %>%
      print()
    
    diff.cmpds <- setdiff(d.mean.sel$Var1, d.mix.sel$Var1)
    d.mean %>%
      filter(Var1 %in% diff.cmpds) %>%
      select(-Metadata_moa.x, -pass) %>%
      left_join(d.mix, by = "Var1") %>%
      select(-pass) %>%
      rename(p.val.mix = p.val.y, p.val.mean = p.val.x, Compound = Var1) %>%
      arrange(p.val.mean) %>%
      knitr::kable() %>%
      print()
  } else {
    l.mean <- lapply(d.mean["pass", ], function(x) sum(x)) 
    l.mix <- lapply(d.mix["pass", ], function(x) sum(x))
    
    D <- data.frame(method = "mean", k = k, tp = (unlist(l.mean)))
    D <- rbind(D, 
               data.frame(method = "mean+cov.", k = k, tp = (unlist(l.mix))))
    g <- ggplot(D, aes(x = k, y = tp, color = method)) + 
      geom_point() + 
      geom_line() + 
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_continuous(breaks = k, minor_breaks = k) +
      ylab("No. of treatment with a \n same MOA/Pathway treatment in their k-NNs")
    print(g)    
    ggsave("classification_comparison.png", g, width = 7, height = 5)
  }
} else {
  d.median <- cmpd_knn_classification(cr.median, metadata, k, not.same.batch = not.same.batch) 
  d.median.cov <- cmpd_knn_classification(cr.median.cov, metadata, k, not.same.batch = not.same.batch) 
  d.median.mad.cov <- cmpd_knn_classification(cr.median.mad.cov, metadata, k, not.same.batch = not.same.batch) 
  d.median.gbg <- cmpd_knn_classification(cr.median.gbg, metadata, k, not.same.batch = not.same.batch) 
  d.median.cov.gbg <- cmpd_knn_classification(cr.median.cov.gbg, metadata, k, not.same.batch = not.same.batch) 
  d.median.mad.cov.gbg <- cmpd_knn_classification(cr.median.mad.cov.gbg, metadata, k, not.same.batch = not.same.batch) 

  if (length(k) == 1) {
    d.mean %>% NROW %>% print
    d.mix %>% NROW %>% print
    d.diff <- d.mean %>%
      full_join(d.mix, by = "Var1") %>%
      mutate(Metadata_moa = ifelse(is.na(Metadata_moa.x.x), Metadata_moa.x.y, Metadata_moa.x.x)) %>%
      select(-Metadata_moa.x.x, -Metadata_moa.x.y) %>%
      mutate(pass.mean = ifelse(is.na(pass.x), "No", "Yes"),
             pass.mix = ifelse(is.na(pass.y), "No", "Yes")) %>%
      select(-pass.x, -pass.y) %>% 
      arrange(pass.mix) 
    
    d.diff %>%  
      filter(pass.mix == "Yes" & pass.mean == "No") %>%
      arrange(Metadata_moa) %>%
      htmlTable::htmlTable() %>%
      print
    
    d.diff %>%  
      filter(pass.mix == "No" & pass.mean == "Yes") %>%
      arrange(Metadata_moa) %>%
      htmlTable::htmlTable() %>%
      print
  } else {
    l.median <- lapply(d.median[3, ], function(x) sum(x))
    l.median.cov <- lapply(d.median.cov[3, ], function(x) sum(x))
    l.median.mad.cov <- lapply(d.median.mad.cov[3, ], function(x) sum(x))
    l.median.gbg <- lapply(d.median.gbg[3, ], function(x) sum(x))
    l.median.cov.gbg <- lapply(d.median.cov.gbg[3, ], function(x) sum(x))
    l.median.mad.cov.gbg <- lapply(d.median.mad.cov.gbg[3, ], function(x) sum(x))
    
    D <- data.frame(method = "median", k = k, tp = (unlist(l.median)))
    D <- rbind(D, 
               data.frame(method = "median+cov.", k = k, tp = (unlist(l.median.cov))))
    D <- rbind(D, 
               data.frame(method = "median+mad+cov.", k = k, tp = (unlist(l.median.mad.cov))))
    D <- rbind(D, 
               data.frame(method = "median+garbage", k = k, tp = (unlist(l.median.gbg))))
    D <- rbind(D, 
               data.frame(method = "median+cov.+garbage", k = k, tp = (unlist(l.median.cov.gbg))))
    D <- rbind(D, 
               data.frame(method = "median+mad+cov.+garbage", k = k, tp = (unlist(l.median.mad.cov.gbg))))

    lvls <- c("median","median+cov.","median+mad+cov.","median+garbage","median+cov.+garbage","median+mad+cov.+garbage")
    D <- D %>% mutate(method = factor(method, levels = lvls))
    
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
  }
}

top.prec <- c(seq(from = 0.98, to = 0.997, by = 0.002))
enrichment_top_conn <- Vectorize(enrichment_top_conn, vectorize.args = "top.perc")

sm.median <- perpare_sm(sm = cr.median, metadata = metadata)
sm.median.cov <- perpare_sm(sm = cr.median.cov, metadata = metadata)
sm.median.mad.cov <- perpare_sm(sm = cr.median.mad.cov, metadata = metadata)
sm.median.gbg <- perpare_sm(sm = cr.median.gbg, metadata = metadata)
sm.median.cov.gbg <- perpare_sm(sm = cr.median.cov.gbg, metadata = metadata)
sm.median.mad.cov.gbg <- perpare_sm(sm = cr.median.mad.cov.gbg, metadata = metadata)

# saveRDS(sm.median, "sm_median.rds")
# saveRDS(sm.median.cov, "sm_median_cov.rds")
# saveRDS(sm.median.mad.cov, "sm_median_mad_cov.rds")
# saveRDS(sm.median.gbg, "sm_median_gbg.rds")
# saveRDS(sm.median.cov.gbg, "sm_median_cov_gbg.rds")
# saveRDS(sm.median.mad.cov.gbg, "sm_median_mad_cov_gbg.rds")

median.res <- enrichment_top_conn(sm = sm.median, metadata = metadata, top.perc = top.perc, not.same.batch = not.same.batch)
median.cov.res <- enrichment_top_conn(sm = sm.median.cov, metadata = metadata, top.perc = top.perc, not.same.batch = not.same.batch)
median.mad.cov.res <- enrichment_top_conn(sm = sm.median.mad.cov, metadata = metadata, top.perc = top.perc, not.same.batch = not.same.batch)
median.gbg.res <- enrichment_top_conn(sm = sm.median.gbg, metadata = metadata, top.perc = top.perc, not.same.batch = not.same.batch)
median.cov.gbg.res <- enrichment_top_conn(sm = sm.median.cov.gbg, metadata = metadata, top.perc = top.perc, not.same.batch = not.same.batch)
median.mad.cov.gbg.res <- enrichment_top_conn(sm = sm.median.mad.cov.gbg, metadata = metadata, top.perc = top.perc, not.same.batch = not.same.batch)

median.res <- median.res[3, ] %>% unlist %>% unname()
median.cov.res <- median.cov.res[3, ] %>% unlist %>% unname()
median.mad.cov.res <- median.mad.cov.res[3, ] %>% unlist %>% unname()
median.gbg.res <- median.gbg.res[3, ] %>% unlist %>% unname()
median.cov.gbg.res <- median.cov.gbg.res[3, ] %>% unlist %>% unname()
median.mad.cov.gbg.res <- median.mad.cov.gbg.res[3, ] %>% unlist %>% unname()

D5 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.res, method = "median")
D6 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.cov.res, method = "median+cov.")
D7 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.cov.res, method = "median+mad+cov.")
D8 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.gbg.res, method = "median+garbage")
D9 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.cov.gbg.res, method = "median+cov.+garbage")
D10 <- data.frame(top.prec = top.prec * 100, odds.ratio = median.mad.cov.gbg.res, method = "median+mad+cov.+garbage")

D <- rbind(D5, D6, D7, D8, D9, D10)

#lvls <- sort(unique(as.character(D$method)))
lvls <- c("median","median+cov.","median+mad+cov.","median+garbage","median+cov.+garbage","median+mad+cov.+garbage")
D <- D %>% mutate(method = factor(method, levels = lvls))
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
