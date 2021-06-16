library(dplyr)
library(ggplot2)
library(mixOmics)
library(FactoMineR)
library(nnTensor)
library(pheatmap)
source("moa_evaluations.R")

k.snf <- 7     # neighborhood size in SNF
t <- 10
k <- 1:10      # k top hits are used for classification
sigma <- .5

mean.na <- function(x) {mean(x, na.rm = T)}

do.PCA <- function(mat, final.dim)
{
  pca.res <- prcomp(mat)
  transformed <- pca.res$x
  as.data.frame(transformed[,1:final.dim])
}

read.and.summarize <- function(profile_dir, plate.list, feat.list, profile.type, metadata.df) {
  fls <- paste0(plate.list, "_covariance.csv")
  if(!is.null(profile_dir))
  {
    fls <- paste0(plate.list, "/", fls)
  }
  feat.list.s <- feat.list
  
  profiles.nrm <- foreach (fl = fls, .combine = rbind) %do% {
    if (profile.type == "cov") {
      cov.profile = ifelse(is.null(profile_dir), paste0("../output/", fl), paste0(profile_dir, fl))
      if (file.exists(cov.profile)) {
        x <- readr::read_csv(cov.profile)    
      } else {
        x <- NULL
      }
    } else if (profile.type == "mean") {
      pl <- str_split(fl, "_")[[1]][1]
      mean.profile = ifelse(is.null(profile_dir), paste0("../input/", pl, "_normalized.csv"), paste0(profile_dir, pl, "_normalized.csv"))
      x <- readr::read_csv(mean.profile)  
      if (!is.null(feat.list)) {
        x <- x %>%
          dplyr::select(matches("Metadata_"), one_of(feat.list))
      }
    } else if (profile.type == "median") {
      pl <- str_split(fl, "_")[[1]][1]
      if(is.null(profile_dir))
      {
        init <- list.dirs("../backend", recursive = F)
        fl.name <- paste0(init, "/", pl, "/", pl, "_normalized_median_mad.csv")
      }else{
        fl.name <- paste0(profile_dir, "/", pl, "_normalized_median_mad.csv")
      }
      
      if (file.exists(fl.name)) {
        x <- readr::read_csv(fl.name)    
      } else {
        x <- NULL
        warning(paste0("Plate ", pl, " is missing."))
      }
      
      if (!is.null(feat.list) & ! is.null(x)) {
        x <- x %>%
          dplyr::select(matches("Metadata_"), one_of(paste0(feat.list, "_median")))
        feat.list.s <- paste0(feat.list, "_median")
      }
    } else if (profile.type == "mad") {
      pl <- str_split(fl, "_")[[1]][1]
      if(is.null(profile_dir))
      {
        init <- list.dirs("../backend", recursive = F)
        fl.name <- paste0(init, "/", pl, "/", pl, "_normalized_median_mad.csv")
      }else{
        fl.name <- paste0(profile_dir, "/", pl, "_normalized_median_mad.csv")
      }
      
      if (file.exists(fl.name)) {
        x <- readr::read_csv(fl.name)    
      } else {
        x <- NULL
        warning(paste0("Plate ", pl, " is missing."))
      }
      
      if (!is.null(feat.list) & ! is.null(x)) {
        x <- x %>%
          dplyr::select(matches("Metadata_"), one_of(paste0(feat.list, "_mad")))
        feat.list.s <- paste0(feat.list, "_mad")
      }
    } else if (profile.type == "location") {
      pl <- str_split(fl, "_")[[1]][1]
      if(is.null(profile_dir))
      {
        init <- list.dirs("../backend", recursive = F)
        fl.name <- paste0(init, "/", pl, "/", pl, "_location.csv")
      }else{
        fl.name <- paste0(profile_dir, "/", pl, "_location.csv")
      }
      
      if (file.exists(fl.name)) {
        x <- readr::read_csv(fl.name)    
      } else {
        x <- NULL
        warning(paste0("Plate ", pl, " is missing."))
      }
      
      if (!is.null(feat.list) & ! is.null(x)) {
        x <- x %>%
          dplyr::select(matches("Metadata_"), one_of(feat.list))
        feat.list.s <- feat.list
      }
    }
    else if (str_detect(profile.type, "\\+")) {
      p1 <- str_split(profile.type, "\\+")[[1]][1]  
      p2 <- str_split(profile.type, "\\+")[[1]][2]  
      
      if (!is.null(feat.list)) {
        p1 <- str_split(profile.type, "\\+")[[1]][1]  
        p2 <- str_split(profile.type, "\\+")[[1]][2]  
        feat.list.s <- c(paste0(feat.list, "_", p1), paste0(feat.list, "_", p2))
      }
      
      pl <- str_split(fl, "_")[[1]][1]
      if(is.null(profile_dir))
      {
        init <- list.dirs("../backend", recursive = F)
        fl.name <- paste0(init, "/", pl, "/", pl, "_normalized_", p1, "_", p2, ".csv")
      }else{
        fl.name <- paste0(profile_dir, "/", pl, "_normalized_", p1, "_", p2, ".csv")
      }
      if (file.exists(fl.name)) {
        x <- readr::read_csv(fl.name)    
      } else {
        x <- NULL
        warning(paste0("Plate ", pl, " is missing."))
      }
    }
    x
  }
  
  variable.names <- colnames(profiles.nrm)
  variable.names <- variable.names[which(!str_detect(variable.names, "Metadata_"))]
  
  # in some special cases of mean profiles (e.g. CDRP), there seems to be Inf, and NA value.
  # this is to treat those cases
  if (profile.type != "cov") {
    meta.cols <- setdiff(colnames(profiles.nrm), variable.names)
    
    if (is.null(feat.list)) {
      ids <- apply(profiles.nrm[,variable.names], 2, function(x) !any(is.na(x) | is.nan(x) | is.infinite(x) | sd(x) > 10)) %>% which
      variable.names <- variable.names[ids]
    } else {
      variable.names <- feat.list.s
    }
    
    profiles.nrm <- profiles.nrm %>% dplyr::select(one_of(c(meta.cols, variable.names)))
  }
  
  print(length(fls))
  
  if (!"Metadata_mmoles_per_liter" %in% colnames(profiles.nrm)) {
    profiles.nrm <- profiles.nrm %>%
      mutate(Metadata_mmoles_per_liter = 10)
  }
  
  prf <- profiles.nrm %>% 
    group_by(Metadata_broad_sample, Metadata_mmoles_per_liter, Metadata_Plate_Map_Name) %>%
    summarise_at(.vars = variable.names, .funs = "mean.na")
  
  prf <- prf %>%
    arrange(abs(Metadata_mmoles_per_liter - 10)) %>%
    group_by(Metadata_broad_sample) %>%
    slice(1) %>%
    ungroup()
  
  if (is.null(metadata.df)) {
    prf <- prf %>% 
      left_join(profiles.nrm %>% dplyr::select(Metadata_broad_sample, Metadata_moa) %>% unique, by = "Metadata_broad_sample")
  } else {
    prf <- prf %>% 
      left_join(metadata.df %>% dplyr::select(Metadata_broad_sample, Metadata_moa) %>% unique, by = "Metadata_broad_sample")
  }
  
  profiles.nrm <- prf
  feats <- colnames(prf)
  feats <- feats[which(!str_detect(feats, "Metadata_"))]
  
  return(list(data = profiles.nrm, feats = feats))
}

profiles2melt <- function(pf, profile.type)
{
  profiles.nrm <- Pf$data
  feats <- Pf$feats
  if (!is.null(metadata.df)) {
    profiles.meta <- metadata.df %>% dplyr::select("Metadata_broad_sample", "Metadata_moa") %>% unique
  } else {
    profiles.meta <- profiles.nrm %>% dplyr::select("Metadata_broad_sample", "Metadata_moa") %>% unique
  }
  
  pm <- profiles.nrm %>% dplyr::select(Metadata_broad_sample, Metadata_Plate_Map_Name) %>% unique 
  profiles.meta <- profiles.meta %>% left_join(pm, by = "Metadata_broad_sample")
  
  cr <- cor(profiles.nrm[, feats] %>% t)
  rownames(cr) <- profiles.nrm$Metadata_broad_sample
  colnames(cr) <- profiles.nrm$Metadata_broad_sample
  
  cr.melt <- cr %>% reshape2::melt()
  
  cr.melt <- cr.melt %>% left_join(profiles.meta, by = c("Var1" = "Metadata_broad_sample")) %>% left_join(profiles.meta, by = c("Var2" = "Metadata_broad_sample")) 
  
  print(paste0("cr_",  profile.type, ifelse(profile.type == "mix", paste0(mix1, mix2), ""), ".rds"))
  saveRDS(cr.melt, paste0("cr_",  profile.type, ifelse(profile.type == "mix", paste0(mix1, mix2), ""), ".rds"))
}

fuse.matrices <- function(matrix.list, method)
{
  affinity.result <- NULL
  transformed <- NULL
  # method must be one of: pseudo-PFA, block.pls, rgcca, MFA, jNMF, iNMF, SNF
  if(method == 'pseudo-PFA')
  {
    comp.cnt <- c(sapply(matrix.list, ncol), 100) %>% min
    transformed <- lapply(matrix.list, function(mat){cat('first phase PCA... '); do.PCA(mat, comp.cnt)})
    print('binding first phase results...')
    transformed <- do.call(cbind, transformed)
    transformed <- do.PCA(transformed, 100)
  }else if(method == 'block.pls')
  {
    
  }else if(method == 'rgcca')
  {
    design <- matrix(1, nrow=length(matrix.list), ncol=length(matrix.list))
    transformed <- wrapper.rgcca(X = matrix.list, design = design, tau = rep(1, length(matrix.list)), ncomp = 5, scheme = "centroid")
    transformed <- do.call(cbind, transformed$variates)
  }else if(method == 'MFA')
  {
    groups <- sapply(matrix.list, ncol)
    name.group <- matrix.list %>% names()
    transformed <- do.call(cbind, matrix.list)
    transformed <- MFA(transformed, group=groups, ncp=100, name.group=name.group, graph=F)
    transformed <- transformed$global.pca$ind$coord
  }else if(method == 'jNMF')
  {
    transformed <- nnTensor::jNMF(matrix.list, J=min(100, nrow(matrix.list[[1]])-1))$W
  }else if(method == 'iNMF')
  {
    
  }else if(method == 'SNF')
  {
    affinities <- lapply(matrix.list, function(mat){
      cor.mat <- mat %>% t %>% cor %>% sim_normalize
      d <- apply(cor.mat, 1, function(x) !(sum(is.na(x)) >= (NROW(cor.mat) -1 )))
      cor.mat <- cor.mat[d, d]
      cor.mat[is.na(cor.mat)] <- 0
      SNFtool::affinityMatrix(Diff = 1 - cor.mat, K = k.snf, sigma = sigma)
    })
    affinity.result <- SNFtool::SNF(affinities, K = k.snf, t = round(3/2 * t), auto_stop=T)
  }
  if(method != 'SNP' & !is.null(transformed))
  {
    names(transformed) <- make.names(names(transformed))
    cor.mat <- transformed %>% scale %>% t %>% cor %>% sim_normalize
    affinity.result <- SNFtool::affinityMatrix(Diff = 1 - cor.mat, K = k.snf, sigma = sigma)
  }
  affinity.result
}

generate_data <- function(rows, cols, rand_gens)
{
  controls <- sort(sample(rows, rows %/% 2))
  cases <- sort(setdiff(1:rows, controls))
  result <- lapply(rand_gens, function(random.gen){
    data1 <- as.data.frame(matrix(get(random.gen)(cols*rows), nrow=rows))
    data1[controls,] <- apply(data1[controls,], 1, function(row){sort(row, decreasing=F)})
    data1[cases,] <- apply(data1[cases,], 1, function(row){sort(row, decreasing=T)})
    scale(data1)
  })
  names(result) <- names(rand_gens)
  data.labels <- rep('case', rows)
  data.labels[controls] <- 'control'
  data.labels <- as.factor(data.labels)
  list(data = result, labels = data.labels)
}

visualize <- function(dataset, labels=rep(1, nrow(dataset)), type='dataset')
{
  if(type == 'dataset')
  {
    dataset <- scale(dataset)
    reduced <- do.PCA(dataset, 2)
    ggplot(reduced, aes(x=PC1, y=PC2, colour=labels)) + geom_point(size=10) + theme_bw()
  }else if(type == 'affinity')
  {
    colnames(dataset) <- as.character(labels)
    rownames(dataset) <- as.character(labels)
    pheatmap(dataset, border_color = NA)
  }
}


# rand_gens <- c(a='runif', b='rexp', c='rnorm')
# sim.data <- generate_data(380, 1800, rand_gens)
# matrix.list <- sim.data$data
# sim.labels <- sim.data$labels

# method must be one of: pseudo-PFA, block.pls, rgcca, MFA, jNMF, iNMF, SNF
# for(method.name in c("pseudo-PFA", "MFA", "jNMF", "SNF", "rgcca")){
#   png(paste0('../../results/test/', method.name, '.png'), height=800, width=800)
#   print(method.name)
#   fused.mat <- fuse.matrices(matrix.list, method.name)
#   visualize(fused.mat, sim.labels, 'affinity')
#   dev.off()
# }

