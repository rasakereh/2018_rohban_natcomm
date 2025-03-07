library(dplyr)
library(magrittr)
library(cytominer)
library(foreach)
library(stringr)
library(readr)
library(doParallel)

source("generate_component_matrix.R")

profile.plate <- function(pl, project.name, batch.name, n.components = 3000, rand.density = 0.1, out.path, in.path = NULL, in.type="sqlite", cores = 4, nrm.column = "Metadata_broad_sample", nrm.value = "DMSO", feat.list = NULL) {
    
  doParallel::registerDoParallel(cores = cores)
  
  sql.path <- ifelse(is.null(in.path), paste0("../backend/", batch.name, "/", pl, "/", pl, ".sqlite"), paste0(in.path, "/", pl, "/", pl, ".sqlite"))
  out.file <- paste0(out.path, "/", pl, ".csv")
  
  if (!file.exists(sql.path)) {
    system(command = paste0("aws s3 cp 's3://cellpainting-datasets/",
                            project.name, 
                            "/workspace/backend/",
                            batch.name,
                            "/", 
                            pl, 
                            "/", 
                            pl, 
                            ".sqlite' ",
                            sql.path))
  }
  
  norm.path = ifelse(is.null(in.path) , paste0("../input/", pl, "_normalized.csv"), paste0(in.path, "/", pl, "/", pl, "_normalized.csv"))
  
  if (!file.exists(norm.path)) {
    system(command = paste0("aws s3 cp 's3://cellpainting-datasets/",
                            project.name, 
                            "/workspace/backend/",
                            batch.name,
                            "/", 
                            pl, 
                            "/", 
                            pl, 
                            "_normalized.csv' ", norm.path))
  }
  
  prf <- readr::read_csv(norm.path)
  
  sites.all <- unique(prf$Metadata_Well)
  
  variables <- colnames(prf)
  variables <- variables[which(!str_detect(variables, "Metadata_"))]
  prf.metadata <- setdiff(colnames(prf), variables)
  if (!is.null(feat.list)) {
    variables <- feat.list
  }
  
  sqlite_file <- sql.path
  db <- DBI::dbConnect(RSQLite::SQLite(), sqlite_file)
  RSQLite::initExtension(db)
  
  image_object_join_columns <- c("TableNumber", "ImageNumber")
  strata <- c("Image_Metadata_Plate", "Image_Metadata_Well")
  
  if(in.type == "sqlite"){
    image <- dplyr::tbl(src = db, "image") %>%
      dplyr::select(all_of(c(image_object_join_columns, strata)))
  }else{
    root.path = paste0(in.path, "/", pl, "/")
    image <- read.csv(paste0(root.path, "Image.csv")) %>%
      dplyr::select(all_of(c(image_object_join_columns, strata)))
    cells <- read.csv(paste0(root.path, "Cells.csv"))
    cytoplasm <- read.csv(paste0(root.path, "Cytoplasm.csv"))
    nuclei <- read.csv(paste0(root.path, "Nuclei.csv"))
  }
  
  image.coll <- image %>%
    dplyr::select(Image_Metadata_Well, ImageNumber, TableNumber)
  if(in.type == "sqlite")
  {
    image.coll <- image.coll %>% dplyr::collect()
  }
  
  t <- proc.time()
  
  profiles <- foreach (sites = sites.all, .combine = rbind) %do% {
    
    #saveRDS(sites, paste0("../tmp_trad/", sites, ".rds"))  
    
    image.sub <- image.coll %>% 
      dplyr::filter(Image_Metadata_Well == sites) 
    
    append_operation_tag <- function(s) stringr::str_c(s, operation, sep = "_")
    
    dt.sub <- foreach (i = 1:NROW(image.sub), .combine = rbind) %dopar% {
      if(in.type == "sqlite"){
        dbPar <- DBI::dbConnect(RSQLite::SQLite(), sqlite_file)
        RSQLite::initExtension(dbPar)
        image <- dplyr::tbl(src = dbPar, "Image") %>%
          dplyr::select(all_of(c(image_object_join_columns, strata)))
        cells <- dplyr::tbl(src = dbPar, "Cells")
        cytoplasm <- dplyr::tbl(src = dbPar, "Cytoplasm")
        nuclei <- dplyr::tbl(src = dbPar, "Nuclei")
      }

      cells.sub <- cells %>% 
        dplyr::filter(ImageNumber == !!image.sub$ImageNumber[i] & TableNumber == !!image.sub$TableNumber[i])
      
      cytoplasm.sub <- cytoplasm %>% 
        dplyr::filter(ImageNumber == !!image.sub$ImageNumber[i] & TableNumber == !!image.sub$TableNumber[i])
      
      nuclei.sub <- nuclei %>% 
        dplyr::filter(ImageNumber == !!image.sub$ImageNumber[i] & TableNumber == !!image.sub$TableNumber[i])
      
      if(in.type == "sqlite")
      {
        cells.sub <- cells.sub %>% dplyr::collect()
        cytoplasm.sub <- cytoplasm.sub %>% dplyr::collect()
        nuclei.sub <- nuclei.sub %>% dplyr::collect()
        DBI::dbDisconnect(dbPar)
      }
      
      dt <- cells.sub %>%
        dplyr::inner_join(cytoplasm.sub, by = "ObjectNumber") %>%
        dplyr::inner_join(nuclei.sub, by = "ObjectNumber") %>%
        dplyr::inner_join(image.sub, by = image_object_join_columns) 
      
      dt
    }
    
    all.variables <- colnames(dt.sub)
    all.variables <- all.variables[which(str_detect(all.variables, "Cells") | str_detect(all.variables, "Cytoplasm") | str_detect(all.variables, "Nuclei"))]
    metadata <- setdiff(colnames(dt.sub), all.variables)
    
    dt.sub <- dt.sub %>%
      select(one_of(c(metadata, variables)))
    
    if(nrow(dt.sub)==1){
      dt.sub[, variables] <- as.list(apply(dt.sub[, variables], 2, function(x) as.numeric(x)))
    }else{
      dt.sub[, variables] <- apply(dt.sub[, variables], 2, function(x) as.numeric(x))
    }
    
    profile <- cytominer::covariance(population = dt.sub, variables = variables)
    profile <- cbind(profile, data.frame(Metadata_Plate = pl, Metadata_Well = sites))
  }
  t2 <- proc.time()
  
  cov.variables <- colnames(profiles)
  cov.variables <- cov.variables[which(str_detect(cov.variables, "Cells_") | str_detect(cov.variables, "Cytoplasm_") | str_detect(cov.variables, "Nuclei_"))]
  saveRDS(cov.variables, "../tmp/cov_var.rds")
  
  cov.metadata <- setdiff(colnames(profiles), cov.variables)
  
  dmso.ids <- prf %>% 
    dplyr::filter((!!rlang::sym(nrm.column)) == nrm.value) %>%
    select(Metadata_Plate, Metadata_Well)
  
  samples.nrm <- profiles %>%
    mutate(Metadata_Plate = as.character(Metadata_Plate)) %>%
    semi_join(dmso.ids %>%
                mutate(Metadata_Plate = as.character(Metadata_Plate)),
              by = c("Metadata_Plate", "Metadata_Well"))
  
  mn <- apply(samples.nrm %>% select(one_of(cov.variables)), 2, function(x) mean(x, na.rm = T))
  sdv <- apply(samples.nrm %>% select(one_of(cov.variables)), 2, function(x) sd(x, na.rm = T))
  
  dt <- profiles[, cov.variables]
  # mn <- apply(dt, 2, function(x) mean(x, na.rm=T))
  # sdv <- apply(dt, 2, function(x) sd(x, na.rm=T))
  
  dt.nrm <- scale(dt, center = mn, scale = sdv)
  
  profiles.nrm <- cbind(dt.nrm, profiles[, cov.metadata])
  
  profiles.nrm.meta <- profiles.nrm[, cov.metadata] %>% 
    mutate(Metadata_Plate = as.character(Metadata_Plate)) %>%
    left_join(prf %>% select(one_of(prf.metadata)) %>%
                mutate(Metadata_Plate = as.character(Metadata_Plate)), by = c("Metadata_Plate", "Metadata_Well"))
  
  profiles.nrm <- cbind(profiles.nrm[, cov.variables], profiles.nrm.meta)
  
  cov.metadata <- setdiff(colnames(profiles.nrm), cov.variables)
  
  rds.path = paste0(out.path, "/random_projection_unified.rds")
  if (!file.exists(rds.path)) {
    rand.proj <- generate_component_matrix(n_features = length(cov.variables), n_components = n.components, density = rand.density)
    saveRDS(rand.proj, rds.path)
  } else {
    rand.proj <- readRDS(rds.path)
  }
  
  dt <- as.matrix(profiles.nrm[, cov.variables])
  dt[is.na(dt)] <- 0
  
  profiles.nrm.red <- dt %*% as.matrix(rand.proj)
  profiles.nrm <- cbind(as.data.frame(profiles.nrm.red), profiles.nrm[, cov.metadata])
  
  readr::write_csv(profiles.nrm, paste0(out.path, "/", pl, "_covariance.csv"))
  
  if(is.null(in.path))
  {
    system(paste0(paste0("rm ", sql.path)))
  }
}
