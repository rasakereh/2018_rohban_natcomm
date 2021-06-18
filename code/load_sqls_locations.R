library(dplyr)
library(magrittr)
library(cytominer)
library(foreach)
library(stringr)
library(readr)
library(doParallel)
library(rdist)
library(infotheo)
library(kit)

bug.detected <- c("BRD-K83508485-001-01-9", "BRD-K83597974-003-05-7", "BRD-K83637872-001-01-7", "BRD-K84175871-003-02-2", "BRD-K84214706-001-02-4", "BRD-K84266862-001-01-6", "BRD-K84421793-001-01-1", "BRD-K84595254-001-03-0", "BRD-K84663978-001-02-3", "BRD-K84709232-001-02-6", "BRD-K85015012-003-01-1", "BRD-K85090592-003-08-7", "BRD-K85104575-001-03-2", "BRD-K85119730-001-06-5", "BRD-K85242180-001-01-1", "BRD-K85266041-304-01-7", "BRD-K85333151-003-04-3", "BRD-K85383046-001-02-5", "BRD-K85555482-001-01-1", "BRD-K85871428-001-01-9", "BRD-K85872723-001-02-7", "BRD-K85880973-001-02-7", "BRD-K85883481-001-02-6", "BRD-K86204871-001-02-2", "BRD-K86301799-001-14-0", "BRD-K86434416-001-02-7", "BRD-K84663978-001-03-1", "BRD-K86509404-001-02-0", "BRD-K86600316-003-01-2", "BRD-K86727142-001-06-6", "BRD-K86873305-236-03-0", "BRD-K86882815-001-01-6", "BRD-K86958018-001-01-9", "BRD-K87048468-003-01-9", "BRD-K87049188-001-03-6", "BRD-K87156652-001-05-1", "BRD-K87194840-001-03-0", "BRD-K87349602-001-01-2", "BRD-K87492696-001-04-1", "BRD-K87510569-001-01-0", "BRD-K87696786-003-01-0", "BRD-K87798455-001-02-6", "BRD-K87817668-001-02-7", "BRD-K87904882-001-02-3", "BRD-K87905482-001-02-4", "BRD-K87990216-001-03-7", "BRD-K87991767-001-02-0", "BRD-K88090157-050-03-3", "BRD-K88172511-310-03-8", "BRD-K88219015-001-02-5", "BRD-K88282786-066-04-9", "BRD-K88358234-003-01-4", "BRD-K88366685-300-03-7", "BRD-K88405679-003-01-5", "BRD-K88429204-001-04-7", "BRD-K88544581-001-01-2", "BRD-K88568253-065-02-1", "BRD-K88611939-001-01-8", "BRD-K88646909-004-03-0", "BRD-K88677950-001-01-3", "BRD-K88682005-001-04-2", "BRD-K88741031-001-01-0", "BRD-K88743730-001-01-0", "BRD-K88759641-003-01-0", "BRD-K88789588-001-03-2", "BRD-K88809146-003-02-1", "BRD-K88849294-001-02-1", "BRD-K88871508-001-01-9", "BRD-K88954184-001-02-2", "BRD-K89046952-001-03-4", "BRD-K89055274-048-03-7", "BRD-K89085489-001-02-6", "BRD-K89125793-001-05-4", "BRD-K89152108-001-03-3", "BRD-K89210380-001-03-8", "BRD-K89225759-001-02-1", "BRD-K89271983-001-01-5", "BRD-K89272762-001-04-4", "BRD-K89274813-001-01-5")

df2numeric <- function(df)
{
  if(nrow(df) == 1)
  {
    return(apply(df, 2, as.numeric) %>% t() %>% as.data.frame())
  }
  apply(df, 2, as.numeric) %>% as.data.frame()
}

profile.plate.location <- function(pl, project.name, batch.name, n.components = 3000, rand.density = 0.1, correlation, out.path, in.path = NULL, in.type="sqlite", cores = 4, nrm.column = "Metadata_broad_sample", nrm.value = "DMSO", feat.list = NULL) {
    
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
    
    if(dt.sub[1,"Metadata_broad_sample"] %in% bug.detected)
    {
      print('bad well')
      print(dt.sub[1,"Metadata_broad_sample"])
      write.csv(paste0('../', dt.sub[1,"Metadata_broad_sample"], '.csv'))
    }
   
    cell.count = nrow(dt.sub)
    nth.topn = as.integer(cell.count**.5)
    
    if(cell.count==1){
      dt.sub[, variables] <- as.list(apply(dt.sub[, variables], 2, function(x) as.numeric(x)))
    }else{
      dt.sub[, variables] <- apply(dt.sub[, variables], 2, function(x) as.numeric(x))
    }
    
    pos.features <- c('Cells_AreaShape_Center_X', 'Cells_AreaShape_Center_Y')
    features <- setdiff(all.variables, pos.features)
    cell.dists <- pdist(dt.sub[pos.features], metric='euclidean')
    threshold <- cell.dists %>% apply(2, function(col.data){
      if(length(col.data[col.data!=0]) == 0) return(0)
      nth.index = (col.data[col.data!=0] %>% topn(nth.topn, decreasing=F))[nth.topn]
      col.data[col.data!=0][nth.index]
    })
    
    indices <- which(0 < cell.dists & cell.dists < threshold)
    valid.rows.all <- ((indices-1) %/% cell.count) + 1
    valid.cols.all <- ((indices-1) %% cell.count) + 1
    indices <- indices[valid.rows.all < valid.cols.all]
    valid.rows <- valid.rows.all[valid.rows.all < valid.cols.all]
    valid.cols <- valid.cols.all[valid.rows.all < valid.cols.all]
    cell.dists <- cell.dists[indices]

    diff <- abs(df2numeric(dt.sub[valid.rows,features]) - df2numeric(dt.sub[valid.cols,features]))
    
    location.info <- cor(diff, cell.dists) %>% abs() %>% t() %>% as.data.frame()
    location.info[is.na(location.info)] <- 0
    location.info <- cbind(location.info, dt.sub[1, metadata])
    
    location.info <- location.info %>%
      select(one_of(c(metadata, variables)))
    
    # if(!(dt.sub.saved))
    # {
    #   write.csv(dt.sub, '../dt_sub.csv')
    #   write.csv(location.info, '../location_info.csv')
    #   dt.sub.saved <<- TRUE
    # }
    profile <- location.info
    profile <- cbind(profile, data.frame(Metadata_Plate = pl, Metadata_Well = sites))
    print("new dim:")
    print(dim(profile))
    print("stepped")
    profile
  }
  t2 <- proc.time()
  print("profiles formed")
  
  cov.variables <- colnames(profiles)
  cov.variables <- cov.variables[which(str_detect(cov.variables, "Cells_") | str_detect(cov.variables, "Cytoplasm_") | str_detect(cov.variables, "Nuclei_"))]
  
  cov.metadata <- setdiff(colnames(profiles), cov.variables)
  
  dmso.ids <- prf %>% 
    dplyr::filter((!!rlang::sym(nrm.column)) == nrm.value) %>%
    select(Metadata_Plate, Metadata_Well)
  
  dt <- profiles[, cov.variables]
  mn <- apply(dt, 2, function(x) mean(x, na.rm=T))
  sdv <- apply(dt, 2, function(x) sd(x, na.rm=T))
  dt.nrm <- scale(dt, center = mn, scale = sdv)
  print("normalized")
  
  profiles.nrm <- cbind(dt.nrm, profiles[, cov.metadata])
  
  profiles.nrm.meta <- profiles.nrm[, cov.metadata] %>% 
    mutate(Metadata_Plate = as.character(Metadata_Plate)) %>%
    left_join(prf %>% select(one_of(prf.metadata)) %>%
                mutate(Metadata_Plate = as.character(Metadata_Plate)), by = c("Metadata_Plate", "Metadata_Well"))
  
  profiles.nrm <- cbind(profiles.nrm[, cov.variables], profiles.nrm.meta)
  
  cov.metadata <- setdiff(colnames(profiles.nrm), cov.variables)
  
  readr::write_csv(profiles.nrm, paste0(out.path, "/", pl, "_location.csv"))
  
  if(is.null(in.path))
  {
    system(paste0(paste0("rm ", sql.path)))
  }
}
