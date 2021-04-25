#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends

'profile
Usage:
profile_trad -n <project_name> -b <batch_name> -t <input_type> -i <input_dir> -r <output_dir> -p <plate_number> -o <operation> -l <norm_column> -v <norm_value> -c <no_cores> [-f <feat_list_file>]

Options:
-h --help                                         Show this screen.
-n <project_name> --name=<project_name>           Project name on s3.
-b <batch_name> --batch=<batch_name>              Batch name. 
-t <input_type> --type=<input_type>               accepts "csv" and "sqlite",
-i <input_dir> --input=<input_dir>                Path to the directory containing plates,
-r <output_dir> --result=<output_dir>             Path to the result directory. 
-p <plate_number> --plate=<plate_number>          Plate number.
-o <operation> --operation=<operation>            Profiling summary function; e.g. mean, mean+sd, etc. 
-f <feat_list_file> --feats=<feat_list_file>      Path to the file containing the list of features. 
-l <norm_column> --col=<norm_column>              Column name to be used to select samples for normalization.
-v <norm_value> --value=<norm_value>              Value of the mentioned column which indicates the sample.
-c <no_cores> --cores=<no_cores>                  Number of cores to be used. 
' -> doc

opts <- docopt::docopt(doc)

source("load_sqls_mean_sd_2.R")
library(readr)
library(cytominer)
library(dplyr)
library(stringr)

pl <- opts[["plate"]]
proj.name <- opts[["name"]]
batch.name <- opts[["batch"]]
input.type <- opts[["type"]]
input.dir <- opts[["input"]]
result.dir <- opts[["result"]]
feat.list <- opts[["feats"]]
operation <- opts[["operation"]]
col.name <- opts[["col"]]
col.val <- opts[["value"]]
cores <- opts[["cores"]]

if (!is.null(feat.list)) {
  feat.list <- readr::read_csv(feat.list, col_names = F)  
  feat.list <- unname(unlist(feat.list))
}

if (!is.null(result.dir)) {
  outp.dir <- paste0(result.dir, "/", batch.name, "/", pl)
}else{
  outp.dir <- paste0("../backend/", batch.name, "/", pl)
}

if(!dir.exists(outp.dir)) {
  dir.create(outp.dir, recursive = T)
}

if(is.null(input.type) || !any(input.type %in% c("sqlite", "csv")))
{
  print("Invalid/Unspecified input type, using `sqlite`")
  input.type <- "sqlite"
}

profile.plate.traditional.2(pl = pl, project.name = proj.name, batch.name = batch.name, operation = operation, nrm.column = col.name, nrm.value = col.val, out.path = outp.dir, in.path = input.dir, in.type = input.type, cores = cores, feat.list = feat.list)