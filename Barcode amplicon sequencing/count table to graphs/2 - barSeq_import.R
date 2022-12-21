#this script will import a csv/tsv file containing a list of barcodes and will generate a table with the count of each barcode. 

#####parameters####
input_dir <- "inputs/PAT_FS"
output_dir <- "outputs/PAT_FS"
files_id <- "cluster.csv" #a string containing a piece of text that should be present in the name of the files to import.
input_type <- "read_count" #options are "read_count" or "barcode_list". "read_count" must be a table with at least two columns,
#one with the barcodes and the other with the read count. "barcode_list must be a table with at least 1 column containing the barcodes. 
barcode_col <- "Center" #name or index of the column containing the barcodes
count_col <- "time_point_1" #only used if input type is "read_count"
sample_info <- "inputs/PAT_FS/sample_info.xlsx"
min_count <- 2 #barcodes with read count lower than this number will be dropped. 
add_missing_barcodes <- FALSE #if TRUE, will add to each sample the barcodes found in all other samples...
#...These will receive a count value and percent value of 0. #Useful for plots, but can lead to very large objects.
save_r_object <- TRUE #if true, will save the list of inported barcodes as and R object in a folder called "outputs"

#libraries
library(readr)
library(tidyverse)

BClist <- import_bc_lists(path = input_dir, 
                          pattern = files_id,
                          input_type = input_type,
                          sample_info = sample_info,
                          barcode_col = barcode_col,
                          count_col = count_col)

if(add_missing_barcodes){
  BClist <- add_missing_bc(BClist)
}

if(save_r_object){
  dir.create(output_dir)
  saveRDS(BClist, paste0(output_dir, "/BClist.rds"))
}

#cleaning environment
rm(files)
rm(files_id)
rm(input_dir)
rm(input_type)
rm(row)
#to do: if a file is not in the sample_info, return an error.

