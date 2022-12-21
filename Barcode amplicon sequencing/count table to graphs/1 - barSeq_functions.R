#import_bc_files####
#this function will inport barcode count tables into a single list of BClist objects.
####inporting files####
import_bc_lists <- function(path, pattern, input_type, sample_info, barcode_col = 1, count_col = 2, min_count = 0){
  library(readr)
  library(readxl)
  #inv_grep is a subfunction that will search a vector for elements that partially match a give pattern
  inv_grep <- function(vector, pattern){
    indices <- c()
    for(i in c(1:length(vector))){
      if(grepl(vector[i], pattern)){
        indices <- c(indices, i)
      }
    }
    return(indices)
  }
  #getting files
  files <- list.files(path = path,
                      pattern = pattern)
  
  #reading samples metadata: to do - make it optional
  sample_info <- read_excel(sample_info)
  
  print(paste("found", length(files), "files"))
  
  BClist <- vector(mode = "list", length = length(files))
  for(i in c(1:length(files))){
    if(input_type == "barcode_list"){
      bc_table <- read_csv(paste(path, files[i], sep = "/"), col_names = TRUE)
      colnames(bc_table)[1] <- "barcode"
      bc_table <- bc_table$barcode
      bc_table <- table(bc_table)
    }else{
      bc_table <- read_csv(paste(path, files[i], sep = "/"))
    }
    bc_table <- data.frame(bc_table[,c(barcode_col, count_col)])
    colnames(bc_table) <- c("barcode", "count")
    bc_table <- filter(bc_table, count >= min_count)
    bc_table$percent <- bc_table$count/sum(bc_table$count)
    bc_table <- arrange(bc_table, desc(percent))
    bc_table$position <- c(1:nrow(bc_table))
    list <- list(bc_table = bc_table, file = files[i])
    #appending information from the sample_info table
    row <- inv_grep(sample_info$id, files[i])
    for(j in c(1:ncol(sample_info))){
      variable <- colnames(sample_info)[j]
      value <- as.character(sample_info[row,j])
      list <- append(list, value)
      names(list)[length(names(list))] <- variable
      
    }
    BClist[[i]] <- list
    names(BClist)[i] <- list$id
  }
  return(BClist)
}


#add_missing_bc####
#this function will take a list of BC_list objects and will ensure that all BC_list objects have the same bc_table.
#bc_table which does not belong to a particular BC_list will have a read count of 0.
add_missing_bc <- function(BC_list, verbose = TRUE){
  library(dplyr)
  #checking input
  if(length(BC_list) < 2){
    stop("a list of BC objects should be provided, not a single BC object.")
  }
  #first it will create a vector containing all bc_table present across all BC objects in the list
  all_bc_table <- as.character(unique(bind_rows(BC_list)$bc_table$barcode))
  if(verbose){
    print(paste("found a total of", length(all_bc_table), "bc_table"))
  }
  
  #now for each BC object in the list, it will check which bc_table are missing and add them.
  for(i in c(1:length(BC_list))){
    data <- BC_list[[i]]$bc_table #getting the BC object
    data$bc_added <- FALSE #creates a column indicating which bc_table were added. As this is the original data, no barcode is added so far, so FALSE
    data$barcode <- as.character(data$barcode) #removing factor levels
    missing_bc_table <- as.character(data$barcode) #making a vector containg all bc_table found in the BClist
    missing_bc_table <- all_bc_table[which(!all_bc_table %in% missing_bc_table)] #check which bc_table in the all_bc_table vector are missing in this BC object.
    #if no barcode is missing, skip to the next BC object
    if(length(missing_bc_table) < 1){ 
      if(verbose){
        print(paste("no missing barcode present in BC object number ", i, ", skipping it...", sep = ""))
      }
      next
    }
    if(verbose){
      print(paste(length(missing_bc_table), "bc_table with read count 0 were added to the BC object number", i))
    }
    missing_data <- as.data.frame(matrix(nrow = length(missing_bc_table), ncol = ncol(data)))
    colnames(missing_data) <- colnames(data)
    missing_data$barcode <- missing_bc_table
    missing_data$count <- 0
    missing_data$percent <- 0
    missing_data$position <- NA
    missing_data$bc_added <- TRUE
    if("bc_length" %in% colnames(missing_data)){
      missing_data <- mutate(missing_data, bc_length = nchar(barcode))
    }
    data <- rbind(data, missing_data)
    data <- arrange(data, desc(percent))
    data$position <- c(1:nrow(data))
    BC_list[[i]]$bc_table <- data #change the reads attribute of the BC object to the new data.frame containing the missing bc_table.
  }
  return(BC_list)
}

#bind_bc_lists####
#this function will take a list of BC_list objects and will merge them in a tidy dataframe.
#this dataframe is directly compatible with ggplot
bind_bc_lists <- function(x){
  library(tidyverse)
  x <- bind_rows(x)
  x <- cbind(x$bc_table, x[,-1])
  return(x)
}


#combine_bc_lists####
#this function will take a list of BClist objects and will combine them into a single BClist object
#different approaches are available for combination: union, intersection and mean.
#by default the combined object will receive the other information such as timepoints, replicates, groups, etc, from the 1st object
combine_bc_lists <- function(x, method = "union"){
    library(tidyverse)
    method <- tolower(method)
    l_length <- length(x)
    l_names <- names(x)
    a <- bind_rows(x)$bc_table
    if(method == "intersection"){
      print(paste("combining", l_length, "objects using method", method))
      a <- a %>%
        group_by(barcode) %>%
        summarise(nsamples = sum(count > 0), count = sum(count)) %>%
        filter(nsamples == max(nsamples)) %>%
        select(barcode, count)
    }else{
      if(method == "union"){
        print(paste("combining", l_length, "objects using method", method))
        a <- a %>%
          group_by(barcode) %>%
          summarise(count = sum(count)) 
      }else{
        if(method %in% c("mean", "average")){
          print(paste("combining", l_length, "objects using method", method))
          a <- a %>%
            group_by(barcode) %>%
            summarise(count = mean(count, na.rm = TRUE)) 
        }else{
        stop("method should be \"union\", \"intersection\", or \"average\"")
      }
      }
    }
    a <- a %>%
      mutate(percent = count/sum(count)) %>%
      arrange(desc(percent))
    a$position <- c(1:nrow(a))
    x <- x[[1]]
    x$bc_table <- a
    x$id <- paste0(x$id, "_combined")
    return(x)
}  


#compare_bc_lists####
compare_bc_lists <- function(to_compare, 
                             reference, 
                             na_to_zero = TRUE,
                             combine_method = "mean",
                             cname_prefix = "ini_"){
  library(tidyverse)
  #first it will chech if comparisons involve single BClist objects or lists of objects
  to_compare_length <- length(to_compare)
  print(paste("comparing", to_compare_length, "BClist objects to reference."))
  
  if(to_compare_length == 1){
    to_compare <- to_compare[[1]]
  }
  if("bc_table" %in% names(to_compare)){
    to_compare_length <- 1
  }
  
  reference_length <- length(reference)
  if(reference_length == 1){
    reference <- reference[[1]]
  }
  if("bc_table" %in% names(reference)){
    reference_length <- 1
  }
  if(reference_length > 1){
    print(paste("reference has multiple BClist objects. Combining them into a single object with method:", combine_method))
    reference <- combine_bc_lists(reference, method = combine_method)
  }
  
  #now it will create a function to compare 2 BClist objects
  compare <- function(x, y){
    cnames <- colnames(x$bc_table)
    indexes <- match(x$bc_table$barcode, y$bc_table$barcode)
    x$bc_table$prev_count <- y$bc_table$count[indexes]
    x$bc_table$prev_percent <- y$bc_table$percent[indexes]
    x$bc_table$prev_position <- y$bc_table$position[indexes]
    if(na_to_zero){
      if(sum(is.na(x$bc_table > 0))){
        x$bc_table <- x$bc_table %>%
          mutate(prev_percent = ifelse(is.na(prev_percent), 0, prev_percent)) %>%
          mutate(prev_count = ifelse(is.na(prev_count), 0, prev_count))
      }
    }
    x$bc_table$prev_fold_change <- log(x$bc_table$percent/x$bc_table$prev_percent)
    indexes <- which(!colnames(x$bc_table) %in% cnames)
    cnames <- colnames(x$bc_table)[indexes]
    cnames <- gsub("prev_", cname_prefix, cnames)
    colnames(x$bc_table)[indexes] <- cnames
    rm(cnames)
    rm(indexes)
    x$comparison <- paste0(x$id, "/", y$id)
    return(x)
  }
  
  if(to_compare_length > 1){
    to_compare <- lapply(to_compare, compare, reference)
  }else{
    to_compare <- compare(to_compare, reference)
    
  }
  return(to_compare)
}

#calculate_HD####
#calculates the minimum Hamming distance between a barcode and any other barcode in a vector.
calculate_HD <- function(x, verbose = FALSE){
  ids <- x
  x <- as.character(x)
  x <- strsplit(x, split = "")
  diff_matrix <- matrix(nrow = length(x), ncol = length(x))
  colnames(diff_matrix) <- ids
  rownames(diff_matrix) <- ids
  for(i in c(1:length(x))){
    barcode_1 <- x[[i]]
    barcode_1_size <- length(barcode_1)
    vec <- vector(length = length(x))
    for(j in c(1:length(x))){
      barcode_2 <- x[[j]]
      barcode_2_size <- length(barcode_2)
      if(barcode_1_size != barcode_2_size){
        stop(paste("barcode", i, "and", j, "have a different number of nucleotides"))
      }
      if(verbose == TRUE){
        print(paste("comparing barcode", i, "with", j))
      }
      vec[j] <- sum(barcode_1 != barcode_2)
    }
    diff_matrix[,i] <- vec
  }
  return(diff_matrix)
}

#minHD_network####
#this function exports a network graph representing the similarities between bc_table in a sample
minHD_network <- function(bclist, sample_name, HD_threshold = 3, min_count = 2, bc_lim = NA, output_path, mass_factor = 0.2, solver = 4){
  library(visNetwork)
  library(igraph)
  
  #the visnetwork package has 4 options for solving the physics of the network. They are listed here.
  vis_solver <- c("barnesHut", "repulsion", "hierarchicalRepulsion", "forceAtlas2Based")
  
  #making the color palette for the graph
  color_vec <- c("#ff0000", "#ff9100", "#fffb00")
  color_vec <- colorRampPalette(color_vec)
  
  #getting the sample name 
  if(missing(sample_name)){
    sample_name <- bclist$id
  }
  #to calculate hamming distances bc_table must have the same length. Thus, aberrant bc_table are removed
  bclist <- bclist$bc_table
  bclist <- bclist %>%
    mutate(bc_length = nchar(as.character(barcode)))%>%
    filter(bc_length == median(bc_length)) %>%
    arrange(desc(count))
  #in case bc_lim is set, will keep just the number of bc_table specified by bc_lim.
  if(is.numeric(bc_lim)){
    if(bc_lim < nrow(bclist)){
      bclist <- bclist[1:bc_lim,]
    }
  }
  #removing bc_table with count lower than min_count
  bclist <- filter(bclist, count >= min_count)
  #selecting the bc_table
  bc_table <- as.character(bclist$barcode)
  #calculating nucleotide distances between bc_table
  dist_matrix <- calculate_HD(bc_table)
  diag(dist_matrix) <- NA
  if(!is.numeric(HD_threshold)){
    HD_threshold <- max(dist_matrix, na.rm = TRUE)
  }
  dist_matrix[dist_matrix > HD_threshold] <- NA
  #creating a minimum spaning tree to represent these distances
  network <- graph_from_adjacency_matrix(dist_matrix, weighted = TRUE, mode = "undirect")
  network <- minimum.spanning.tree(network)
  edges <- toVisNetworkData(network)$edges
  colnames(edges)[3] <- "steps"
  #creating nodes and edge lists
  nodes <- data.frame(id = bclist$barcode, label = bclist$position, title = bclist$barcode, size = sqrt(bclist$count))
  #removing edges with distance greater than the HD_threshold
  edges <- filter(edges, steps <= HD_threshold)
  #creating the color vector for the steps in the network
  color_vec <- color_vec(HD_threshold)
  
  #formating nodelist
  #first it will create a sub function to format the labels in the nodes so numbers do not affect the size of the nodes.
  #this is done by adding space characters in numbers with less characters.
  format <- function(x){
    x <- as.character(x)
    max_size <- max(nchar(x), na.rm = TRUE)
    x <- data.frame(id = x)
    x <- x %>%
      mutate(size_diff = max_size-nchar(id)) %>%
      mutate(size_diff = size_diff *1.01)%>%
      mutate(size_diff = round(size_diff/2))
    x$label <- NA
    for(i in c(1:nrow(x))){
      id <- x$id[i]
      size_diff <- x$size_diff[i]
      if(is.na(size_diff)){
        next
      }
      spaces <- paste0(rep(" ", times = size_diff), collapse = "")
      x$label[i] <- paste0(spaces, id, spaces)
    }
    return(as.character(x$label))
  }
  
  nodes$label <- format(nodes$label)
  nodes <- nodes %>%
    mutate(mass = (size/min(size))*mass_factor) %>%
    mutate(mass = ifelse(mass < 1, 1, mass))%>%
    mutate(font.size = round(size/2))
  #formating edgelist
  edges <- edges %>%
    mutate(label = ifelse(steps == 1, NA, as.character(steps))) %>%
    mutate(color = color_vec[steps])
  #generating and exporting the network
  if(missing(output_path)){
    output_path <- getwd()
  }
  file_name <- paste0(sample_name, "_Network.html")
  visNetwork(nodes, edges, main = paste0(sample_name, ": HD threshold = ", HD_threshold, " and maximum bc_table limitted to ", bc_lim), height = "95vh", width = "100%") %>%
    visIgraphLayout(physics = TRUE, layout = "layout_with_dh") %>%
    visPhysics(solver = vis_solver[solver], stabilization = FALSE) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE , collapse = FALSE) %>%
    visEdges(width = 10) %>%
    visNodes(color = list(border = "black", background = "#383838"),
             font = list(face = "consolas", color = "white"), 
             shadow = TRUE, 
             shape = "circle", 
             scaling = list(label = list(enabled = T, drawThreshold = 1))) %>%
    visExport() %>% 
    visSave(paste0(output_path, "/", file_name), selfcontained = FALSE)
  #file.copy(file_name, paste0(output_path, "/", file_name))
  #file.remove(file_name)
}

#seqlogo_plot####
#this function generates a sequence logo representing the nucleotide diversity of each library
seqlogo_plot <- function(x, save_plot = FALSE, output_path){
  library(tidyverse)
  library(ggplot2)
  library(ggseqlogo)
  sample_name <- as.character(x$id)
  x <- x$bc_table
  if(missing(output_path)){
    output_path <- getwd()
  }
  x$barcode <- as.character(x$barcode)
  n_bc_table_ini <- nrow(x)
  x <- x %>%
    mutate(bc_length = nchar(barcode)) %>%
    filter(bc_length == median(bc_length))
  n_bc_table_fin <- nrow(x)
  diff <- n_bc_table_ini - n_bc_table_fin
  x <- x$barcode
  p <- ggseqlogo(data = x, method = "prob")
  p <- p + labs(title = paste("Sequence logo for sample:", sample_name), 
                subtitle = paste("Number of sequences:", n_bc_table_ini, "of which", diff, "where removed due to wrong sequence size."), 
                x = "Nucleotide Position")
  if(save_plot){
    ggsave(paste0(output_path,"/",sample_name, "seqlogo.pdf"), width = 9, height = 3, units = "in")
  }else{
    return(p)
  }
  
}

#bc_length_plot####
bc_length_plot <- function(x, save_plot = FALSE, output_path){
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  
  sample_name <- as.character(x$id)
  x <- x$bc_table
  if(missing(output_path)){
    output_path <- getwd()
  }
  
  if(!"bc_length" %in% colnames(x)){
    x <- mutate(x, bc_length = nchar(barcode))
  }
  
  x <- x %>%
    group_by(bc_length) %>%
    summarize(count = sum(count), number_of_bc_table = n())
  scale <- c(0, max(x$count, na.rm = TRUE))
  scale[2] <- ceiling(scale[2]/10)*10
  coeff <- scale[2]/max(x$number_of_bc_table, na.rm = T)
  p <- x %>%
    ggplot(aes(x = factor(bc_length), y = count))+
    geom_col()+
    geom_point(aes(y = number_of_bc_table*coeff), color = "blue")+
    geom_line(aes(y = number_of_bc_table*coeff), color = "blue", group = 1)+
    scale_y_continuous(sec.axis = sec_axis(~./coeff, name="Number of different bc_table"))+
    labs(title = sample_name, x = "Barcode length", y = "number of reads")
  if(save_plot){
    ggsave(paste0(output_path,"/",sample_name, "bc_length.pdf"), width = 9, height = 3, units = "in")
  }else{
    return(p)
  }
}

#bc_dist_plot####
bc_dist_plot <- function(x, top, save_plot = FALSE, output_path){
  library(tidyverse)
  library(ggplot2)
  sample_name <- as.character(x$id)
  x <- x$bc_table
  x$cummulative <- cumsum(x$percent)
  
  if(missing(top)){
    top <- max(x$position)
  }
  if(missing(output_path)){
    output_path <- getwd()
  }
  
  scale <- c(0, max(x$percent, na.rm = TRUE))
  scale[2] <- ceiling(scale[2]*100)/100
  coeff <- 1/scale[2]
  p <- x %>%
    filter(position <= top) %>%
    mutate(cummulative = cummulative/coeff)%>%
    ggplot(aes(x = position, y = percent))+
    geom_col()+
    geom_line(color = "blue", aes(y = cummulative))+
    scale_x_continuous()+
    scale_y_continuous(limits = scale,
                       sec.axis = sec_axis(~.*coeff, name="Cummulative"))+
    labs(title = sample_name, y = "Frequency", x = "Barcode")
  if(save_plot){
    ggsave(paste0(output_path,"/",sample_name, "bc_distribution.pdf"), width = 9, height = 3, units = "in")
  }else{
    return(p)
  }
}

#similar_bc_calc####
#this function will count the number of similar bc_table in a given sample: not working yet, do not use!
similar_bc_calc <- function(x){
  if(!"bc_length" %in% colnames(x)){
    x <- x %>%
      mutate(bc_length = nchar(barcode))
  }
  #create a vector with all bc_table
  bcs <- as.character(x$barcode)
  #find the normal size of the bc_table
  norm_len <- median(x$bc_length)
  #determine how many nucleotides will be removed from the extremeties
  to_cut <- round(norm_len/10)
  #remove the nucleotides from the extremeties
  bcs <- substr(bcs, to_cut, norm_len-to_cut)
  #create an extra column where the number of simiilar bc_table will be stored
  x$n_similar_bc_table <- NA
  #for each barcode, count how many similar bc_table are found
  for(i in c(1:nrow(x))){
    bc <- bcs[i]
    x$n_similar_bc_table[i] <- sum(grepl(bc, bc_char$barcode))
  }
  return(x)
}

