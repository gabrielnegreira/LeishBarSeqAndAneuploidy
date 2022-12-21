#Set of custom functions to work with BGlist objects.
#BGlist objects are lists consisting of a bedgraph file in dataframe format (a BG_table) and additional elements with meta data information
#it must contain at least the folowing columns: chromosome, start, end, read_count

#import_bg_files####
#this function will import bedgraph tables into a single list of BGlist objects.

####inporting files####
import_bg_lists <- function(inputs_path, pattern, input_type, column_names, sample_info){
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
  files <- list.files(path = inputs_path,
                      pattern = pattern)
  
  #to do, change sample_info format to accept excel files instead of csvs
  sample_info <- read_excel(sample_info)
  
  print(paste("found", length(files), "files"))
  
  BGlist <- vector(mode = "list", length = length(files))
  for(i in c(1:length(files))){
    bg_table <- read.delim(paste(inputs_path, files[i], sep = "/"), header = TRUE)
    bg_table <- data.frame(bg_table)
    colnames(bg_table) <- column_names
    bg_table$bin <- c(1:nrow(bg_table))
    bg_table <- bg_table %>%
      group_by(chromosome) %>%
      mutate(chromo_bin = abs(max(bin) - bin - (max(bin) - min(bin)))+1) %>% #there's probably an easier way of doing this but my head is hurting already
      ungroup()
    list <- list(bg_table = bg_table, file = files[i])
    #appending information from the sample_info table
    row <- inv_grep(sample_info$id, files[i])
    for(j in c(1:ncol(sample_info))){
      variable <- colnames(sample_info)[j]
      value <- as.character(sample_info[row,j])
      list <- append(list, value)
      names(list)[length(names(list))] <- variable
    }
    BGlist[[i]] <- list
    names(BGlist)[i] <- list$id
  }
  return(BGlist)
}


#remove_outliers####
#remove_outlier_bins will remove outlier bins based on read count (in a per-chromosome base) and/or gc content (whole genome)
remove_outlier_bins <- function(x, gc_outlier_method = "boxplot", read_count_outlier_method = "boxplot", save_plot = FALSE, plot_path){
  
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  
  if(missing(plot_path)){
    plot_path <- getwd()
  }
  
  if(is.na(gc_outlier_method) | tolower(gc_outlier_method) %in% c("none", "na", "null")){
    gc_outlier_method <- "none"
  }
  
  if(is.na(read_count_outlier_method) | tolower(read_count_outlier_method) %in% c("none", "na", "null")){
    read_count_outlier_method <- "none"
  }
  
  bg_table <- x$bg_table
  #ensuring all plots will be done in the same scale...
  
  if(save_plot){
    ylimits <- c(0, max(bg_table$read_count)*1.1)
    #getting extra information 
    bins_position <- bg_table %>%
      group_by(chromosome) %>%
      summarize(first_bin = min(bin), last_bin = max(bin), medium_bin = median(bin))
    to_plot <- bg_table%>%
      group_by(chromosome) %>%
      mutate(median_read_count = median(read_count))
    
    p1 <- ggplot(to_plot, aes(x = bin, y = read_count, col = chromosome, group = chromosome))+
      geom_point(size = 0.75, alpha = 0.5)+
      geom_line(aes(y = median_read_count), color = "red", size = 1)+
      guides(col = "none")+
      labs(title = "Before removal of outlier bins", x = "Chromosome", y = "Read Count")+
      scale_x_continuous(label = bins_position$chromosome, breaks = bins_position$medium_bin)+
      scale_y_continuous(limits = ylimits)+
      scale_color_manual(values = rep(c("black", "orange"), times = 36))+  
      theme_linedraw()+
      theme(title = element_text(size = 20, face = "bold"), 
            axis.title = element_text(size = 20, face = "bold"), 
            axis.text = element_text(size = 19, angle = 45, hjust = 1), 
            panel.grid.major.x = element_blank(), 
            panel.grid.major = element_line(color = "grey"), 
            panel.grid.minor = element_blank())
  }
  
  
  #now will remove bins with outlier gc content
  if(gc_outlier_method != "none"){
    gc_out <- TRUE
    if(length(gc_outlier_method) == 2){
      qmin <- quantile(bg_table$gc, gc_outlier_method)
      qmax <- qmin[length(qmin)]
      qmin <- qmin[1]
      print(paste("removing bins with gc content lower than ", qmin, " or higher than ", qmax, "...", sep = ""))
      before <- nrow(bg_table)
      bg_table <- filter(bg_table, gc >= qmin & gc <= qmax)
      after <- nrow(bg_table)
      print(paste(before-after, "bins removed"))
    }else{
      if(tolower(gc_outlier_method) == "boxplot"){
        print(paste("removing bins with outlier gc content using the boxplot.stats() function..."))
        bg_table <- filter(bg_table, ! gc %in% boxplot.stats(bg_table$gc)$out)
      }
    }
  }
  
  #removing bins with outlier read count
  if(read_count_outlier_method != "none"){
    rc_out <- TRUE
    if(length(read_count_outlier_method) == 2){
      print("removing bins with oulier read counts using quantiles...")
      before <- nrow(bg_table)
      bg_table <- bg_table %>%
        group_by(chromosome) %>%
        mutate(keep = ifelse(read_count >= quantile(read_count, read_count_outlier_method)[1] & read_count <= quantile(read_count, read_count_outlier_method)[2], TRUE, FALSE))%>%
        filter(keep) %>%
        mutate(keep = NULL)
      after <- nrow(bg_table)
      print(paste(before-after, "bins removed"))
    }else{
      if(tolower(read_count_outlier_method) == "boxplot"){
        print("removing bins with oulier read counts using the boxplot.stats() function...")
        before <- nrow(bg_table)
        bg_table <- bg_table %>%
          group_by(chromosome) %>%
          mutate(keep = ifelse(!read_count %in% boxplot.stats(read_count)$out, TRUE, FALSE))%>%
          filter(keep) %>%
          mutate(keep = NULL)
        after <- nrow(bg_table)
        print(paste(before-after, "bins removed"))
      }
    }
  }
  
  if(save_plot){
    plot_title <- ifelse(rc_out & gc_out, "After removal of outlier bins based on gc and read count",
                         ifelse(rc_out & !gc_out, "After removal of outlier bins based on read count only",
                                ifelse(!rc_out & gc_out, "After removal of outlier bins based on gc only",
                                       "no removal of outlier bins")))
    if(rc_out | gc_out){
      plot_title <- paste0(plot_title, ": (", before-after, " bins removed)")
    }
    #saving plot after removal of outlier bins
    bins_position <- bg_table %>%
      group_by(chromosome) %>%
      summarize(first_bin = min(bin), last_bin = max(bin), medium_bin = median(bin))
    to_plot <- bg_table%>%
      group_by(chromosome) %>%
      mutate(median_read_count = median(read_count))
    
    p2 <- ggplot(to_plot, aes(x = bin, y = read_count, col = chromosome, group = chromosome))+
      geom_point(size = 0.75, alpha = 0.5)+
      geom_line(aes(y = median_read_count), color = "red", size = 1)+
      guides(col = "none")+
      labs(title = plot_title, x = "Chromosome", y = "Read Count")+
      scale_x_continuous(label = bins_position$chromosome, breaks = bins_position$medium_bin)+
      scale_y_continuous(limits = ylimits)+
      scale_color_manual(values = rep(c("black", "orange"), times = 36))+  theme_linedraw()+
      theme(title = element_text(size = 20, face = "bold"), 
            axis.title = element_text(size = 20, face = "bold"), 
            axis.text = element_text(size = 19, angle = 45, hjust = 1), 
            panel.grid.major.x = element_blank(), 
            panel.grid.major = element_line(color = "grey"), 
            panel.grid.minor = element_blank())
    p <- plot_grid(p1,p2, ncol = 1, align = "v")
    
    sample_name <- ifelse("name" %in% names(x), x$name, x$id)
    
    ggsave(paste0(paste(plot_path, sample_name, sep = "/"), "_outlier_removal.png"),
           plot = p, width = 15, height = 15, units = "in", dpi = 300, device = "png")
  }
  x$bg_table <- bg_table
  return(x)
}

correct_gc_bias <- function(x, save_plot = FALSE, plot_path, baseline_ploidy = 2){
  #defaults
  if(missing(plot_path)){
    plot_path <- getwd()
  }
  
  #subfunctions
  Mode <- function(x) {
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    ux[tab == max(tab)]
  }
  
  #getting the bg_table
  bg_table <- x$bg_table
  #finding the baseline chromosomes based on rough estimation of somies
  rough_somy <- bg_table %>%
    group_by(chromosome) %>%
    summarise(median = median(read_count, na.rm = TRUE)) %>% 
    mutate(normalized = median/median(median)) %>% 
    mutate(scaled = normalized * baseline_ploidy) %>% 
    mutate(rounded = round(scaled)) 
  baseline <- Mode(rough_somy$rounded)
  bg_table$rough_somy <- rough_somy$rounded[match(bg_table$chromosome, rough_somy$chromosome)]
  
  baseline <- rough_somy %>%
    filter(rounded == baseline) %>%
    select(chromosome)
  
  baseline <- c(baseline$chromosome)
  print(paste("correcting for gc content bias based on chromosomes", paste(baseline, collapse = " "), "...", collapse = ""))
  
  #subseting the baseline chromosomes
  gc_model_bins <- bg_table %>%
    filter(chromosome %in% baseline) %>%
    filter(!gc %in% boxplot.stats(gc_model$gc)$out) %>%
    filter(!read_count %in% boxplot.stats(gc_model$read_count)$out)
  
  #building a loess model
  gc_model <- loess(read_count ~ gc, 
                    data = gc_model_bins, 
                    span = 0.25,
                    control = loess.control(surface = "direct"))
  
  gc_model_bins$fitted <- gc_model$fitted
  correction_factor <- predict(gc_model, bg_table$gc)
  correction_factor <- correction_factor/mean(correction_factor)
  correction_factor <- 1/correction_factor
  
  #correcting read counts
  bg_table$correction_factor <- correction_factor
  bg_table <- bg_table %>%
    mutate(read_count_corrected = read_count * correction_factor) %>%
    mutate(correction_difference = read_count - read_count_corrected)

  # making plot
  if(save_plot){
    bins_position <- bg_table %>%
      group_by(chromosome) %>%
      summarize(first_bin = min(bin), last_bin = max(bin), medium_bin = median(bin))
    to_plot <- bg_table%>%
      group_by(chromosome) %>%
      mutate(median_read_count = median(read_count_corrected))
    
    ylimits <- c(0, max(c(to_plot$read_count, to_plot$read_count_corrected))*1.1)
    p1 <- ggplot(to_plot, aes(x = bin, y = read_count, col = chromosome, group = chromosome))+
      geom_point(size = 0.75, alpha = 0.5)+
      geom_line(aes(y = median_read_count), color = "red", size = 1)+
      guides(col = "none")+
      labs(title = "Before GC bias correction", x = "Chromosome", y = "Read Count")+
      scale_x_continuous(label = bins_position$chromosome, breaks = bins_position$medium_bin)+
      scale_y_continuous(limits = ylimits)+
      scale_color_manual(values = rep(c("black", "orange"), times = 36))+  
      theme_linedraw()+
      theme(title = element_text(size = 20, face = "bold"), 
            axis.title = element_text(size = 20, face = "bold"), 
            axis.text = element_text(size = 19, angle = 45, hjust = 1), 
            panel.grid.major.x = element_blank(), 
            panel.grid.major = element_line(color = "grey"))
    
    
    p2 <- ggplot(to_plot, aes(x = bin, y = read_count_corrected, col = chromosome, group = chromosome))+
      geom_point(size = 0.75, alpha = 0.5)+
      geom_line(aes(y = median_read_count), color = "red", size = 1)+
      guides(col = "none")+
      labs(title = "After GC bias correction", x = "Chromosome", y = "Read Count")+
      scale_x_continuous(label = bins_position$chromosome, breaks = bins_position$medium_bin)+
      scale_y_continuous(limits = ylimits)+
      scale_color_manual(values = rep(c("black", "orange"), times = 36))+  
      theme_linedraw()+
      theme(title = element_text(size = 20, face = "bold"), 
            axis.title = element_text(size = 20, face = "bold"), 
            axis.text = element_text(size = 19, angle = 45, hjust = 1), 
            panel.grid.major.x = element_blank(), 
            panel.grid.major = element_line(color = "grey"))
    
    p1 <- plot_grid(p1,p2, ncol = 1, align = "v")
    
    scale <- log(to_plot$correction_factor, 2)
    scale <- max(abs(c(min(scale, na.rm = T), max(scale, na.rm = T))))
    scale <- scale *1.1
    scale <- c(-scale, scale)
    p4 <- to_plot %>%
      ggplot(aes(x = gc, y = log(correction_factor, 2)))+
      geom_line(group = 1, show.legend = FALSE)+
      scale_y_continuous(limits = scale)+
      labs(title = "correction factor vs GC content", x = "GC percentage", y = "Correction factor (log2)")+
      theme(text = element_text(size = 18))
    
    p5 <- to_plot %>%
      ggplot(aes(x = gc, y = read_count, color = factor(rough_somy)))+
      geom_point(shape = 1)+
      geom_smooth()+
      labs(title = "Norm. depth Vs GC content - before correction", x = "GC percentage", y = "read count", color = "Rough\nsomy")+
      scale_color_brewer(palette = "RdYlBu", direction = -1)+
      theme(text = element_text(size = 18))
    
    p6 <- to_plot %>%
      ggplot(aes(x = gc, y = read_count_corrected, color = factor(rough_somy)))+
      geom_point(shape = 1)+
      geom_smooth()+
      labs(title = "Norm. depth Vs GC content - after correction", x = "GC percentage", y = "corrected read count", color = "Rough\nsomy")+
      scale_color_brewer(palette = "RdYlBu", direction = -1)+
      theme(text = element_text(size = 18))
    
    to_plot <- gc_model_bins
    p3 <- to_plot %>%
      ggplot(aes(x = gc, y = read_count))+
      geom_point(shape = 1)+
      geom_line(aes(y = fitted), col = "red", size = 2)+
      labs(title = "loess model of GC content vs read count", x = "GC percentage", y = "read count", color = "Baseline\nChromosome")+
      theme(text = element_text(size = 18))
    
    p3 <- plot_grid(p3, p4, p5, p6, ncol =2, align = "v", axis = "tblr")
    p <- plot_grid(p1,p3, nrow = 1, rel_widths = c(0.4,0.6))  
    #saving the plot
    sample_name <- ifelse("name" %in% names(x), x$name, x$id)
    ggsave(paste0(paste(plot_path, sample_name, sep = "/"), "_gc_correction.png"),
           plot = p, width = 30, height = 15, units = "in", dpi = 300, device = "png")
  }
  bg_table$rough_somy <- NULL
  bg_table$read_count <- bg_table$read_count_corrected
  bg_table$read_count_corrected <- NULL
  x$bg_table <- bg_table
  return(x)
}

#calc_somy####
#calc_somy will take a BGlist obbject and will calculate the somy

calc_somy <- function(x, baseline_ploidy = 2){
  library(tidyverse)
  bg_table <- x$bg_table

  somy_table <- bg_table %>%
    group_by(chromosome) %>%
    summarise(min_bin = min(bin), 
              max_bin = max(bin), 
              min_chromo_bin = min(chromo_bin),
              max_chromo_bin = max(chromo_bin),
              median_read_count = median(read_count)) %>%
              mutate(normalized_depth = median_read_count/median(median_read_count)) %>% 
              mutate(somy = normalized_depth * baseline_ploidy) %>% 
              mutate(integer_somy = round(somy)) 
  x <- c(list(somy_table), x)
  names(x)[1] <- "somy_table"
 return(x) 
}

#bind_bg_lists####
bind_bg_lists <- function(x, table_to_bind = "bg_table"){
  library(tidyverse)
  if(!table_to_bind %in% c("bg_table", "somy_table")){
    stop("table_to_bind must be either \"bg_table\" or \"somy_table\"")
  }
  
  if(table_to_bind == "bg_table"){
    x <- lapply(x, function(y){
      if("somy_table" %in% names(y)){
        warning("somy_table will be removed.")
        y$somy_table <-NULL
      }
      return(y)
    })
  }else{
    if(table_to_bind == "somy_table"){
      x <- lapply(x, function(y){
        if("bg_table" %in% names(y)){
          warning("bg_table will be removed.")
          y$bg_table <-NULL
        }
        return(y)
      })
    }
  }
  x <- bind_rows(x)
  names(x)[1] <- "first"
  x <- cbind(x$first, x[,-1])
  return(x)
}
