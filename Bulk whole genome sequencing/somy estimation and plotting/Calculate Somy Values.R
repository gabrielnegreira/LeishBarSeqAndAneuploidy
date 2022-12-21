#Author: Gabriel H. Negreira

#code to estimate somy values from binned bedgraph files

#inputs####
inputs_path <- "inputs/PAT_FS"
outputs_path <- "outputs/PAT_FS"
sample_info <- "inputs/PAT_FS/sample_info.xlsx"
file_id <- ".binned_depth_gc.csv" #string which is present in all files which should be analyzed. 
column_names <-  c("chromosome", "start", "end", "read_count", "gc") #name of the columns of the datasets. 
save_to_disk <- TRUE #if true will save the BGlist object to disk.

#parameters####
gc_correction <- TRUE # will attempt to correct gc bias. Check gc plot to see if it worked.
remove_outliers <- FALSE #if true will remove bins with outlier read count.
gc_outlier_method <- "boxplot" #removes bins with outlier gc content. Use "none", "boxplot" to use the outliers determined by the boxplot.stats() function or provide a vector with low and high quantiles such as "c(0.01, 0.99)".
read_count_outlier_method <- "boxplot" #same as above, but removes outlier read count values in a per-chromosome manner. 
baseline_ploidy <- 2 #for triploid samples set this to 3. 
color_option <- 1 #two color pallets available. 
smooth <- FALSE #if true will use loess to smooth the depth of bins in each chromosome before calculating somies. Might be useful for noisy data or data from ATACseq.
smooth_span <- 0.1 # only used if smooth is TRUE.
plot_resolution <- 72 #resolution for the plots
#libraries####
library(readr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(hues)
###################start#########################

#setting color pallets#### 
if(color_option == 1){
  integer_col <- c("#001221","#002342", "#014175","#00c3ff", "#33ff00", "#fffa00", "#ffa600", "#D73027", "#A50026", "#541b1b", "#4d0600")
  heat_col <- c("#001221", "#002342", "#002342", "#014175", "#035ba3", "#00c3ff", "#00ffee", "#33ff00", "#ccff00", "#fffa00","#ffa600", "#D73027", "#A50026", "#541b1b", "#4d0600")
}else{
  integer_col <- c("#001221", "#4575B4", "#a4cbe0", "#E0F3F8", "#FFFFBF", "#FDAE61", "#F46D43", "#D73027", "#A50026", "#541b1b")
  heat_col <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026", "#541b1b")
}
color_palette <- colorRampPalette(heat_col)



#creating output folder
gc_correction_plots <- paste0(outputs_path,"/GC_correction")
bin_removal_plots <- paste0(outputs_path,"/outlier_bins_removal")
dir.create(outputs_path)
dir.create(bin_removal_plots)
dir.create(gc_correction_plots)
#inporting files
BGlist <- import_bg_lists(inputs_path = inputs_path,
                          pattern = file_id,
                          column_names = column_names,
                          sample_info = sample_info)
#removing kDNA
BGlist <- lapply(BGlist, function(x){
  x$bg_table <- filter(x$bg_table, chromosome != "Ld37")
  return(x)
})

#removing outlier bins and correcting for gc bias.
if(remove_outliers){
  BGlist <- BGlist %>%
    lapply(remove_outlier_bins,
           gc_outlier_method = gc_outlier_method,
           read_count_outlier_method = read_count_outlier_method,
           save_plot = TRUE,
           plot_path = bin_removal_plots) 
}

#correcting GC bias
if(gc_correction){
  BGlist <- BGlist %>%
    lapply(correct_gc_bias,
           save_plot = TRUE,
           plot_path = gc_correction_plots)
}

#calculating somies. 
BGlist <- BGlist %>%
  lapply(calc_somy)

#saving BGlist to disk
if(save_to_disk){
  write_rds(BGlist, paste0(outputs_path,"/BGlist.rds"))
}

#saving a summary of somy values for each sample
BGlist %>%
  bind_bg_lists(table_to_bind = "somy_table") %>%
  select(name, chromosome, somy) %>%
  pivot_wider(names_from = chromosome, values_from = somy)%>%
  write.csv(file = paste0(outputs_path, "/somies.csv"))

#summary plot
BGlist %>%
  bind_bg_lists(table_to_bind = "somy_table") %>%
  #mutate(sample_group = ifelse(sample_group == "control", "PBS", "382 µM SbIII"))%>%
  #mutate(chromosome = factor(chromosome))%>%
  mutate(chromosome = factor(chromosome, levels = rev(levels(factor(chromosome)))))%>%
ggplot(aes(x = name, y = chromosome, fill = somy, label = round(somy, 2)))+
  geom_tile()+
  #geom_text(color = "white", size = 2)+
  scale_fill_gradientn(colors = heat_col, limits = c(0,8))+
  #scale_y_discrete(breaks = c("Ld01", "Ld06", "Ld11", "Ld16", "Ld21", "Ld26", "Ld31", "Ld36"))+
  facet_grid(cols = vars(sample_group), space = "free", scales = "free")+
  labs(y = "Chromosome", x = "", fill = "Somy")+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        panel.grid = element_blank())

#MIL_FS plot
BGlist %>%
  bind_bg_lists(table_to_bind = "somy_table") %>%
  mutate(chromosome = gsub("Ld", "", chromosome))%>%
  mutate(chromosome = factor(chromosome, levels = rev(levels(factor(chromosome)))))%>%
  mutate(timepoint = paste("passage", timepoint))%>%
  mutate(sample_group = gsub("MIL_C", "McPOP", sample_group))%>%
  mutate(sample_group = gsub("MIL_P1_", "", sample_group))%>%
  mutate(sample_group = gsub("MIL_P9_", "", sample_group))%>%
  mutate(concentration = ifelse(grepl("McPOP", sample_group), NA, sample_group))%>%
  mutate(sample_group = ifelse(grepl("McPOP", sample_group), "PBS (control)", paste(sample_group, "µM")))%>%
  mutate(sample_group = factor(sample_group, levels = c("PBS (control)", "25 µM", "50 µM", "100 µM")))%>%
  mutate(population_name = ifelse(grepl("PBS", sample_group), "McPOP", paste0("Me", concentration, "POP")))%>%
  mutate(population_name = paste0(population_name, replicate))%>%
  ggplot(aes(x = population_name, y = chromosome, fill = somy, label = round(somy, 2)))+
  geom_tile()+
  #geom_text(color = "white", size = 2)+
  scale_fill_gradientn(colors = heat_col, limits = c(0,8))+
  #scale_y_discrete(breaks = c("Ld01", "Ld06", "Ld11", "Ld16", "Ld21", "Ld26", "Ld31", "Ld36"))+
  facet_grid(cols = vars(sample_group), rows = vars(timepoint), space = "free", scales = "free")+
  labs(y = "Chromosome", x = "", fill = "Somy")+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        panel.grid = element_blank())


# #PAT_FS plot
# BGlist %>%
#   bind_bg_lists(table_to_bind = "somy_table") %>%
#   #mutate(sample_group = ifelse(sample_group == "control", "PBS", "382 µM SbIII"))%>%
#   #mutate(chromosome = factor(chromosome))%>%
#   mutate(chromosome = gsub("Ld", "", chromosome))%>%
#   mutate(chromosome = factor(chromosome, levels = rev(levels(factor(chromosome)))))%>%
#   filter(!population_name == "BD")%>%
#   ggplot(aes(x = timepoint, y = chromosome, fill = somy, label = round(somy, 2)))+
#   geom_tile()+
#   #geom_text(color = "white", size = 2)+
#   scale_fill_gradientn(colors = heat_col, limits = c(0,8))+
#   #scale_y_discrete(breaks = c("Ld01", "Ld06", "Ld11", "Ld16", "Ld21", "Ld26", "Ld31", "Ld36"))+
#   facet_grid(cols = vars(population_name), space = "free", scales = "free")+
#   labs(y = "Chromosome", x = "Passage", fill = "Somy")+
#   theme(axis.title = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         #axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         strip.text = element_text(size = 16),
#         panel.grid = element_blank())


# 
BGlist %>%
  bind_bg_lists(table_to_bind = "bg_table")%>%
  group_by(file, chromosome)%>%
  #filter(group == "382 uM SbIII")%>%
  filter(chromosome == "Ld23")%>%
  filter(timepoint == 5)%>%
  mutate(read_count = read_count/median(read_count))%>%
  mutate(name = gsub("P5D", "SePOP", name))%>%
  mutate(group = factor(group, levels = c("control", "382 uM SbIII")))%>%
  mutate(name = gsub("P5C", "cPOP", name))%>%
  mutate(chromosome = gsub("Ld", "Chromosome", chromosome))%>%
  ggplot(aes(x = chromo_bin, y = read_count, color = name))+
  geom_line(linewidth = 0.5)+
  scale_color_manual(values = iwanthue(8, cmin = 40, cmax = 70, lmin = 15, lmax = 85))+
  labs(x = "20kb bin", y = "normalized read count", color = "")+
  facet_grid(cols = vars(chromosome), rows = vars(group))

# #trajectory pca OBS: to be moved to another script later
# 
# sample_info <- read_excel(sample_info)
# 
# pca <- BGlist %>%
#   bind_bg_lists(table_to_bind = "somy_table") %>%
#   select(id, chromosome, somy) %>%
#   pivot_wider(names_from = chromosome, values_from = somy) %>%
#   as.data.frame
# 
# rownames(pca) <- pca$id
# pca <- pca[,-1]
# to_plot <- prcomp(pca)
# pca_metrics <- as.data.frame(summary(to_plot)$importance)[2,1:2]
# pca_metrics <- round(pca_metrics*100, digits = 2)
# to_plot <- as.data.frame(to_plot$x)
# to_plot$id <- rownames(pca)
# # to_plot <- as.data.frame(cbind(rownames(to_plot), to_plot))
# # colnames(to_plot)[1] <- "id"
# 
# to_plot <- cbind(to_plot, sample_info[match(to_plot$id, sample_info$id),-1])
# 
# 
# to_plot %>%
#   filter(group != "Before SbIII")%>%
#   mutate(path = paste(group, "rep.", replicate))%>%
#   arrange(timepoint)%>%
#   ggplot(aes(x = PC1, y = PC2, label = timepoint, color = path, group = path, fill = path))+
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey")+
#   geom_path(size = 1.5,
#             arrow = arrow(type = "open"))+
#   geom_text_repel(size = 7,
#                   #color = "black",
#                   show.legend = FALSE,
#                   min.segment.length = 0,
#                   nudge_x = .2)+
#   labs(x = paste0("PC1 (", as.numeric(pca_metrics[1]), "%)"),
#        y = paste0("PC2 (", as.numeric(pca_metrics[2]), "%)"),
#        color = "Sample",
#        fill = "Sample")+
#   scale_color_manual(values = c(RColorBrewer::brewer.pal(4, "Spectral"), RColorBrewer::brewer.pal(4, "Greys")))+
#   # scale_x_continuous(limits = c(-1,1))+
#   # scale_y_continuous(limits = c(-1,1))+
#   theme_bw()+
#   theme(panel.grid = element_blank(), text = element_text(size = 16))


# 
