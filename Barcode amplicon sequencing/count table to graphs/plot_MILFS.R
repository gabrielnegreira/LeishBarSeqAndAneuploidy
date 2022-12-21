#libraries
library(tidyverse)
#library(ggplot2)
library(ggrepel)
library(ggridges)
library(rlist)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(ggalluvial)

#inputs
#to_plot must be a list of BClist objects containing samples with group, timepoint and replicates. 
to_plot <- import_bc_lists(path = "inputs/MIL_FS", 
                           pattern = "cluster.csv",
                           input_type = "read_count",
                           sample_info = "inputs/MIL_FS/sample_info.xlsx",
                           barcode_col = "Center",
                           count_col = "time_point_1",
                           min_count = 0) %>%
  add_missing_bc()

#parameters
initial_only <- TRUE #if true, will only consider lineages that are found in the initial population.
drug_group <- "25 uM MIL"
control_group <- "control"
top_bc <- 1000000000
fold_change_threshold <- 2


#remove barcodes which are not in the initial population####
#it will first compare all dataframes agains the initial population to see which are present at the time point 0 (t0).
#obs: in this particular experiments all samples at t0 are identical so we can compare everyone against one of the t0 samples.
to_plot <- compare_bc_lists(to_plot, to_plot[["MILFS-P0C1"]], cname_prefix = "ini_") %>%
  bind_bc_lists() %>%
  filter(ini_count > 0)

#code####
#plotting dynamics####
#adding a random numeric identifier to each barcode
to_plot$bc_id <- as.factor(match(to_plot$barcode, sample(unique(to_plot$barcode))))
#selecting which barcodes to color. This will select only barcodes that reach a percentage of at least 1% in one of the replicates.
bc_to_color <- to_plot %>%
  filter(group == drug_group & timepoint == max(timepoint, na.rm = T))%>%
  group_by(bc_id)%>%
  summarise(max_percent = max(percent, na.rm = T))%>%
  ungroup()%>%
  filter(max_percent > 0.1) %>%
  select(bc_id) %>%
  unlist() %>%
  as.character() %>%
  unique()

#setting color palette
col_palette <- c("#55386b", "#b9a1cc", "#315195", "#E0F3F8", "#38526b" , "#386b67" ,"#81ba68", "#FEE090", "#FDAE61", "#F46D43", "#6b4c38" ,"#D73027", "#A50026", "#541b1b")
col_palette <- colorRampPalette(col_palette)

col_vec <- col_palette(length(bc_to_color))
#col_vec <- sample(col_vec) #to shuffle the colors each time the code is run
col_vec <- c(col_vec, "grey")
names(col_vec) <- c(bc_to_color, "other")

#choosing which lineages to add percentage labels in the following plot
to_plot$to_label <- FALSE
to_plot <- to_plot %>%
  select(barcode, bc_id, group, replicate, timepoint, percent)%>%
  pivot_wider(names_from = "timepoint", values_from = "percent")%>%
  mutate(to_label = `1` > 0.1) %>%
  pivot_longer(cols = c(`0`, `1`), names_to = "timepoint", values_to = "percent")
  

to_plot %>%
  #filter(group == drug_group)%>%
  mutate(color = ifelse(bc_id %in% names(col_vec), bc_id, "other"))%>%
  mutate(timepoint = ifelse(timepoint == 0, "day 0", "day 17"))%>%
  #mutate(replicate = ifelse(group == "control", paste0("McPOP", replicate),paste0("MePOP", replicate)))%>%
  mutate(label = ifelse(to_label, paste0(round(percent*100, 1), "%"), NA))%>%
  mutate(label = ifelse(group == control_group, NA, label))%>%
  mutate(group = ifelse(group == control_group, "McPOP", "MePOP"))%>%
  ggplot(aes(x= timepoint, y = percent, group = barcode, fill = color, alluvium = barcode, stratum = barcode, label = label))+
  geom_alluvium(alpha = 0.8, decreasing = FALSE)+
  geom_label_repel(stat = "stratum", 
                   #color = "#f2f2f2",
                   decreasing = FALSE, 
                   direction = "y", 
                   nudge_x = -0.15, 
                   nudge_y = 0.1,
                   show.legend = FALSE, 
                   #aes(size = abs(log(percent * 100, 2))),
                   min.segment.length = unit(0, 'lines'))+
  guides(fill = "none")+
  guides(fill=guide_legend(ncol=2))+
  scale_fill_manual(values = col_vec)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  facet_grid(cols = vars(replicate), rows = vars(group))+
  labs(x="", fill = "Lineage", y = "Frequency")+
  theme(#axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(3, "mm"), 
        legend.key.height = unit(1, "mm"), 
        text = element_text(size = 16))






