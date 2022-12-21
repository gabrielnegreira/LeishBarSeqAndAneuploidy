#libraries
library(tidyverse)
library(ggrepel)
library(ggridges)
library(rlist)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(ggalluvial)
library(patchwork)

#inputs
#to_plot must be a list of BClist objects containing samples with group, timepoint and replicates. 
to_plot <- readRDS("outputs/PAT_FS/BClist.rds") %>%
  add_missing_bc()

#parameters
initial_only <- TRUE #if true, will only consider lineages that are found in the initial population.
drug_group <- "382 uM SbIII"
control_group <- "control"
top_bc <- 1000000000
fold_change_threshold <- 2

#code
#first it will take all timepoints of the BClist objects of each group and each replicate and will compare them to the first timepoint
conditions <- lapply(to_plot,function(x){
  x <- x[-1]
}) %>%
  bind_rows()%>%
  select(id, group, replicate, timepoint)%>%
  group_by(id)%>%
  summarise(group = unique(group), replicate = unique(replicate), timepoint = unique(timepoint))

#comparing to the initial condition (t0)
list_to_plot <- compare_bc_lists(to_plot, to_plot[["PATFS_P0C1"]], cname_prefix = "ini_")

#comparing to control
for(t in unique(conditions$timepoint)){
    list <- list.filter(list_to_plot, timepoint == t)
    reference <- list.filter(list_to_plot, timepoint == t & group == control_group)
    indexes <- list.findi(list_to_plot, timepoint == t, n = length(list_to_plot))
    list_to_plot[indexes] <- compare_bc_lists(list, reference, cname_prefix = "di_")
}


#keeping only barcodes which are found in the initial population
if(initial_only){
  list_to_plot <- list_to_plot %>%
    lapply(function(x){
      x$bc_table <- filter(x$bc_table, ini_count > 0)
      x$bc_table$percent <- x$bc_table$count/sum(x$bc_table$count)
      return(x)
    })
}

#then it will convert the list into a single dataframe for plotting
to_plot <- bind_bc_lists(list_to_plot)
ini_fold_change_range <- c(min(to_plot$ini_fold_change[!is.infinite(to_plot$ini_fold_change)], na.rm = T),
                           max(to_plot$ini_fold_change[!is.infinite(to_plot$ini_fold_change)], na.rm = T))

di_fold_change_range <- c(min(to_plot$di_fold_change[!is.infinite(to_plot$di_fold_change)], na.rm = T),
                          max(to_plot$di_fold_change[!is.infinite(to_plot$di_fold_change)], na.rm = T))


to_plot <- to_plot %>%
  mutate(ini_fold_change = ifelse(ini_percent == 0 & percent == 0, NA, 
                                  ifelse(ini_percent > 0 & percent == 0, ini_fold_change_range[1]-5, 
                                         ifelse(ini_percent == 0 & percent > 0, ini_fold_change_range[2]+5, ini_fold_change)))) %>%
  mutate(di_fold_change = ifelse(di_percent == 0 & percent == 0, 0, 
                                  ifelse(di_percent > 0 & percent == 0, di_fold_change_range[1]-5, 
                                         ifelse(di_percent == 0 & percent > 0, di_fold_change_range[2], di_fold_change))))

#rm(list_to_plot)


to_plot <- filter(to_plot, group %in% c(control_group, drug_group))


#lineage diversity plot####
#this will plot the number of different lineages with a frequency higher than a given threshold in each condition and in each timepoint.
threshold <- round((1/453)/20, 4)
threshold <- 0
#threshold <- 0.0001
y_label <- ifelse(threshold == 0, "Number of barcodes", paste("Number of barcodes at frequency >", threshold))
to_plot %>%
  filter(count > 0) %>%
  filter(percent > threshold) %>%
  group_by(group, timepoint, replicate) %>%
  summarize(n_lineages = n()) %>%
  ggplot(aes(x = timepoint, y = n_lineages, color = replicate, group = replicate))+
  geom_jitter(width = 0.01, height = 0)+
  geom_line()+
  #labs(x = "Passage", y = paste0("number of lineages at > ", threshold*100, "%"))+
  labs(x = "Passage", y = y_label, color = "Replicate")+
  facet_grid(rows = vars(group))+
  scale_color_brewer(palette = "Spectral")


#skew plots####


#ridge plot: fold change####
to_plot %>%
  mutate(di_fold_change = ifelse(is.infinite(ini_fold_change), NA, ini_fold_change))%>%
  #filter(timepoint != 0) %>%
  ggplot(aes(x = ini_fold_change, y = timepoint, group = id))+
  geom_density_ridges(fill = "#CD534CFF", alpha = 0.75, scale = 3)+
  #geom_vline(xintercept = c(-skew_threshold,skew_threshold), color = "black", linetype = "dashed")+
  labs(y = "Passage", x = "Fold-change of barcode frequency (log2)")+
  guides(fill = "neutral")+
  scale_y_discrete(expand = expansion(mult = c(0.01, .5))) +
  scale_x_continuous(breaks = c(min(to_plot$ini_fold_change), -5, 0, 5), labels = c("-inf", "-5", "0", "5"))+
  facet_grid(cols = vars(replicate), rows = vars(group), space = "free")+
  #theme_ridges()+
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 16))

#plotting fold_change lines
to_plot %>%  
  filter(ini_position <= top_bc)%>%
  ggplot(aes(x = timepoint, y = ini_fold_change, group = barcode, color = factor(barcode)))+
  geom_line()+
  geom_hline(yintercept = -14, color = "red", linetype = "dashed")+
  scale_color_viridis_d()+
  guides(col = "none")+
  labs(x = "Passage", y = "Fold-change of barcode frequency (log 2)")+
  facet_grid(cols = vars(replicate), rows = vars(group))

#plotting the drug-induced (di) fold-change
#di-fold change is the fold change of a barcode in the drug group compared to the average of that barcode in the control group in the same timepoint
to_plot <- to_plot %>%
  filter(group %in% c(control_group, drug_group))%>%
  #filter(ini_count > 10)%>%
  mutate(effect = ifelse(is.na(di_fold_change), "similar",
                  ifelse(di_fold_change > fold_change_threshold, "higher", 
                  #ifelse(di_fold_change < di_fold_change_range[1], "lower",
                  ifelse(di_fold_change < -fold_change_threshold, "lower", "similar")))) %>%
  mutate(effect = factor(effect, levels = c("higher", 
                                            "similar", 
                                            #"lower", 
                                            "lower")))

my_breaks <- round(di_fold_change_range*0.2)/0.2
my_breaks <- c(min(to_plot$di_fold_change, na.rm = T),
               my_breaks[1], 
               -my_breaks[2],
               0, 
               my_breaks[2])
my_labels <- my_breaks
my_labels[1] <- "-inf"

p1 <- to_plot %>%
  filter(group == drug_group)%>%
  filter(di_fold_change != max(di_fold_change, na.rm = F))%>%
  mutate(replicate = paste0("SePOP", replicate))%>%
  ggplot(aes(x = timepoint, y = di_fold_change, color = effect, group = timepoint))+
    geom_jitter(aes(size = percent), shape = 1)+
    scale_color_manual(values = c("#a84832", "grey", "#3289a8", "#022763"))+
    scale_y_continuous(breaks = my_breaks, labels = my_labels)+
    labs(y = "fold-change (Log2) relative to control", x = "Passage", color = "Frequency\nCompared to\nControl", size = "Barcode\nFrequency")+
    #geom_line(aes(group = barcode))+
    #geom_violin(aes(y = violin), color = "black", alpha = 0.2, draw_quantiles = c(0.05, 0.95))+
    facet_grid(cols = vars(replicate), rows = vars(group))+
    theme(text = element_text(size = 18), axis.title = element_text(size = 16))

#checking which are the common enriched lineages
to_plot %>%
  filter(timepoint >0)%>%
  filter(ini_count > 0)%>%
  filter(group == drug_group) %>%
  group_by(barcode, group, timepoint, effect) %>%
  summarise(n_replicates = n()) %>%
  group_by(group, timepoint, n_replicates, effect) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = timepoint, y = count, fill = factor(n_replicates), label = count))+
  geom_col(position = "fill")+
  #geom_col()+
  geom_text(position = position_fill(vjust = 0.5))+
  labs(fill = "Number of\nreplicates", x = "Passage", y = "Percentage of barcodes")+
  scale_fill_brewer(palette = "RdYlBu", direction = -1)+
  facet_grid(rows = vars(effect), cols = vars(group))

#heatmap plot
#this will plot a heatmap showing the drug-induced fold change of the top x barcodes in the last timepoint

heat_plot <- to_plot %>%
  #filter(group == drug_group)%>%
  filter(timepoint == max(timepoint))%>%
  mutate(replicate = paste(group, replicate, sep = "_"))%>%
  select(barcode, replicate, di_fold_change)%>%
  pivot_wider(names_from = replicate, values_from = di_fold_change)%>%
  column_to_rownames("barcode")

paletteLength <- 100
myColor <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(paletteLength)
myBreaks <- c(seq(min(heat_plot), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heat_plot)/paletteLength, max(heat_plot), length.out=floor(paletteLength/2)))

heat_plot %>%
  pheatmap(show_rownames = F,
           #color = c("#a84832", "grey", "#3289a8", "#022763"),
           #breaks = c(0:4),
           clustering_method = "ward.D2",
           #cutree_rows = 4,
           cutree_cols = 2,
           breaks = myBreaks,
           color = myColor)

#estimating the initial percentage of barcodes of each effect
p2 <- to_plot %>%
  group_by(effect, group, timepoint, replicate) %>%
  summarize(ini_percent = sum(ini_percent))%>%
  mutate(label = ifelse(ini_percent < 0.01, NA, paste0(round(ini_percent*100,1), "%")))%>%
  filter(group == drug_group)%>%
  ggplot(aes(x = timepoint, y = ini_percent, fill = effect, label = label))+
  geom_col()+
  geom_text(position = position_stack(vjust = 0.5), size = 2.5, show.legend = F)+
  labs(x = "Passage", y = "Fraction of the initial population", fill = "Frequency\nCompared to\nControl")+
  scale_fill_manual(values = c("#a84832", "grey", "#3289a8", "#022763"))+
  facet_grid(cols = vars(replicate), rows = vars(group))#+
  #theme(text = element_text(size = 24))

plot_grid(p1,p2, ncol = 1, align = "v")

#comparing the frequency of each group of barcodes. 
to_plot$final_effect <- NA
to_plot$index <- c(1:nrow(to_plot))
for(g in unique(to_plot$group)){
  for(r in unique(to_plot$replicate)){
    final <- filter(to_plot, group == g & replicate == r & timepoint == max(to_plot$timepoint))
    for(t in unique(to_plot$timepoint)){
      print(paste(g,r,t))
      data <- filter(to_plot, group == g & timepoint == t & replicate == r)
      rows <- as.numeric(data$index)
      to_plot$final_effect[rows] <- as.character(final$effect[match(data$barcode, final$barcode)])
    }
  }
}
to_plot$final_effect <- factor(to_plot$final_effect, levels = levels(to_plot$effect))

to_plot %>%
  #filter(timepoint %in% c(0,5))%>%
  filter(group == drug_group)%>%
  group_by(group, replicate, timepoint, final_effect)%>%
  summarize(percent = sum(percent)) %>%
  mutate(label = ifelse(percent < 0.01, NA, paste0(round(percent*100,1), "%")))%>%
  mutate(label = ifelse(timepoint %in% c(0,5), label, NA))%>%
  mutate(replicate = paste0("SePOP", replicate))%>%
  ggplot(aes(x = as.numeric(timepoint), y = percent, fill = final_effect, label = label))+
  #geom_col()+
  geom_area()+
  geom_label_repel(position = position_stack(vjust = 0.5), size = 3.5, show.legend = F)+
  scale_fill_manual(values = c("#a84832", "grey", "#3289a8", "#022763"), labels = c("higher than control", "similar to control", "lower than control", "lower"))+
  facet_grid(cols = vars(replicate), rows = vars(group))+
  labs(x= "Passage", y = "Percentage", fill = "Final Frequency (P5)")+
  theme(text = element_text(size = 14))



#plotting dynamics####
#adding a random numeric identifier to each barcode
to_plot$bc_id <- as.factor(match(to_plot$barcode, sample(unique(to_plot$barcode))))
bc_to_color <- to_plot %>%
  filter(group == drug_group & timepoint == max(timepoint, na.rm = T))%>%
  group_by(bc_id)%>%
  summarise(max_percent = max(percent, na.rm = T))%>%
  ungroup()%>%
  filter(max_percent > 0.01) %>%
  select(bc_id) %>%
  unlist() %>%
  as.character() %>%
  unique()


#TO REMOVE LATER
# a <- to_plot %>%
#   filter(timepoint == 5, group == drug_group)%>%
#   group_by(bc_id, final_effect)%>%
#   summarise(replicates = n())%>%
#   pivot_wider(names_from = final_effect, values_from = replicates)
# a[is.na(a)] <- 0
# 
# bc_to_color <- a %>%
#   mutate(sum = lower+lower)%>%
#   filter(higher >= 1)%>%
#   select(bc_id)%>%
#   unlist()%>%
#   as.character()

#setting color palette
col_palette <- c("#55386b", "#b9a1cc", "#315195", "#E0F3F8", "#38526b" , "#386b67" ,"#81ba68", "#FEE090", "#FDAE61", "#F46D43", "#6b4c38" ,"#D73027", "#A50026", "#541b1b")
col_palette <- colorRampPalette(col_palette)

col_vec <- col_palette(length(bc_to_color))
#col_vec <- sample(col_vec) #to shuffle the colors each time the code is run
col_vec <- c(col_vec, "grey")
names(col_vec) <- c(bc_to_color, "other")

#choosing which lineages to add percentage labels in the following plot
to_plot2 <- to_plot
to_plot2$to_label <- FALSE
to_plot2 <- to_plot2 %>%
  select(barcode, bc_id, group, replicate, timepoint, percent)%>%
  pivot_wider(names_from = "timepoint", values_from = "percent")%>%
  mutate(to_label = `5` > 0.10) %>%
  pivot_longer(cols = c(`0`, `1`, `2`, `3`, `4`, `5`), names_to = "timepoint", values_to = "percent")

plot2 <- to_plot2 %>%
  #filter(group == drug_group)%>%
  mutate(color = ifelse(bc_id %in% names(col_vec), bc_id, "other"))%>%
  mutate(label = ifelse(to_label, paste0(round(percent*100, 1), "%"), NA))%>%
  mutate(label = ifelse(group == control_group, NA, label))%>%
  mutate(label = ifelse(timepoint %in% c(0,5), label, NA))%>%
  mutate(percent = ifelse(is.na(percent), 0, percent))%>%
  mutate(group = ifelse(group == control_group, "ScPOP", "SePOP"))%>%
  ggplot(aes(x= timepoint, y = percent, group = barcode, fill = color, alluvium = barcode, stratum = barcode, label = label))+
  guides(fill = "none")+
  guides(fill=guide_legend(ncol=2))+
  scale_fill_manual(values = col_vec)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  facet_grid(cols = vars(replicate), rows = vars(group))+
  labs(x="", fill = "Clone", y = "Frequency")+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.key.width = unit(3, "mm"), 
        legend.key.height = unit(1, "mm"), 
        text = element_text(size = 16))

#bump plot
plot2 <- plot2+geom_alluvium(alpha = 0.8, decreasing = FALSE)+
  geom_label_repel(stat = "stratum", 
                   #color = "#f2f2f2",
                   decreasing = FALSE, 
                   direction = "y", 
                   nudge_x = -0.15, 
                   nudge_y = 0.1,
                   show.legend = FALSE, 
                   #aes(size = abs(log(percent * 100, 2))),
                   min.segment.length = unit(0, 'lines'))
plot2

#combining with aneuploidy plot (to remove later)
load("pat_fs_bc_somies.RData")
plot <- plot+
  theme(axis.text.y = element_text(size = 6))+
  scale_x_discrete(expand = c(0,0))
plot2/plot+plot_layout(heights = c(1,1))

#alluvial comparing replicates at last timepoint
to_plot %>%
  filter(timepoint == max(timepoint))%>%
  #filter(timepoint == 5)%>%
  #filter(group == control_group)%>%
  filter(group == drug_group)%>%
  #mutate(percent = 1-percent)%>%
  #mutate(percent = sqrt(percent))%>%
  #mutate(replicate = factor(replicate, levels = c(1,4,2,3)))%>%
  group_by(replicate, effect)%>%
  mutate(percent = 1/n())%>%
  ggplot(aes(x = replicate,
             #y = percent,
             alluvium = bc_id, 
             stratum = effect, 
             fill = effect, 
             label = effect))+
  #geom_alluvium()+
  geom_flow()+
  geom_stratum()+
  geom_text(stat = "stratum", size = 4, color = "#f5f5f5")+
  scale_fill_manual(values = c("#a84832", "grey", "#3289a8", "#022763"))+
  #scale_fill_manual(values = col_vec)+
  guides(fill = "none")+
  theme_linedraw()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        )


#ven diagram for each effect
effects <- as.character(unique(to_plot$effect))

a <- to_plot %>%
  filter(group == drug_group)%>%
  filter(timepoint == max(timepoint))

b <- list()
for(i in c(1:length(effects))){
  effect_to_plot <- effects[i]
  print(effect_to_plot)
  b[[i]] <- a %>%
    filter(effect == effect_to_plot)%>%
    select(barcode, replicate)%>%
    pivot_wider(names_from = replicate, values_from = barcode, values_fn = "list")%>%
    as.list()%>%
    lapply(unlist)%>%
    fromList()%>%
    upset(sets.x.label = "Clones per replicate", mainbar.y.label = "Number of Clones", empty.intersections = "on")
}

b


#coordiagram at last timepoint
#this assumes that there is one common initial population that is later divided in multiple replicates
a <- to_plot %>%
  filter(group == drug_group)%>%
  #filter(count >0)%>%
  filter(timepoint  == max(timepoint) | replicate == 1 & timepoint == 0) %>%
  select(bc_id, timepoint, replicate, percent)%>%
  mutate(direction = ifelse(timepoint == 0, 0, replicate))%>%
  group_by(bc_id)%>%
  summarise(from = min(direction), to = direction, value_from = percent[which(direction == min(direction))] ,value_to = percent)%>%
  filter(to != from)%>%
  mutate(value_from = value_from/4)%>%
  #mutate(value_to = value_to/value_from)%>%
  mutate(color = col_vec[match(bc_id, names(col_vec))])%>%
  mutate(color = ifelse(is.na(color), "light grey", color))

gap <- c(50, 10, 10, 10, 50)
circos.par(start.degree = 293, gap.after = gap)
chordDiagram(a[,-c(1, ncol(a))],
             scale = FALSE, 
             col = a$color,
             #reduce = 0,
             #group = group,
             #grid.col = to_plot$color,
             #directional = 2,
             link.sort = TRUE,
             link.decreasing = TRUE)
circos.clear()






