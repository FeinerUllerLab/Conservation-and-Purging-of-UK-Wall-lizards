rm(list = ls())

##### Load libraries ######

library(openxlsx)
library(dplyr)
library(tidyr)
library(tibble)
library(geosphere) # need this package
library(ggplot2)
library(tmap)
library(sf)
library(raster)
library(rnaturalearth)
library(RColorBrewer)
library(ggforce)
library(scatterpie)
library(ggiraph)
library(PieGlyph)
library(ggspatial)
library(tinytable)
library(ggrepel)
library(patchwork)
library(ape)
library(stringr)
library(paletteer)
library(scales)
library(effsize)
library(lme4)
library(ggtree)

### Plot maps for Fig. 1


# Read the Excel file - Get one sample per location ----
LizardsData <- readWorkbook("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Samples_Info.xlsx", sheet = 1)
LizardsData$Group <- paste(LizardsData$Origin,"-",LizardsData$NativeVSIntro, sep="")
UKLizards <- subset(LizardsData, NativeVSIntro=="Intro")
UKLizards <- UKLizards %>%
  group_by(Longitude, Latitude) %>%
  slice(1)
NativeLizards <- subset(LizardsData, NativeVSIntro=="Native")

# Create the maps ----

# Species distribution:
world <- ne_countries(scale = "large", returnclass = "sf") #change to large for final plot!
sf_use_s2(FALSE)
europe_cropped <- st_crop(world, xmin = -10.00, xmax = 35.00, ymin = 35.00, ymax = 52.00)
europe_cropped_native <- st_crop(world, xmin = -5, xmax = 14, ymin = 41, ymax = 49)
europe_cropped_intro <- st_crop(world, xmin = -6, xmax = 2.00, ymin = 49.5, ymax = 52.00)

#read in distribution map of P. muralis from IUCN
Pmur_dist <- st_read("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Pmuralis_distribution/data_0.shp")
map_overview <- ggplot(data = europe_cropped) + 
  geom_sf(data = europe_cropped, fill= "white", color = "grey", size = 0.1) +
  geom_sf(data = Pmur_dist, fill = "lightgrey", color = NA) +
  coord_sf(xlim = c(-10.00, 35.00), ylim = c(35.00, 52.00), expand = F) +
  annotate("rect", xmin = -5, xmax = 14, ymin = 41, ymax = 49,color = "black", fill=NA) +
  annotate("rect", xmin = -6, xmax = 2, ymin = 49.5, ymax = 52,color = "black", fill=NA) +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(height=unit(0.2,"cm"), bar_cols=c("black","white"), pad_x=unit(0.2,"cm"), pad_y=unit(0.3,"cm")) +
  labs(subtitle = "A") +
  theme(panel.background = element_rect(fill = "aliceblue"), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.25))

# W-EUR + C-ITA
map_native <- ggplot(data = europe_cropped_native) +
  geom_sf(fill = "white", color = "grey", size=0.1) +
  geom_sf(data = st_crop(Pmur_dist, xmin = -5, xmax = 14, ymin = 41, ymax = 49), fill = "lightgrey", color = NA) +
  geom_point(data = NativeLizards, aes(x = Longitude, y = Latitude, colour = Origin),
             size = 5, shape = 20) +
  coord_sf(xlim = c(-5, 14), ylim = c(41, 49), expand = FALSE) +
  scale_colour_manual(values = c("W-EUR" = "#DD6A27", "C-ITA" = "#0673B3")) +
  geom_text_repel(data = NativeLizards, aes(label = Abbpop, x = Longitude, y = Latitude),
                  box.padding = 0.7,
                  point.padding = 0,
                  size = 5,
                  segment.color = 'grey50') +
  annotation_scale(height=unit(0.2,"cm"), bar_cols=c("black","white"), pad_x=unit(0.2,"cm"), pad_y=unit(0.3,"cm")) +
  labs(subtitle = "B") +
  theme(panel.background = element_rect(fill = "aliceblue"), panel.border = element_rect(colour = "black", 
       fill=NA, linewidth=0.25), axis.title = element_blank(), panel.grid.major = element_blank(),
       axis.text  = element_blank(),
       axis.ticks = element_blank()) 

#UK
map_UK <- ggplot(data = europe_cropped_intro) +
  geom_sf(fill = "white", color = "grey", size = 0.1) +
  geom_point(data = UKLizards, aes(x = Longitude, y = Latitude, colour = Origin),
             size = 5, shape = 20) +
  coord_sf(xlim = c(-6.0, 2.0), ylim = c(49.5, 52.0), expand = FALSE) +
  scale_colour_manual(values = c("W-EUR" = "#EFE808", "C-ITA" = "#6BCBDA")) +
  geom_text_repel(data = UKLizards, aes(label = Abbpop, x = Longitude, y = Latitude),
                  box.padding = 0.7,
                  point.padding = 0,
                  size = 5,
                  segment.color = 'grey50') +
  annotation_scale(height=unit(0.2,"cm"), bar_cols=c("black","white"), pad_x=unit(0.2,"cm"), pad_y=unit(0.3,"cm")) +
  labs(subtitle = "C") +
  theme(panel.background = element_rect(fill = "aliceblue"), panel.border = element_rect(colour = "black", 
        fill=NA, linewidth=0.25), axis.title = element_blank(), panel.grid.major = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank()) 

# Combine the maps and display them together
maps <- map_overview + map_native + map_UK + plot_layout(ncol = 1)
pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/DistributionMap.pdf", height=10, width=10, useDingbats = F)
print(maps)
dev.off()


### Plot PCA and phylogeny for Fig. 2

# Using the PCA outputs from Plink
eval = read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/All_Origins_PCA_Final.eigenval", header=F)
evec = read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/All_Origins_PCA_Final.eigenvec", header=F)

# Merge the data with the origin. 
colnames(evec) <- c("ID", "RM", paste0("PC", 1:(ncol(evec) - 2)))
# Remove or ignore the redundant ID column
evec <- evec[, -2]  # removes the second column

# Calculate percentage variance explained for each PC
pve <- data.frame(PC = 1:20, pve = eval / sum(eval) * 100)
evec1.pc <- round(eval[1,1]/sum(eval)*100,digits=2)
evec2.pc <- round(eval[2,1]/sum(eval)*100,digits=2)
evec3.pc <- round(eval[3,1]/sum(eval)*100,digits=2)
evec4.pc <- round(eval[4,1]/sum(eval)*100,digits=2)
evec5.pc <- round(eval[5,1]/sum(eval)*100,digits=2)

# Merge PCA data and location
pca_data <- merge(evec, LizardsData, by = "ID")

# 1.1.2 Plots including individual names and origin ----
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 6, alpha=0.7, shape=16) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC2 (", evec2.pc, "%)", sep = "")) +
  scale_color_manual(
    values = c("C-ITA-Native" = "#0673B3", 
               "W-EUR-Native" = "#DD6A27", 
               "C-ITA-Intro" = "#6BCBDA", 
               "W-EUR-Intro" = "#EFE808"),
    labels = c("C-ITA-Native" = "Native C-ITA", 
               "W-EUR-Native" = "Native W-EUR", 
               "C-ITA-Intro" = "Non-native C-ITA", 
               "W-EUR-Intro" = "Non-native W-EUR")) +
  labs(subtitle = "A") +
  theme_bw() +
  theme(legend.position = "bottom", axis.title = element_text(size = 12))

# Create p2 without legend but updated labels
p2 <- ggplot(pca_data, aes(x = PC1, y = PC3, color = Group)) +
  geom_point(size = 6, alpha=0.7, shape=16) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC3 (", evec3.pc, "%)", sep = "")) +
  scale_color_manual(
    values = c("C-ITA-Native" = "#0673B3", 
               "W-EUR-Native" = "#DD6A27", 
               "C-ITA-Intro" = "#6BCBDA", 
               "W-EUR-Intro" = "#EFE808"),
    labels = c("C-ITA-Native" = "Native C-ITA", 
               "W-EUR-Native" = "Native W-EUR", 
               "C-ITA-Intro" = "Non-native C-ITA", 
               "W-EUR-Intro" = "Non-native W-EUR")) +
  labs(subtitle = "B") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 12))

# Combine the plots with a shared label
PCA_plot <- p1 + p2 + 
  plot_layout(guides = "collect") + 
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, size = 15))) &
  theme(legend.position = "bottom")

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/PCA.pdf", height=4, width=8, useDingbats = F)
print(PCA_plot)
dev.off()

## 1.2. TREE ----

# 1.2.1 Data preparation ----
#Load the tree data
tree <- read.tree("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/Final_tree.treefile")

# Convert tree to tibble format
tree_data <- as_tibble(tree)

# Get correct node numbers for internal nodes
n_tips <- length(tree$tip.label)
internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)

# Create bootstrap data frame with proper node numbering
bootstrap_data <- data.frame(node = internal_nodes, bootstrap = as.numeric(tree$node.label))# Convert to numeric

# Merge with main tree data
tree_data <- tree_data %>% 
  left_join(bootstrap_data, by = "node") %>% 
  # Merge with lizard metadata (tips only)
  left_join(LizardsData, by = c("label" = "ID"))

# Create the tree plot
treeplot <- ggtree(tree, layout = "unrooted") %<+% tree_data +
  geom_tiplab2(aes(label = ID_alias, color = Group), size = 3, offset = 0.0, align = FALSE, show.legend = FALSE) +
  geom_tippoint(aes(color = Group), size = 2) +
  geom_nodelab(aes(label = round(bootstrap)), 
               color = "black", 
               hjust = -0.05,
               size = 2, 
               na.rm = TRUE,) +
  scale_color_manual(values = c("C-ITA-Native" = "#0673B3", 
                                "W-EUR-Native" = "#DD6A27", 
                                "C-ITA-Intro" = "#6BCBDA", 
                                "W-EUR-Intro" = "#EFE808"),
                     labels = c("Non-native W-EUR", "Non-native C-ITA", 
                                "Native W-EUR", "Native C-ITA")) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right", panel.grid.major = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)  )

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Tree_V1.pdf", height=10, width=10, useDingbats = F)
print(treeplot)
dev.off()

## 1.3. ADMIXTURE for Supplementary Fig. S1

# 1.3.1 Plot the CV error from the ADMIXTURE analysis----
cv_error <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/CV_errors_summary_Final.txt", h=F)
colnames(cv_error)<- c("rm","rm2","K","CV_error")
cv_error <- cv_error %>%
  mutate(K = str_extract(K, "\\d+"))
cv_error <- cv_error %>% dplyr::select(-rm, -rm2)
plot(cv_error) # Best supported K = 2 , 3 and 4


# 1.3.2 Organize the Data of K = 2,3 and 4 ----

# Load Admix data and merge with organized population data 

Lizards_admix<- readWorkbook("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Samples_Info.xlsx", sheet = 3)
K2<- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/K2_Final.Q", header=F)
K3<- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/K3_Final.Q", header=F)
K4<- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/K4_Final.Q", header=F)

# Create proper names for the columns in the Ks data sets and ID per sample.
K2<- cbind(Lizards_admix, K2)
colnames(K2) <- c("ID","ID_alias","Origin","Abbpop","Q1","Q2")
K3<- cbind(Lizards_admix, K3)
colnames(K3) <- c("ID","ID_alias","Origin","Abbpop","Q1","Q2","Q3")
K4<- cbind(Lizards_admix, K4)
colnames(K4) <- c("ID","ID_alias","Origin","Abbpop","Q1","Q2","Q3","Q4")

# Organize samples based on K2 (Best supported)
ordered_individuals <- K2 %>%
  arrange(-Q1) %>%
  pull(ID_alias)  

# Apply ordering before pivoting
K2$ID_alias <- factor(K2$ID_alias, levels = ordered_individuals)
K3$ID_alias <- factor(K3$ID_alias, levels = ordered_individuals)
K4$ID_alias <- factor(K4$ID_alias, levels = ordered_individuals)

# Convert to long format
K2_long <- K2 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")
K3_long <- K3 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")
K4_long <- K4 %>% pivot_longer(cols = starts_with("Q"), names_to = "Ancestry", values_to = "Q_value")

# Add K identifier
K2_long$K <- "K2"
K3_long$K <- "K3"
K4_long$K <- "K4"

# Combine after ordering
admix_data <- bind_rows(K2_long, K3_long, K4_long)

# Re-apply factor ordering
admix_data$ID_alias <- factor(admix_data$ID_alias, levels = ordered_individuals)

# 1.3.3 Plot the ADMIXTURE ----
pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Admixture.pdf", height=6, width=10, useDingbats = F)
ggplot(admix_data, aes(x = ID_alias, y = Q_value, fill = Ancestry)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~K, ncol = 1) +
  scale_fill_manual(values = as.vector(paletteer_d("ggsci::default_jco"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Admixture Proportion", x = "Individual")
dev.off()


#### 2. GENETIC DIVERSITY AND INBREEDING ####

## 2.1. HETEROZYGOSITY + ROH ----

RoH_C_ITA <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_bcftools_C-ITA_2Mb.txt", header = T)
RoH_C_ITA <- merge(RoH_C_ITA, Lizards_admix, by.x = "Sample", by.y = "ID")

geno_Het_C_ITA <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/heterozygosity_summary_genomewide_C-ITA.tsv", header = T)
geno_Het_C_ITA$O.HET <- geno_Het_C_ITA$NumberVariableSites-geno_Het_C_ITA$ObservedHomozygous
geno_Het_C_ITA$geno_F.HET <- geno_Het_C_ITA$O.HET/(geno_Het_C_ITA$GenotypedSites-geno_Het_C_ITA$MissingSites)
geno_Het_C_ITA <- geno_Het_C_ITA[,-(2:6)]
RoH_HET_C_ITA <- merge(RoH_C_ITA, geno_Het_C_ITA, by = "Sample")

nonRoH_Het_C_ITA <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/heterozygosity_summary_nonRoH_C-ITA_2Mb.tsv", header = T)
nonRoH_Het_C_ITA$O.HET <- nonRoH_Het_C_ITA$NumberVariableSites-nonRoH_Het_C_ITA$ObservedHomozygous
nonRoH_Het_C_ITA$nonRoH_F.HET <- nonRoH_Het_C_ITA$O.HET/(nonRoH_Het_C_ITA$GenotypedSites-nonRoH_Het_C_ITA$MissingSites)
nonRoH_Het_C_ITA$TotalNumberSites <- nonRoH_Het_C_ITA$GenotypedSites-nonRoH_Het_C_ITA$MissingSites
nonRoH_Het_C_ITA <- nonRoH_Het_C_ITA[,-c(2:6,8)]
RoH_Het_C_ITA <- merge(RoH_HET_C_ITA, nonRoH_Het_C_ITA, by = "Sample")
RoH_Het_C_ITA$IDRisk <- RoH_Het_C_ITA$FRoH*RoH_Het_C_ITA$nonRoH_F.HET*1000

RoH_Het_C_ITA %>%
  group_by(Origin) %>%
  summarise(
    mean_het   = mean(geno_F.HET, na.rm = TRUE),
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n()
  )

p1 <- ggplot(RoH_HET_C_ITA, aes(y = geno_F.HET, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0.001,0.0045) + theme_bw()
p2 <- ggplot(RoH_Het_C_ITA, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  geom_text_repel(data = subset(RoH_Het_C_ITA, Origin == "Int-C-ITA"), aes(label = ID_alias), color="black") +
  theme_bw()
p3 <- ggplot(RoH_HET_C_ITA, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0,0.125) + theme_bw()
p1RISK <- ggplot(RoH_Het_C_ITA, aes(y = nonRoH_F.HET, x = FRoH, color=IDRisk)) +
  geom_point(size=3) + 
  scale_color_gradient2(
    low      = "#FFCE03",   # very light salmon
    mid      = "orange",   # medium pink‑red
    high     = "#de2d26",   # dark red
    midpoint = median(RoH_Het_C_ITA$IDRisk, na.rm = TRUE),   # or any value you choose
    space    = "Lab"       # smoother perceptual interpolation
  ) +
  theme_bw() +theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))

# Same for W-EUR
RoH_W_EUR <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_bcftools_W-EUR_2Mb.txt", header = T)
RoH_W_EUR <- merge(RoH_W_EUR, Lizards_admix, by.x = "Sample", by.y = "ID")

geno_Het_W_EUR <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/heterozygosity_summary_genomewide_W-EUR.tsv", header = T)
geno_Het_W_EUR$O.HET <- geno_Het_W_EUR$NumberVariableSites-geno_Het_W_EUR$ObservedHomozygous
geno_Het_W_EUR$geno_F.HET <- geno_Het_W_EUR$O.HET/(geno_Het_W_EUR$GenotypedSites-geno_Het_W_EUR$MissingSites)
geno_Het_W_EUR <- geno_Het_W_EUR[,-(2:6)]
RoH_HET_W_EUR <- merge(RoH_W_EUR, geno_Het_W_EUR, by = "Sample")

nonRoH_Het_W_EUR <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/heterozygosity_summary_nonRoH_W-EUR_2Mb.tsv", header = T)
nonRoH_Het_W_EUR$O.HET <- nonRoH_Het_W_EUR$NumberVariableSites-nonRoH_Het_W_EUR$ObservedHomozygous
nonRoH_Het_W_EUR$nonRoH_F.HET <- nonRoH_Het_W_EUR$O.HET/(nonRoH_Het_W_EUR$GenotypedSites-nonRoH_Het_W_EUR$MissingSites)
nonRoH_Het_W_EUR$TotalNumberSites <- nonRoH_Het_W_EUR$GenotypedSites-nonRoH_Het_W_EUR$MissingSites
nonRoH_Het_W_EUR <- nonRoH_Het_W_EUR[,-c(2:6,8)]
RoH_Het_W_EUR <- merge(RoH_HET_W_EUR, nonRoH_Het_W_EUR, by = "Sample")
RoH_Het_W_EUR$IDRisk <- RoH_Het_W_EUR$FRoH*RoH_Het_W_EUR$nonRoH_F.HET*1000

RoH_Het_W_EUR %>%
  group_by(Origin) %>%
  summarise(
    mean_het   = mean(geno_F.HET, na.rm = TRUE),
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n()
  )

p4 <- ggplot(RoH_HET_W_EUR, aes(y = geno_F.HET, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0.001,0.0045) + theme_bw()
p5 <- ggplot(RoH_Het_W_EUR, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  geom_text_repel(data = subset(RoH_Het_W_EUR, Origin == "Int-W-EUR"), aes(label = ID_alias), color="black") +
  theme_bw()
p6 <- ggplot(RoH_HET_W_EUR, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0,0.125) + theme_bw()
p2RISK <- ggplot(RoH_Het_W_EUR, aes(y = nonRoH_F.HET, x = FRoH, color=IDRisk)) +
  geom_point(size=3) + 
  scale_color_gradient2(
    low      = "#FFCE03",   # very light salmon
    mid      = "orange",   # medium pink‑red
    high     = "#de2d26",   # dark red
    midpoint = median(RoH_Het_W_EUR$IDRisk, na.rm = TRUE),   # or any value you choose
    space    = "Lab"       # smoother perceptual interpolation
  ) +
  theme_bw() +theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/HetRoH_V3.pdf", height=12, width=10, useDingbats = F)
(p1+p4)/(p2+p5)/(p3+p6)
dev.off()

### Statistics
## Heterozygosity
# Wilcoxon test
wilcox.test(geno_F.HET ~ Origin, data = RoH_HET_C_ITA)
# Cliff's delta
cliff.delta(subset(RoH_HET_C_ITA, Origin == "Nat-C-ITA")$geno_F.HET, subset(RoH_HET_C_ITA, Origin == "Int-C-ITA")$geno_F.HET)
# Wilcoxon test
wilcox.test(geno_F.HET ~ Origin, data = RoH_HET_W_EUR)
# Cliff's delta
cliff.delta(subset(RoH_HET_W_EUR, Origin == "Nat-W-EUR")$geno_F.HET, subset(RoH_HET_W_EUR, Origin == "Int-W-EUR")$geno_F.HET)

## ROH Length
# Length of ROHs comparison 
C_ITA_BCF_model <- glmer(Length ~ Origin + (1 | Sample), 
                       family = Gamma(link = "log"), 
                       data = RoH_HET_W_EUR)
summary(C_ITA_BCF_model)

# Get the proportion of the effect 
exp(fixef(C_ITA_BCF_model)[-1]) #  ~ meaning 0.7X difference.

# Confidence intervals for the model 
# Function for bootstrapping 1000 replicates 
get_boot_ci <- function(model, term, nsim = 1000) {
  boot_fun <- function(.) {
    fixef(.)[term]
  }
  boot_ci <- bootMer(model, 
                     FUN = boot_fun, 
                     nsim = nsim,
                     type = "parametric")
  ci <- boot.ci(boot_ci, type = "perc", conf = 0.95)$percent[4:5]
  return(exp(ci))
}

## FROH
wilcox.test(FRoH ~ Origin, data = RoH_HET_C_ITA)
wilcox.test(FRoH ~ Origin, data = RoH_HET_W_EUR)

#####
# ROH for Supplements (Plink + different length)
#####

RoH_C_ITA_Supp <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_bcftools_C-ITA_500kb.txt", header = T)
RoH_C_ITA_Supp <- merge(RoH_C_ITA_Supp, Lizards_admix, by.x = "Sample", by.y = "ID")
RoH_C_ITA_Supp$Length <- RoH_C_ITA_Supp$Length/1000000
p1S <- ggplot(RoH_C_ITA_Supp, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  geom_text_repel(data = subset(RoH_C_ITA_Supp, Origin == "Int-C-ITA"), aes(label = ID_alias), color="black") +
  theme_bw()
p2S <- ggplot(RoH_C_ITA_Supp, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0,0.31) + theme_bw()

RoH_W_EUR_Supp <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_bcftools_W-EUR_500kb.txt", header = T)
RoH_W_EUR_Supp <- merge(RoH_W_EUR_Supp, Lizards_admix, by.x = "Sample", by.y = "ID")
RoH_W_EUR_Supp$Length <- RoH_W_EUR_Supp$Length/1000000
p3S <- ggplot(RoH_W_EUR_Supp, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  geom_text_repel(data = subset(RoH_W_EUR_Supp, Origin == "Int-W-EUR"), aes(label = ID_alias), color="black") +
  theme_bw()
p4S <- ggplot(RoH_W_EUR_Supp, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0,0.31) + theme_bw()

RoH_C_ITA_plink_2Mb<- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_plink_C-ITA_2Mb.hom.indiv", header = T)
RoH_C_ITA_plink_2Mb$FRoH <- RoH_C_ITA_plink_2Mb$KB/RoH_C_ITA_plink_2Mb$KB[RoH_C_ITA_plink_2Mb$FID == "Pxx"]
RoH_C_ITA_plink_2Mb <- subset(RoH_C_ITA_plink_2Mb, FID != "Pxx")
RoH_C_ITA_plink_2Mb$Length <- RoH_C_ITA_plink_2Mb$KB/1000
colnames(RoH_C_ITA_plink_2Mb) <- c("Sample","IID","PHE","nRoH","KB","KBAVG","FRoH","Length")
RoH_C_ITA_plink_2Mb <- merge(RoH_C_ITA_plink_2Mb, Lizards_admix, by.x = "Sample", by.y = "ID")
p_pl_C_ITA_2Mb_1 <- ggplot(RoH_C_ITA_plink_2Mb, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  geom_text_repel(data = subset(RoH_C_ITA_plink_2Mb, Origin == "Int-C-ITA"), aes(label = ID_alias), color="black") +
  theme_bw()
p_pl_C_ITA_2Mb_2 <- ggplot(RoH_C_ITA_plink_2Mb, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0,0.18) + theme_bw()

RoH_W_EUR_plink_2Mb <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_plink_W-EUR_2Mb.hom.indiv", header = T)
RoH_W_EUR_plink_2Mb$FRoH <- RoH_W_EUR_plink_2Mb$KB/RoH_W_EUR_plink_2Mb$KB[RoH_W_EUR_plink_2Mb$FID == "Pxx"]
RoH_W_EUR_plink_2Mb <- subset(RoH_W_EUR_plink_2Mb, FID != "Pxx")
RoH_W_EUR_plink_2Mb$Length <- RoH_W_EUR_plink_2Mb$KB/1000
colnames(RoH_W_EUR_plink_2Mb) <- c("Sample","IID","PHE","nRoH","KB","KBAVG","FRoH","Length")
RoH_W_EUR_plink_2Mb <- merge(RoH_W_EUR_plink_2Mb, Lizards_admix, by.x = "Sample", by.y = "ID")
p_pl_W_EUR_2Mb_1 <- ggplot(RoH_W_EUR_plink_2Mb, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  geom_text_repel(data = subset(RoH_W_EUR_plink_2Mb, Origin == "Int-W-EUR"), aes(label = ID_alias), color="black") +
  theme_bw()
p_pl_W_EUR_2Mb_2 <- ggplot(RoH_W_EUR_plink_2Mb, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0,0.18) + theme_bw()

RoH_C_ITA_plink_500kb <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_plink_C-ITA_500kb.hom.indiv", header = T)
RoH_C_ITA_plink_500kb$FRoH <- RoH_C_ITA_plink_500kb$KB/RoH_C_ITA_plink_500kb$KB[RoH_C_ITA_plink_500kb$FID == "Pxx"]
RoH_C_ITA_plink_500kb <- subset(RoH_C_ITA_plink_500kb, FID != "Pxx")
RoH_C_ITA_plink_500kb$Length <- RoH_C_ITA_plink_500kb$KB/1000
colnames(RoH_C_ITA_plink_500kb) <- c("Sample","IID","PHE","nRoH","KB","KBAVG","FRoH","Length")
RoH_C_ITA_plink_500kb <- merge(RoH_C_ITA_plink_500kb, Lizards_admix, by.x = "Sample", by.y = "ID")
p_pl_C_ITA_500kb_1 <- ggplot(RoH_C_ITA_plink_500kb, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  geom_text_repel(data = subset(RoH_C_ITA_plink_500kb, Origin == "Int-C-ITA"), aes(label = ID_alias), color="black") +
  theme_bw()
p_pl_C_ITA_500kb_2 <- ggplot(RoH_C_ITA_plink_500kb, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0,0.35) + theme_bw()

RoH_W_EUR_plink_500kb <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_plink_W-EUR_500kb.hom.indiv", header = T)
RoH_W_EUR_plink_500kb$FRoH <- RoH_W_EUR_plink_500kb$KB/RoH_W_EUR_plink_500kb$KB[RoH_W_EUR_plink_500kb$FID == "Pxx"]
RoH_W_EUR_plink_500kb <- subset(RoH_W_EUR_plink_500kb, FID != "Pxx")
RoH_W_EUR_plink_500kb$Length <- RoH_W_EUR_plink_500kb$KB/1000
colnames(RoH_W_EUR_plink_500kb) <- c("Sample","IID","PHE","nRoH","KB","KBAVG","FRoH","Length")
RoH_W_EUR_plink_500kb <- merge(RoH_W_EUR_plink_500kb, Lizards_admix, by.x = "Sample", by.y = "ID")
p_pl_W_EUR_500kb_1 <- ggplot(RoH_W_EUR_plink_500kb, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  geom_text_repel(data = subset(RoH_W_EUR_plink_500kb, Origin == "Int-W-EUR"), aes(label = ID_alias), color="black") +
  theme_bw()
p_pl_W_EUR_500kb_2 <- ggplot(RoH_W_EUR_plink_500kb, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0,0.35) + theme_bw()

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/HetRoH_Supp_V2.pdf", height=10, width=14, useDingbats = F)
(p1S | p2S | p3S | p4S) /
(p_pl_C_ITA_2Mb_1 | p_pl_C_ITA_2Mb_2 | p_pl_W_EUR_2Mb_1 | p_pl_W_EUR_2Mb_2) /
(p_pl_C_ITA_500kb_1 | p_pl_C_ITA_500kb_2 | p_pl_W_EUR_500kb_1 | p_pl_W_EUR_500kb_2)
dev.off()

RoH_C_ITA_plink_2Mb %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_W_EUR_plink_2Mb %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_C_ITA_Supp %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_W_EUR_Supp %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_C_ITA_plink_500kb %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_W_EUR_plink_500kb %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())


#### Position on chromosome

# 2.2 ROHs ----
# 2.2.1 BCF tools ----
C_ITA_BCF <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/roh.pseudo.qual.C-ITA.2Mb", h=F)

# Naming the columns
C_ITA_BCF <- C_ITA_BCF %>%
  dplyr::select(-1) %>% # Remove the first column - a bunch of RG
  rename(ID = V2, CHR = V3, POS1 = V4, POS2 = V5, BP = V6, MARKERS = V7, QUALITY = V8) 

# Merge the clean data with the origin
C_ITA_BCF<- merge(Lizards_admix,C_ITA_BCF,by="ID", all.x=T)
C_ITA_BCF <- C_ITA_BCF %>%
  filter(!grepl("W-EUR", Origin))

# Define names  of chromosome names to  numbers
C_ITA_BCF <- C_ITA_BCF %>%
  mutate(CHR = case_when(
    CHR == "CM014743.1" ~ "1",
    CHR == "CM014744.1" ~ "2",
    CHR == "CM014745.1" ~ "3",
    CHR == "CM014746.1" ~ "4",
    CHR == "CM014747.1" ~ "5",
    CHR == "CM014748.1" ~ "6",
    CHR == "CM014749.1" ~ "7",
    CHR == "CM014750.1" ~ "8",
    CHR == "CM014751.1" ~ "9",
    CHR == "CM014752.1" ~ "10",
    CHR == "CM014753.1" ~ "11",
    CHR == "CM014754.1" ~ "12",
    CHR == "CM014755.1" ~ "13",
    CHR == "CM014756.1" ~ "14",
    CHR == "CM014757.1" ~ "15",
    CHR == "CM014758.1" ~ "16",
    CHR == "CM014759.1" ~ "17",
    CHR == "CM014760.1" ~ "18",))

# W-EUR 
W_EUR_BCF <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/roh.pseudo.qual.W-EUR.2Mb",h=F)

# Naming the columns and filtering for quality (Phred score) and minimum length (500k)
W_EUR_BCF <- W_EUR_BCF %>%
  dplyr::select(-1) %>% # Remove the first column - a bunch of RG
  rename(ID = V2, CHR = V3, POS1 = V4, POS2 = V5, BP = V6, MARKERS = V7, QUALITY = V8) 

# Merge the clean data with the origin
W_EUR_BCF<- merge(Lizards_admix,W_EUR_BCF,by="ID", all.x=T)
W_EUR_BCF <- W_EUR_BCF %>%
  filter(!grepl("C-ITA", Origin))

# Change chromosome names 
W_EUR_BCF <- W_EUR_BCF %>%
  mutate(CHR = case_when(
    CHR == "CM014743.1" ~ "1",
    CHR == "CM014744.1" ~ "2",
    CHR == "CM014745.1" ~ "3",
    CHR == "CM014746.1" ~ "4",
    CHR == "CM014747.1" ~ "5",
    CHR == "CM014748.1" ~ "6",
    CHR == "CM014749.1" ~ "7",
    CHR == "CM014750.1" ~ "8",
    CHR == "CM014751.1" ~ "9",
    CHR == "CM014752.1" ~ "10",
    CHR == "CM014753.1" ~ "11",
    CHR == "CM014754.1" ~ "12",
    CHR == "CM014755.1" ~ "13",
    CHR == "CM014756.1" ~ "14",
    CHR == "CM014757.1" ~ "15",
    CHR == "CM014758.1" ~ "16",
    CHR == "CM014759.1" ~ "17",
    CHR == "CM014760.1" ~ "18",)) 


# Positions in chr 1 
C_ITA_BCF <- C_ITA_BCF %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))
W_EUR_BCF <- W_EUR_BCF %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))

# Filter for Chromosome 1
C_ITA_BCF_chr1 <- C_ITA_BCF %>% filter(CHR == 1)

# Get all unique IDs
all_IDs <- unique(C_ITA_BCF$ID_alias)

# Create a lookup table for Origin and Abbpop per ID
ID_lookup <- C_ITA_BCF %>%
  dplyr::select(ID_alias, Origin, Abbpop) %>%
  distinct()

# Fill missing IDs for Chr1, keeping Origin and Abbpop
C_ITA_BCF_chr1_filled <- C_ITA_BCF_chr1 %>%
  right_join(
    tibble(ID_alias = all_IDs),
    by = "ID_alias"
  ) %>%
  # add Origin/Abbpop from lookup if NA
  left_join(ID_lookup, by = "ID_alias", suffix = c("", ".lookup")) %>%
  mutate(
    Origin = coalesce(Origin, Origin.lookup),
    Abbpop = coalesce(Abbpop, Abbpop.lookup),
    POS1 = ifelse(is.na(POS1), 0, POS1),
    POS2 = ifelse(is.na(POS2), 0, POS2),
    CHR  = ifelse(is.na(CHR), 0, CHR),
    BP = ifelse(is.na(BP), 0, BP),
    MARKERS = ifelse(is.na(MARKERS), 0, MARKERS),
    QUALITY = ifelse(is.na(QUALITY), 0, QUALITY)
  ) %>%
  dplyr::select(-Origin.lookup, -Abbpop.lookup)  # remove temporary columns

# Reorder ID by Origin
C_ITA_BCF_chr1_filled <- C_ITA_BCF_chr1_filled %>%
  mutate(
    ID_alias = factor(ID_alias, levels = unique(ID_alias[order(Origin)])))

# Plot only Chr 1
chr1_C_ITA <- ggplot(C_ITA_BCF_chr1_filled, aes(x=POS1, xend=POS2, y=ID_alias, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID_alias), linewidth =3) +  
  scale_color_manual(values= c("Nat-C-ITA" = "#0673B3", "Int-C-ITA" = "#6BCBDA")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), limits = c(0, 130727322)) +
  theme_minimal() + 
  labs(x="Genomic Position", y="Sample", 
       color="Origin") +
  theme(strip.text = element_text(size=12), 
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))

### same for France

# Filter for Chromosome 1
W_EUR_BCF_chr1 <- W_EUR_BCF %>% filter(CHR == 1)

# Get all unique IDs
all_IDs <- unique(W_EUR_BCF$ID_alias)

# Create a lookup table for Origin and Abbpop per ID
ID_lookup <- W_EUR_BCF %>%
  dplyr::select(ID_alias, Origin, Abbpop) %>%
  distinct()

# Fill missing IDs for Chr1, keeping Origin and Abbpop
W_EUR_BCF_chr1_filled <- W_EUR_BCF_chr1 %>%
  right_join(
    tibble(ID_alias = all_IDs),
    by = "ID_alias"
  ) %>%
  # add Origin/Abbpop from lookup if NA
  left_join(ID_lookup, by = "ID_alias", suffix = c("", ".lookup")) %>%
  mutate(
    Origin = coalesce(Origin, Origin.lookup),
    Abbpop = coalesce(Abbpop, Abbpop.lookup),
    POS1 = ifelse(is.na(POS1), 0, POS1),
    POS2 = ifelse(is.na(POS2), 0, POS2),
    CHR  = ifelse(is.na(CHR), 0, CHR),
    BP = ifelse(is.na(BP), 0, BP),
    MARKERS = ifelse(is.na(MARKERS), 0, MARKERS),
    QUALITY = ifelse(is.na(QUALITY), 0, QUALITY)
  ) %>%
  dplyr::select(-Origin.lookup, -Abbpop.lookup)  # remove temporary columns

# Reorder ID by Origin
W_EUR_BCF_chr1_filled <- W_EUR_BCF_chr1_filled %>%
  mutate(ID_alias = factor(ID_alias, levels = unique(ID_alias[order(Origin)])))

# Plot only Chr 1
chr1_W_EUR <- ggplot(W_EUR_BCF_chr1_filled, aes(x=POS1, xend=POS2, y=ID_alias, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID_alias), linewidth =3) +  
  scale_color_manual(values= c("Nat-W-EUR" = "#DD6A27", "Int-W-EUR" = "#EFE808")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), limits = c(0, 130727322)) +
  theme_minimal() + 
  labs(x="Genomic Position", y="Sample", 
       color="Origin") +
  theme(strip.text = element_text(size=12), 
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))


pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Chr1_V2.pdf", height=6, width=10, useDingbats = F)
chr1_C_ITA/chr1_W_EUR
dev.off()

###########################
#### IDrisk
###########################

Other_IDRisk <- readWorkbook("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/Kyriazis2025_TREE_RoHs_Supplements.xlsx", sheet = 1)
Other_IDRisk <- Other_IDRisk[c(2,6,7,11,24),c(1,5)]
colnames(Other_IDRisk) <- c("Sample","IDRisk")
Other_IDRisk$Origin <- "Lit"

RoH_Het_C_ITA <- RoH_Het_C_ITA %>%
  arrange(desc(IDRisk)) %>%  # sort by IDRisk descending
  mutate(Sample = factor(ID_alias, levels = ID_alias))

RoH_Het_W_EUR <- RoH_Het_W_EUR %>%
  arrange(desc(IDRisk)) %>%  # sort by IDRisk descending
  mutate(Sample = factor(ID_alias, levels = ID_alias))

common_cols <- intersect(colnames(RoH_Het_C_ITA), colnames(Other_IDRisk))

df_combined <- rbind(
  RoH_Het_C_ITA[common_cols],
  RoH_Het_W_EUR[common_cols],
  Other_IDRisk[common_cols]
)

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/IDRisk_V2.pdf", height=6, width=6, useDingbats = F)
ggplot(df_combined, aes(x = Sample, y = IDRisk, fill = Origin)) +
  geom_col() +
  scale_fill_manual(values = c("Int-C-ITA" = "#6BCBDA", "Nat-C-ITA" = "#0673B3","Nat-W-EUR" = "#DD6A27", "Int-W-EUR" = "#EFE808")) +
  theme_bw() +
  labs(x = "ID Risk", y = "Value", title = "Barplot of ID Risk by Origin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7))
dev.off()

###########################
# 3. PURGING  -----
###########################
# 3.1 Calculate Rxy for C-ITA-origin samples ----

# Load C-ITA Purging dataset - Relative frequencies and Jackknifing information
C_ITA_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Central-Italy_RefSeq_freq.tsv", h=T)

# Calculate Lxy and Lyx ratios
C_ITA_Purging <- C_ITA_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH / Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH / Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE / Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE / Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW / Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW / Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER / Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER / Sum_fyx_INTERGENIC,
         
         # Compute Rxy per mutation impact category
         Rxy_High = Lxy_HIGH / Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE / Lyx_MODERATE,
         Rxy_Low = Lxy_LOW / Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER / Lyx_MODIFIER)

# Select only relevant Rxy values
Rxy_C_ITA <- C_ITA_Purging %>%
  dplyr::select(Rxy_High, Rxy_Moderate, Rxy_Low, Rxy_Modifier)

# Reshape data into long format for analysis
Rxy_C_ITA <- Rxy_C_ITA %>% 
  pivot_longer(cols = starts_with("Rxy_"),   # Select Rxy columns
               names_to = "Impact",          # Create column for categories
               values_to = "Rxy") %>%        # Assign values to Rxy
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean category names

# Function to retrieve summary of the results after jackknifing 
calculate_minmax <- function(data) { 
  data %>%
    group_by(Impact) %>% 
    summarise(
      mean_Rxy = mean(Rxy, na.rm = TRUE), #Retrieve mean value 
      min_Rxy = min(Rxy, na.rm = TRUE), #Retrieve minimum value
      max_Rxy = max(Rxy, na.rm = TRUE), #Retrieve maximum value 
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), 
                  ~ formatC(., format = "f", digits = 3)))}

calculate_minmax(Rxy_C_ITA)

# Visualize C-ITA-origin 
# Boxplot
C_ITA_RefSeq <- ggplot(Rxy_C_ITA, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#0673B3",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle =  "C-ITA",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.97, 1.04))


# 3.2 Calculate Rxy for W-EUR-origin samples ----
W_EUR_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_West-Europe_RefSeq_freq.tsv", h=T)

# Get Lxy and Lyx 
W_EUR_Purging <- W_EUR_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH/Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH/Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE/Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE/Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW/Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW/Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER/Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER/Sum_fyx_INTERGENIC,
         
         # Get Rxy per impact 
         Rxy_High = Lxy_HIGH/Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE/Lyx_MODERATE,
         Rxy_Low = Lxy_LOW/Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER/Lyx_MODIFIER)

# Filter just Rxy values
Rxy_W_EUR<- W_EUR_Purging %>%
  dplyr::select(Rxy_High,
                Rxy_Moderate,
                Rxy_Low,
                Rxy_Modifier)

# Transform data 
Rxy_W_EUR<- Rxy_W_EUR %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Visualize W-EUR-origin 
# BoxPlot 
W_EUR_RefSeq <- ggplot(Rxy_W_EUR, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#DD6A27",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle = "W-EUR",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.77, 1.04))

qq<- C_ITA_RefSeq + W_EUR_RefSeq
print(qq)  



# Confirm data as factors and check min max values after jackknifing
calculate_minmax(Rxy_W_EUR)



# 3.3 Calculate Rxy for each population ----

# Define populations
c_ita_pops <- c("BB", "DL", "NF", "SH", "VT", "WS")
w_eur_pops <- c("BU", "WB", "WE")

# List input files
freq_files <- list.files("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/IndividualPopulations", pattern = "*_freq\\.tsv$", full.names = TRUE)

# Function to calculate Rxy per file
calculate_Rxy <- function(file_path) {
  pop <- stringr::str_extract(basename(file_path), ".*(?=_freq.tsv)")
  df <- read.table(file_path, header = TRUE)
  
  df <- df %>%
    mutate(Population = pop,
           Lxy_HIGH = Sum_fxy_HIGH / Sum_fxy_INTERGENIC,
           Lyx_HIGH = Sum_fyx_HIGH / Sum_fyx_INTERGENIC,
           Lxy_MODERATE = Sum_fxy_MODERATE / Sum_fxy_INTERGENIC,
           Lyx_MODERATE = Sum_fyx_MODERATE / Sum_fyx_INTERGENIC,
           Lxy_LOW = Sum_fxy_LOW / Sum_fxy_INTERGENIC,
           Lyx_LOW = Sum_fyx_LOW / Sum_fyx_INTERGENIC,
           Lxy_MODIFIER = Sum_fxy_MODIFIER / Sum_fxy_INTERGENIC,
           Lyx_MODIFIER = Sum_fyx_MODIFIER / Sum_fyx_INTERGENIC,
           Rxy_High = Lxy_HIGH / Lyx_HIGH,
           Rxy_Moderate = Lxy_MODERATE / Lyx_MODERATE,
           Rxy_Low = Lxy_LOW / Lyx_LOW,
           Rxy_Modifier = Lxy_MODIFIER / Lyx_MODIFIER) %>%
    dplyr::select(Population, starts_with("Rxy_")) %>%
    pivot_longer(cols = starts_with("Rxy_"),
                 names_to = "Impact", values_to = "Rxy") %>%
    mutate(Impact = stringr::str_remove(Impact, "Rxy_"))
  
  return(df)
}

# Process all files
Rxy_all <- bind_rows(lapply(freq_files, calculate_Rxy)) %>%
  mutate(Group = case_when(
    Population %in% c_ita_pops ~ "C-ITA",
    Population %in% w_eur_pops ~ "W-EUR",
    TRUE ~ "Native" ))

# Gather all the results 
population_results <- lapply(unique(Rxy_all$Population), function(pop) {
  Rxy_all %>%
    filter(Population == pop) %>%
    calculate_minmax() %>%
    mutate(Population = pop) %>%
    dplyr::select(Population, everything())
})

# Name results properly
names(population_results) <- unique(Rxy_all$Population)

# Print formatted results
for (pop in names(population_results)) {
  cat("\n=== Results for", pop, "===\n")
  print(as_tibble(population_results[[pop]]))
}

# Generate plots per population
# Create C-ITA plot with green color
boxplot_C_ITA <- Rxy_all %>%
  filter(Group == "C-ITA") %>%
  ggplot(aes(x = Rxy, y = Impact)) +
  geom_boxplot(fill = "#0673B3",  # Green for C-ITA
               outlier.shape = NA, 
               alpha = 0.6) +
  facet_wrap(~ Population, nrow = 2) +
  theme_bw() +
  scale_x_continuous(limits = c(0.75, 1.1)) +
  labs(x = "Rxy", y = "Impact Category") +
  theme(legend.position = "none")

# Create W-EUR plot with orange color
boxplot_W_EUR <- Rxy_all %>%
  filter(Group == "W-EUR") %>%
  ggplot(aes(x = Rxy, y = Impact)) +
  geom_boxplot(fill = "#DD6A27",  # Orange for W-EUR
               outlier.shape = NA, 
               alpha = 0.6) +
  facet_wrap(~ Population, nrow = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(0.75, 1.1)) +
  labs(x = "Rxy", y = "Impact Category") +
  theme(legend.position = "none")


qq<- (C_ITA_RefSeq + W_EUR_RefSeq) / boxplot_C_ITA / boxplot_W_EUR + plot_layout(heights = c(1, 1, 0.5))
pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Purging_V1.pdf", height=12, width=10, useDingbats = F)
print(qq)
dev.off()


#### Purging for Supplementary Material

# Load C-ITA Purging dataset - Relative frequencies and Jackknifing information
C_ITA_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Central-Italy_Ensemble_freq.tsv", h=T)

# Calculate Lxy and Lyx ratios
C_ITA_Purging <- C_ITA_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH / Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH / Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE / Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE / Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW / Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW / Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER / Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER / Sum_fyx_INTERGENIC,
         
         # Compute Rxy per mutation impact category
         Rxy_High = Lxy_HIGH / Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE / Lyx_MODERATE,
         Rxy_Low = Lxy_LOW / Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER / Lyx_MODIFIER)

# Select only relevant Rxy values
Rxy_C_ITA <- C_ITA_Purging %>%
  dplyr::select(Rxy_High, Rxy_Moderate, Rxy_Low, Rxy_Modifier)

# Reshape data into long format for analysis
Rxy_C_ITA <- Rxy_C_ITA %>% 
  pivot_longer(cols = starts_with("Rxy_"),   # Select Rxy columns
               names_to = "Impact",          # Create column for categories
               values_to = "Rxy") %>%        # Assign values to Rxy
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean category names

calculate_minmax(Rxy_C_ITA)

# Visualize C-ITA-origin 
# Boxplot
C_ITA_Ensemble <- ggplot(Rxy_C_ITA, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#0673B3",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle =  "C-ITA",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.90, 1.04))

# Load C-ITA Purging dataset - Relative frequencies and Jackknifing information
C_ITA_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Central-Italy_Tiberius_freq.tsv", h=T)

# Calculate Lxy and Lyx ratios
C_ITA_Purging <- C_ITA_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH / Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH / Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE / Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE / Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW / Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW / Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER / Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER / Sum_fyx_INTERGENIC,
         
         # Compute Rxy per mutation impact category
         Rxy_High = Lxy_HIGH / Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE / Lyx_MODERATE,
         Rxy_Low = Lxy_LOW / Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER / Lyx_MODIFIER)

# Select only relevant Rxy values
Rxy_C_ITA <- C_ITA_Purging %>%
  dplyr::select(Rxy_High, Rxy_Moderate, Rxy_Low, Rxy_Modifier)

# Reshape data into long format for analysis
Rxy_C_ITA <- Rxy_C_ITA %>% 
  pivot_longer(cols = starts_with("Rxy_"),   # Select Rxy columns
               names_to = "Impact",          # Create column for categories
               values_to = "Rxy") %>%        # Assign values to Rxy
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean category names

calculate_minmax(Rxy_C_ITA)

# Visualize C-ITA-origin 
# Boxplot
C_ITA_Tiberius <- ggplot(Rxy_C_ITA, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#0673B3",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle =  "C-ITA",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.90, 1.04))

# 3.2 Calculate Rxy for W-EUR-origin samples ----
W_EUR_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_West-Europe_Ensemble_freq.tsv", h=T)

# Get Lxy and Lyx 
W_EUR_Purging <- W_EUR_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH/Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH/Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE/Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE/Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW/Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW/Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER/Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER/Sum_fyx_INTERGENIC,
         
         # Get Rxy per impact 
         Rxy_High = Lxy_HIGH/Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE/Lyx_MODERATE,
         Rxy_Low = Lxy_LOW/Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER/Lyx_MODIFIER)

# Filter just Rxy values
Rxy_W_EUR<- W_EUR_Purging %>%
  dplyr::select(Rxy_High,
                Rxy_Moderate,
                Rxy_Low,
                Rxy_Modifier)

# Transform data 
Rxy_W_EUR<- Rxy_W_EUR %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Visualize W-EUR-origin 
# BoxPlot 
W_EUR_Ensemble <- ggplot(Rxy_W_EUR, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#DD6A27",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle = "W-EUR",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.76, 1.04))

W_EUR_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_West-Europe_Tiberius_freq.tsv", h=T)

# Get Lxy and Lyx 
W_EUR_Purging <- W_EUR_Purging %>%
  mutate(Lxy_HIGH = Sum_fxy_HIGH/Sum_fxy_INTERGENIC,
         Lyx_HIGH = Sum_fyx_HIGH/Sum_fyx_INTERGENIC,
         Lxy_MODERATE = Sum_fxy_MODERATE/Sum_fxy_INTERGENIC,
         Lyx_MODERATE = Sum_fyx_MODERATE/Sum_fyx_INTERGENIC,
         Lxy_LOW = Sum_fxy_LOW/Sum_fxy_INTERGENIC,
         Lyx_LOW = Sum_fyx_LOW/Sum_fyx_INTERGENIC,
         Lxy_MODIFIER = Sum_fxy_MODIFIER/Sum_fxy_INTERGENIC,
         Lyx_MODIFIER = Sum_fyx_MODIFIER/Sum_fyx_INTERGENIC,
         
         # Get Rxy per impact 
         Rxy_High = Lxy_HIGH/Lyx_HIGH,
         Rxy_Moderate = Lxy_MODERATE/Lyx_MODERATE,
         Rxy_Low = Lxy_LOW/Lyx_LOW,
         Rxy_Modifier = Lxy_MODIFIER/Lyx_MODIFIER)

# Filter just Rxy values
Rxy_W_EUR<- W_EUR_Purging %>%
  dplyr::select(Rxy_High,
                Rxy_Moderate,
                Rxy_Low,
                Rxy_Modifier)

# Transform data 
Rxy_W_EUR<- Rxy_W_EUR %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Visualize W-EUR-origin 
# BoxPlot 
W_EUR_Tiberius <- ggplot(Rxy_W_EUR, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#DD6A27",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle = "W-EUR",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.76, 1.04))

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Purging_Supp.pdf", height=8, width=10, useDingbats = F)
(C_ITA_Ensemble + W_EUR_Ensemble) / (C_ITA_Tiberius + W_EUR_Tiberius)
dev.off()



#####################################################
### Plot the position of Recombination HotSpots
#####################################################

# ---- 1. Read chromosome sizes ----
chrom <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Chr_Length_GenBank.txt", h=F)
chrom$start <- 1
chrom <- chrom[,c(1,3,2)]
colnames(chrom) <- c("chr", "start", "end")

# ---- 2. Read VCF (skip header lines starting with ##) ----
#HotSpots <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/RecombHotSpots/All_Merged_Hotspots.bed", header = T)
HotSpots <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/RecombHotSpots/allchrom_w2kb_s1kb_f40kb_fd5_sw0_sf0_overlaps-mergeall.txt", header = T)

# ---- 3. Plot ----
pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/RecombHotspots.pdf", height=5, width=13, useDingbats = F)
ggplot() + geom_segment(data = chrom, aes(x = start, xend = end, y = chr, yend = chr), size = 3, color = "grey70") + 
  geom_point(data = HotSpots, aes(x = start, y = chr),  alpha=0.4, color = "red", shape=16) +
  theme_bw() + ggtitle("Recombination hotspots (N=3297)") + labs(x = "Genomic position", y = "Chromosome") + theme(panel.grid = element_blank())
dev.off()


##########################
### Position of HIGH impact variants
##########################

# ---- 1. Read chromosome sizes ----
chrom <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Chromosome_Length.txt", h=F)
colnames(chrom) <- c("chr", "start", "end")

# ---- 2. Read VCF (skip header lines starting with ##) ----
vcf <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Polarised_high_RefSeq.vcf.gz", comment.char = "#", header = FALSE)

# Extract relevant columns
colnames(vcf)[1:2] <- c("chr", "pos")
vcf <- vcf[, c("chr", "pos")]

# ---- 3. Plot ----
p_RefSeq <- ggplot() + geom_segment(data = chrom, aes(x = start, xend = end, y = chr, yend = chr), size = 3, color = "grey70") + 
  geom_point(data = vcf, aes(x = pos, y = chr), size = 2, alpha = 0.7, color = "red") +
  theme_bw() + ggtitle("RefSeq (N=382)") + labs(x = "Genomic position", y = "Chromosome") + theme(panel.grid = element_blank())


# ---- 2. Read VCF (skip header lines starting with ##) ----
vcf <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Polarised_high_Tiberius.vcf.gz", comment.char = "#", header = FALSE)

# Extract relevant columns
colnames(vcf)[1:2] <- c("chr", "pos")
vcf <- vcf[, c("chr", "pos")]

# ---- 3. Plot ----
p_Tiberius <- ggplot() + geom_segment(data = chrom, aes(x = start, xend = end, y = chr, yend = chr), size = 3, color = "grey70") + 
  geom_point(data = vcf, aes(x = pos, y = chr), size = 2, alpha = 0.7, color = "red") +
  theme_bw() + ggtitle("Tiberius (N=287)") + labs(x = "Genomic position", y = "Chromosome") + theme(panel.grid = element_blank())


# ---- 2. Read VCF (skip header lines starting with ##) ----
vcf <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Polarised_high_Ensemble.vcf.gz", comment.char = "#", header = FALSE)

# Extract relevant columns
colnames(vcf)[1:2] <- c("chr", "pos")
vcf <- vcf[, c("chr", "pos")]

vcf <- vcf %>%
  mutate(chr = case_when(
    chr == "1" ~ "NC_041312.1",
    chr == "2" ~ "NC_041313.1",
    chr == "3" ~ "NC_041314.1",
    chr == "4" ~ "NC_041315.1",
    chr == "5" ~ "NC_041316.1",
    chr == "6" ~ "NC_041317.1",
    chr == "7" ~ "NC_041318.1",
    chr == "8" ~ "NC_041319.1",
    chr == "9" ~ "NC_041320.1",
    chr == "10" ~ "NC_041321.1",
    chr == "11" ~ "NC_041322.1",
    chr == "12" ~ "NC_041323.1",
    chr == "13" ~ "NC_041324.1",
    chr == "14" ~ "NC_041325.1",
    chr == "15" ~ "NC_041326.1",
    chr == "16" ~ "NC_041327.1",
    chr == "17" ~ "NC_041328.1",
    chr == "18" ~ "NC_041329.1",)) 

# ---- 3. Plot ----
p_Ensemble <- ggplot() + geom_segment(data = chrom, aes(x = start, xend = end, y = chr, yend = chr), size = 3, color = "grey70") + 
  geom_point(data = vcf, aes(x = pos, y = chr), size = 2, alpha = 0.7, color = "red") +
  theme_bw() + ggtitle("Ensemble (N=697)") + labs(x = "Genomic position", y = "Chromosome") + theme(panel.grid = element_blank())

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/PositionHigh.pdf", height=6, width=16, useDingbats = F)
print(p_RefSeq+p_Tiberius+ p_Ensemble)
dev.off()



## Genetic load (RefSeq High impact variants; N=382)

###
### Genetic load (i.e., High impact variants)
High_alleles <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/RefSeq_genotypes_high.tsv", h=T)
sample_cols <- colnames(High_alleles)[3:ncol(High_alleles)]

# Initialize a results data frame
results <- data.frame(
  Sample = sample_cols,
  Hom_Ancestral = 0,
  Heterozygous = 0,
  Hom_Derived = 0,
  stringsAsFactors = FALSE
)

# Function to classify genotypes
classify_gt <- function(gt) {
  if (gt %in% c("0/0","0|0")) {
    return("hom_ancestral")
  } else if (gt %in% c("0/1","1/0","0|1","1|0")) {
    return("heterozygous")
  } else if (gt %in% c("1/1","1|1")) {
    return("hom_derived")
  } else {
    return(NA)  # for missing or other genotypes
  }
}

# Count genotypes per sample
for (s in sample_cols) {
  gt_class <- sapply(High_alleles[[s]], classify_gt)
  results$Hom_Ancestral[results$Sample == s] <- sum(gt_class == "hom_ancestral", na.rm = TRUE)
  results$Heterozygous[results$Sample == s]   <- sum(gt_class == "heterozygous", na.rm = TRUE)
  results$Hom_Derived[results$Sample == s]   <- sum(gt_class == "hom_derived", na.rm = TRUE)
}

# View the summary
tt(results)
write.csv(results, "C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/Load_per_individual.csv")
Lizards <- readWorkbook("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Samples_Info.xlsx", sheet = 2)
results <- merge(results, Lizards, by.x="Sample", by.y="ID")

# Summarize by Origin
summary_origin <- results %>%
  group_by(Origin) %>%
  summarise(
    mean_Hom_Ancestral = mean(Hom_Ancestral),
    mean_Heterozygous = mean(Heterozygous),
    mean_Hom_Derived  = mean(Hom_Derived),
    .groups = 'drop'
  )

tt(summary_origin)
write.csv(summary_origin, "C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/Load_per_origin.csv")

# Summarize by Abbpop
summary_abb <- results %>%
  group_by(Abbpop) %>%
  summarise(
    mean_Hom_Ancestral = mean(Hom_Ancestral),
    mean_Heterozygous = mean(Heterozygous),
    mean_Hom_Derived  = mean(Hom_Derived),
    .groups = 'drop'
  )

tt(summary_abb)
write.csv(summary_abb, "C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/Load_per_abbpop.csv")
