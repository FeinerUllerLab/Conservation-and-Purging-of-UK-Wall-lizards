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

#France + Italy
#world <- ne_countries(scale = "medium", returnclass = "sf")
#france_italy <- world %>% filter(admin %in% c("France", "Italy","Spain","Germany","Austria","Switzerland","Croatia"))

map_native <- ggplot(data = europe_cropped_native) +
  geom_sf(fill = "white", color = "grey", size=0.1) +
  geom_sf(data = st_crop(Pmur_dist, xmin = -5, xmax = 14, ymin = 41, ymax = 49), fill = "lightgrey", color = NA) +
  geom_point(data = NativeLizards, aes(x = Longitude, y = Latitude, colour = Origin),
             size = 5, shape = 20) +
  coord_sf(xlim = c(-5, 14), ylim = c(41, 49), expand = FALSE) +
  scale_colour_manual(values = c("FRA" = "#DD6A27", "ITA" = "#0673B3")) +
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
#UK<- ne_countries(country = "united kingdom", returnclass = "sf", scale="large")

map_UK <- ggplot(data = europe_cropped_intro) +
  geom_sf(fill = "white", color = "grey", size = 0.1) +
  geom_point(data = UKLizards, aes(x = Longitude, y = Latitude, colour = Origin),
             size = 5, shape = 20) +
  coord_sf(xlim = c(-6.0, 2.0), ylim = c(49.5, 52.0), expand = FALSE) +
  scale_colour_manual(values = c("FRA" = "#EFE808", "ITA" = "#6BCBDA")) +
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
    values = c("ITA-Native" = "#0673B3", 
               "FRA-Native" = "#DD6A27", 
               "ITA-Intro" = "#6BCBDA", 
               "FRA-Intro" = "#EFE808"),
    labels = c("ITA-Native" = "Native Italian", 
               "FRA-Native" = "Native French", 
               "ITA-Intro" = "Non-native Italian", 
               "FRA-Intro" = "Non-native French")) +
  labs(subtitle = "A") +
  theme_bw() +
  theme(legend.position = "bottom", axis.title = element_text(size = 12))

# Create p2 without legend but updated labels
p2 <- ggplot(pca_data, aes(x = PC1, y = PC3, color = Group)) +
  geom_point(size = 6, alpha=0.7, shape=16) +
  labs(x = paste("PC1 (", evec1.pc, "%)", sep = ""),
       y = paste("PC3 (", evec3.pc, "%)", sep = "")) +
  scale_color_manual(
    values = c("ITA-Native" = "#0673B3", 
               "FRA-Native" = "#DD6A27", 
               "ITA-Intro" = "#6BCBDA", 
               "FRA-Intro" = "#EFE808"),
    labels = c("ITA-Native" = "Native Italian", 
               "FRA-Native" = "Native French", 
               "ITA-Intro" = "Non-native Italian", 
               "FRA-Intro" = "Non-native French")) +
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
bootstrap_data <- data.frame(node = internal_nodes,
                             bootstrap = as.numeric(tree$node.label))# Convert to numeric

# Merge with main tree data
tree_data <- tree_data %>% 
  left_join(bootstrap_data, by = "node") %>% 
  # Merge with lizard metadata (tips only)
  left_join(LizardsData, by = c("label" = "ID"))

# Create the tree plot
treeplot <- ggtree(tree, layout = "unrooted") %<+% tree_data +
  geom_tiplab2(aes(color = Group), size = 3, offset = 0.0, align = FALSE, show.legend = FALSE) +
  geom_tippoint(aes(color = Group), size = 2) +
  geom_nodelab(aes(label = round(bootstrap)), 
               color = "black", 
               hjust = -0.05,
               size = 2, 
               na.rm = TRUE,) +
  scale_color_manual(values = c("ITA-Native" = "#0673B3", 
                                "FRA-Native" = "#DD6A27", 
                                "ITA-Intro" = "#6BCBDA", 
                                "FRA-Intro" = "#EFE808"),
                     labels = c("Non-native France", "Non-native Italy", 
                                "Native France", "Native Italy")) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right", panel.grid.major = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)  )

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Tree.pdf", height=10, width=10, useDingbats = F)
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
colnames(K2) <- c("ID","Origin","Abbpop","Q1","Q2")
K3<- cbind(Lizards_admix, K3)
colnames(K3) <- c("ID","Origin","Abbpop","Q1","Q2","Q3")
K4<- cbind(Lizards_admix, K4)
colnames(K4) <- c("ID","Origin","Abbpop","Q1","Q2","Q3","Q4")

# Organize samples based on K2 (Best supported)
ordered_individuals <- K2 %>%
  arrange(-Q1) %>%
  pull(ID)  

# Apply ordering before pivoting
K2$ID <- factor(K2$ID, levels = ordered_individuals)
K3$ID <- factor(K3$ID, levels = ordered_individuals)
K4$ID <- factor(K4$ID, levels = ordered_individuals)

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
admix_data$ID <- factor(admix_data$ID, levels = ordered_individuals)


# 1.3.3 Plot the ADMIXTURE ----
pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Admixture.pdf", height=6, width=10, useDingbats = F)
ggplot(admix_data, aes(x = ID, y = Q_value, fill = Ancestry)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~K, ncol = 1) +
  scale_fill_manual(values = as.vector(paletteer_d("ggsci::default_jco"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Admixture Proportion", x = "Individual")
dev.off()


#### 2. GENETIC DIVERSITY AND INBREEDING ####

## 2.1. HETEROZYGOSITY + ROH ----

RoH_IT <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_bcftools_IT_2Mb.txt", header = T)
RoH_IT <- merge(RoH_IT, Lizards, by.x = "Sample", by.y = "ID")

geno_Het_IT <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/heterozygosity_summary_genomewide_IT.tsv", header = T)
geno_Het_IT$O.HET <- geno_Het_IT$NumberVariableSites-geno_Het_IT$ObservedHomozygous
geno_Het_IT$geno_F.HET <- geno_Het_IT$O.HET/(geno_Het_IT$GenotypedSites-geno_Het_IT$MissingSites)
geno_Het_IT <- geno_Het_IT[,-(2:6)]
RoH_HET_IT <- merge(RoH_IT, geno_Het_IT, by = "Sample")

nonRoH_Het_IT <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/heterozygosity_summary_nonRoH_IT_2Mb.tsv", header = T)
nonRoH_Het_IT$O.HET <- nonRoH_Het_IT$NumberVariableSites-nonRoH_Het_IT$ObservedHomozygous
nonRoH_Het_IT$nonRoH_F.HET <- nonRoH_Het_IT$O.HET/(nonRoH_Het_IT$GenotypedSites-nonRoH_Het_IT$MissingSites)
nonRoH_Het_IT$TotalNumberSites <- nonRoH_Het_IT$GenotypedSites-nonRoH_Het_IT$MissingSites
nonRoH_Het_IT <- nonRoH_Het_IT[,-c(2:6,8)]
RoH_Het_IT <- merge(RoH_HET_IT, nonRoH_Het_IT, by = "Sample")
RoH_Het_IT$IDRisk <- RoH_Het_IT$FRoH*RoH_Het_IT$nonRoH_F.HET*1000

RoH_Het_IT %>%
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
100*3223382/114021718

p1 <- ggplot(RoH_HET_IT, aes(y = geno_F.HET, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0.001,0.0045) + theme_bw()
p2 <- ggplot(RoH_Het_IT, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  geom_text_repel(data = subset(RoH_Het_IT, Origin == "Int-ITA"), aes(label = Abbpop), color="black") +
  theme_bw()
p3 <- ggplot(RoH_HET_IT, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0,0.125) + theme_bw()
p1RISK <- ggplot(RoH_Het_IT, aes(y = nonRoH_F.HET, x = FRoH, color=IDRisk)) +
  geom_point(size=3) + 
  scale_color_gradient2(
    low      = "#FFCE03",   # very light salmon
    mid      = "orange",   # medium pink‑red
    high     = "#de2d26",   # dark red
    midpoint = median(RoH_Het_IT$IDRisk, na.rm = TRUE),   # or any value you choose
    space    = "Lab"       # smoother perceptual interpolation
  ) +
  theme_bw() +theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))

# Same for FR
RoH_FR <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_bcftools_FR_2Mb.txt", header = T)
RoH_FR <- merge(RoH_FR, Lizards, by.x = "Sample", by.y = "ID")

geno_Het_FR <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/heterozygosity_summary_genomewide_FR.tsv", header = T)
geno_Het_FR$O.HET <- geno_Het_FR$NumberVariableSites-geno_Het_FR$ObservedHomozygous
geno_Het_FR$geno_F.HET <- geno_Het_FR$O.HET/(geno_Het_FR$GenotypedSites-geno_Het_FR$MissingSites)
geno_Het_FR <- geno_Het_FR[,-(2:6)]
RoH_HET_FR <- merge(RoH_FR, geno_Het_FR, by = "Sample")

nonRoH_Het_FR <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/heterozygosity_summary_nonRoH_FR_2Mb.tsv", header = T)
nonRoH_Het_FR$O.HET <- nonRoH_Het_FR$NumberVariableSites-nonRoH_Het_FR$ObservedHomozygous
nonRoH_Het_FR$nonRoH_F.HET <- nonRoH_Het_FR$O.HET/(nonRoH_Het_FR$GenotypedSites-nonRoH_Het_FR$MissingSites)
nonRoH_Het_FR$TotalNumberSites <- nonRoH_Het_FR$GenotypedSites-nonRoH_Het_FR$MissingSites
nonRoH_Het_FR <- nonRoH_Het_FR[,-c(2:6,8)]
RoH_Het_FR <- merge(RoH_HET_FR, nonRoH_Het_FR, by = "Sample")
RoH_Het_FR$IDRisk <- RoH_Het_FR$FRoH*RoH_Het_FR$nonRoH_F.HET*1000

RoH_Het_FR %>%
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

p4 <- ggplot(RoH_HET_FR, aes(y = geno_F.HET, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0.001,0.0045) + theme_bw()
p5 <- ggplot(RoH_Het_FR, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  geom_text_repel(data = subset(RoH_Het_FR, Origin == "Int-FRA"), aes(label = Abbpop), color="black") +
  theme_bw()
p6 <- ggplot(RoH_HET_FR, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0,0.125) + theme_bw()
p2RISK <- ggplot(RoH_Het_FR, aes(y = nonRoH_F.HET, x = FRoH, color=IDRisk)) +
  geom_point(size=3) + 
  scale_color_gradient2(
    low      = "#FFCE03",   # very light salmon
    mid      = "orange",   # medium pink‑red
    high     = "#de2d26",   # dark red
    midpoint = median(RoH_Het_FR$IDRisk, na.rm = TRUE),   # or any value you choose
    space    = "Lab"       # smoother perceptual interpolation
  ) +
  theme_bw() +theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/HetRoH_V2.pdf", height=12, width=10, useDingbats = F)
(p1+p4)/(p2+p5)/(p3+p6)
dev.off()

### Statistics
## Heterozygosity
# Wilcoxon test
wilcox.test(geno_F.HET ~ Origin, data = RoH_HET_IT)
# Cliff's delta
cliff.delta(subset(RoH_HET_IT, Origin == "Nat-ITA")$geno_F.HET, subset(RoH_HET_IT, Origin == "Int-ITA")$geno_F.HET)
# Wilcoxon test
wilcox.test(geno_F.HET ~ Origin, data = RoH_HET_FR)
# Cliff's delta
cliff.delta(subset(RoH_HET_FR, Origin == "Nat-FRA")$geno_F.HET, subset(RoH_HET_FR, Origin == "Int-FRA")$geno_F.HET)

## ROH Length
# Length of ROHs comparison 
Ita_BCF_model <- glmer(Length ~ Origin + (1 | Sample), 
                       family = Gamma(link = "log"), 
                       data = RoH_HET_IT)
summary(Ita_BCF_model)

# Get the proportion of the effect 
exp(fixef(Ita_BCF_model)[-1]) #  ~ meaning 0.7X difference.

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
wilcox.test(FRoH ~ Origin, data = RoH_HET_IT)
wilcox.test(FRoH ~ Origin, data = RoH_HET_FR)

#####
# ROH for Supplements (Plink + different length)
#####

RoH_IT_Supp <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_bcftools_IT_500kb.txt", header = T)
RoH_IT_Supp <- merge(RoH_IT_Supp, Lizards, by.x = "Sample", by.y = "ID")
RoH_IT_Supp$Length <- RoH_IT_Supp$Length/1000000
p1S <- ggplot(RoH_IT_Supp, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  geom_text_repel(data = subset(RoH_IT_Supp, Origin == "Int-ITA"), aes(label = Abbpop), color="black") +
  theme_bw()
p2S <- ggplot(RoH_IT_Supp, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0,0.31) + theme_bw()

RoH_FR_Supp <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_bcftools_FR_500kb.txt", header = T)
RoH_FR_Supp <- merge(RoH_FR_Supp, Lizards, by.x = "Sample", by.y = "ID")
RoH_FR_Supp$Length <- RoH_FR_Supp$Length/1000000
p3S <- ggplot(RoH_FR_Supp, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  geom_text_repel(data = subset(RoH_FR_Supp, Origin == "Int-FRA"), aes(label = Abbpop), color="black") +
  theme_bw()
p4S <- ggplot(RoH_FR_Supp, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0,0.31) + theme_bw()

RoH_IT_plink_2Mb<- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_plink_IT_2Mb.hom.indiv", header = T)
RoH_IT_plink_2Mb$FRoH <- RoH_IT_plink_2Mb$KB/RoH_IT_plink_2Mb$KB[RoH_IT_plink_2Mb$FID == "Pxx"]
RoH_IT_plink_2Mb <- subset(RoH_IT_plink_2Mb, FID != "Pxx")
RoH_IT_plink_2Mb$Length <- RoH_IT_plink_2Mb$KB/1000
colnames(RoH_IT_plink_2Mb) <- c("Sample","IID","PHE","nRoH","KB","KBAVG","FRoH","Length")
RoH_IT_plink_2Mb <- merge(RoH_IT_plink_2Mb, Lizards, by.x = "Sample", by.y = "ID")
p_pl_IT_2Mb_1 <- ggplot(RoH_IT_plink_2Mb, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  geom_text_repel(data = subset(RoH_IT_plink_2Mb, Origin == "Int-ITA"), aes(label = Abbpop), color="black") +
  theme_bw()
p_pl_IT_2Mb_2 <- ggplot(RoH_IT_plink_2Mb, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0,0.18) + theme_bw()

RoH_FR_plink_2Mb <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_plink_FR_2Mb.hom.indiv", header = T)
RoH_FR_plink_2Mb$FRoH <- RoH_FR_plink_2Mb$KB/RoH_FR_plink_2Mb$KB[RoH_FR_plink_2Mb$FID == "Pxx"]
RoH_FR_plink_2Mb <- subset(RoH_FR_plink_2Mb, FID != "Pxx")
RoH_FR_plink_2Mb$Length <- RoH_FR_plink_2Mb$KB/1000
colnames(RoH_FR_plink_2Mb) <- c("Sample","IID","PHE","nRoH","KB","KBAVG","FRoH","Length")
RoH_FR_plink_2Mb <- merge(RoH_FR_plink_2Mb, Lizards, by.x = "Sample", by.y = "ID")
p_pl_FR_2Mb_1 <- ggplot(RoH_FR_plink_2Mb, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  geom_text_repel(data = subset(RoH_FR_plink_2Mb, Origin == "Int-FRA"), aes(label = Abbpop), color="black") +
  theme_bw()
p_pl_FR_2Mb_2 <- ggplot(RoH_FR_plink_2Mb, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0,0.18) + theme_bw()

RoH_IT_plink_500kb <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_plink_IT_500kb.hom.indiv", header = T)
RoH_IT_plink_500kb$FRoH <- RoH_IT_plink_500kb$KB/RoH_IT_plink_500kb$KB[RoH_IT_plink_500kb$FID == "Pxx"]
RoH_IT_plink_500kb <- subset(RoH_IT_plink_500kb, FID != "Pxx")
RoH_IT_plink_500kb$Length <- RoH_IT_plink_500kb$KB/1000
colnames(RoH_IT_plink_500kb) <- c("Sample","IID","PHE","nRoH","KB","KBAVG","FRoH","Length")
RoH_IT_plink_500kb <- merge(RoH_IT_plink_500kb, Lizards, by.x = "Sample", by.y = "ID")
p_pl_IT_500kb_1 <- ggplot(RoH_IT_plink_500kb, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  geom_text_repel(data = subset(RoH_IT_plink_500kb, Origin == "Int-ITA"), aes(label = Abbpop), color="black") +
  theme_bw()
p_pl_IT_500kb_2 <- ggplot(RoH_IT_plink_500kb, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#EFE808","#DD6A27")) + scale_fill_manual(values=c("#EFE808","#DD6A27")) +
  ylim(0,0.35) + theme_bw()

RoH_FR_plink_500kb <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/froh_summary_plink_FR_500kb.hom.indiv", header = T)
RoH_FR_plink_500kb$FRoH <- RoH_FR_plink_500kb$KB/RoH_FR_plink_500kb$KB[RoH_FR_plink_500kb$FID == "Pxx"]
RoH_FR_plink_500kb <- subset(RoH_FR_plink_500kb, FID != "Pxx")
RoH_FR_plink_500kb$Length <- RoH_FR_plink_500kb$KB/1000
colnames(RoH_FR_plink_500kb) <- c("Sample","IID","PHE","nRoH","KB","KBAVG","FRoH","Length")
RoH_FR_plink_500kb <- merge(RoH_FR_plink_500kb, Lizards, by.x = "Sample", by.y = "ID")
p_pl_FR_500kb_1 <- ggplot(RoH_FR_plink_500kb, aes(y = nRoH, x = Length, color=Origin)) +
  geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  geom_text_repel(data = subset(RoH_FR_plink_500kb, Origin == "Int-FRA"), aes(label = Abbpop), color="black") +
  theme_bw()
p_pl_FR_500kb_2 <- ggplot(RoH_FR_plink_500kb, aes(y = FRoH, x = Origin, color=Origin)) +
  geom_boxplot(aes(fill=Origin, alpha=0.8), show.legend = F) + geom_point(size=3, show.legend = F) + scale_color_manual(values=c("#6BCBDA","#0673B3")) + scale_fill_manual(values=c("#6BCBDA","#0673B3")) +
  ylim(0,0.35) + theme_bw()

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/HetRoH_Supp_V1.pdf", height=10, width=14, useDingbats = F)
(p1S | p2S | p3S | p4S) /
(p_pl_IT_2Mb_1 | p_pl_IT_2Mb_2 | p_pl_FR_2Mb_1 | p_pl_FR_2Mb_2) /
(p_pl_IT_500kb_1 | p_pl_IT_500kb_2 | p_pl_FR_500kb_1 | p_pl_FR_500kb_2)
dev.off()

RoH_IT_plink_2Mb %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_FR_plink_2Mb %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_IT_Supp %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_FR_Supp %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_IT_plink_500kb %>%
  group_by(Origin) %>%
  summarise(
    mean_nRoH  = mean(nRoH, na.rm = TRUE),
    min_nRoH   = min(nRoH, na.rm = TRUE),
    max_nRoH   = max(nRoH, na.rm = TRUE),
    total_nRoH = sum(nRoH, na.rm = TRUE),
    mean_Length = mean(Length, na.rm = TRUE),
    mean_FRoH = mean(FRoH, na.rm = TRUE),
    n          = n())
RoH_FR_plink_500kb %>%
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
Ita_BCF <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen//roh.pseudo.qual.IT.2Mb", h=F)

# Naming the columns
Ita_BCF <- Ita_BCF %>%
  dplyr::select(-1) %>% # Remove the first column - a bunch of RG
  rename(ID = V2, CHR = V3, POS1 = V4, POS2 = V5, BP = V6, MARKERS = V7, QUALITY = V8) 

# Merge the clean data with the origin
Ita_BCF<- merge(Lizards,Ita_BCF,by="ID", all.x=T)
Ita_BCF <- Ita_BCF %>%
  filter(!grepl("FRA", Origin))

# Define names  of chromosome names to  numbers
Ita_BCF <- Ita_BCF %>%
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

# French 
Fra_BCF <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/roh.pseudo.qual.FR.2Mb",h=F)

# Naming the columns and filtering for quality (Phred score) and minimum length (500k)
Fra_BCF <- Fra_BCF %>%
  dplyr::select(-1) %>% # Remove the first column - a bunch of RG
  rename(ID = V2, CHR = V3, POS1 = V4, POS2 = V5, BP = V6, MARKERS = V7, QUALITY = V8) 

# Merge the clean data with the origin
Fra_BCF<- merge(Lizards,Fra_BCF,by="ID", all.x=T)
Fra_BCF <- Fra_BCF %>%
  filter(!grepl("ITA", Origin))

# Change chromosome names 
Fra_BCF <- Fra_BCF %>%
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
Ita_BCF <- Ita_BCF %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))
Fra_BCF <- Fra_BCF %>%
  mutate(CHR = factor(CHR, levels = sort(as.numeric(as.character(unique(CHR))))))

# Filter for Chromosome 1
Ita_BCF_chr1 <- Ita_BCF %>% filter(CHR == 1)

# Get all unique IDs
all_IDs <- unique(Ita_BCF$ID)

# Create a lookup table for Origin and Abbpop per ID
ID_lookup <- Ita_BCF %>%
  dplyr::select(ID, Origin, Abbpop) %>%
  distinct()

# Fill missing IDs for Chr1, keeping Origin and Abbpop
Ita_BCF_chr1_filled <- Ita_BCF_chr1 %>%
  right_join(
    tibble(ID = all_IDs),
    by = "ID"
  ) %>%
  # add Origin/Abbpop from lookup if NA
  left_join(ID_lookup, by = "ID", suffix = c("", ".lookup")) %>%
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
Ita_BCF_chr1_filled <- Ita_BCF_chr1_filled %>%
  mutate(
    ID = factor(ID, levels = unique(ID[order(Origin)])))

# Plot only Chr 1
chr1_IT <- ggplot(Ita_BCF_chr1_filled, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), linewidth =3) +  
  scale_color_manual(values= c("Nat-ITA" = "#0673B3", "Int-ITA" = "#6BCBDA")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), limits = c(0, 130727322)) +
  theme_minimal() + 
  labs(x="Genomic Position", y="Sample", 
       color="Origin") +
  theme(strip.text = element_text(size=12), 
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))

### same for France

# Filter for Chromosome 1
Fra_BCF_chr1 <- Fra_BCF %>% filter(CHR == 1)

# Get all unique IDs
all_IDs <- unique(Fra_BCF$ID)

# Create a lookup table for Origin and Abbpop per ID
ID_lookup <- Fra_BCF %>%
  dplyr::select(ID, Origin, Abbpop) %>%
  distinct()

# Fill missing IDs for Chr1, keeping Origin and Abbpop
Fra_BCF_chr1_filled <- Fra_BCF_chr1 %>%
  right_join(
    tibble(ID = all_IDs),
    by = "ID"
  ) %>%
  # add Origin/Abbpop from lookup if NA
  left_join(ID_lookup, by = "ID", suffix = c("", ".lookup")) %>%
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
Fra_BCF_chr1_filled <- Fra_BCF_chr1_filled %>%
  mutate(
    ID = factor(ID, levels = unique(ID[order(Origin)])))

# Plot only Chr 1
chr1_FR <- ggplot(Fra_BCF_chr1_filled, aes(x=POS1, xend=POS2, y=ID, color=as.factor(Origin))) +
  geom_segment(aes(yend=ID), linewidth =3) +  
  scale_color_manual(values= c("Nat-FRA" = "#DD6A27", "Int-FRA" = "#EFE808")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1), limits = c(0, 130727322)) +
  theme_minimal() + 
  labs(x="Genomic Position", y="Sample", 
       color="Origin") +
  theme(strip.text = element_text(size=12), 
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))


pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Chr1_V1.pdf", height=6, width=10, useDingbats = F)
chr1_IT/chr1_FR
dev.off()

###########################
#### IDrisk
###########################

Other_IDRisk <- readWorkbook("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/PopGen/Kyriazis2025_TREE_RoHs_Supplements.xlsx", sheet = 1)
Other_IDRisk <- Other_IDRisk[c(2,6,7,11,24),c(1,5)]
colnames(Other_IDRisk) <- c("Sample","IDRisk")
Other_IDRisk$Origin <- "Lit"

RoH_Het_IT <- RoH_Het_IT %>%
  arrange(desc(IDRisk)) %>%  # sort by IDRisk descending
  mutate(Sample = factor(Sample, levels = Sample))

RoH_Het_FR <- RoH_Het_FR %>%
  arrange(desc(IDRisk)) %>%  # sort by IDRisk descending
  mutate(Sample = factor(Sample, levels = Sample))

common_cols <- intersect(colnames(RoH_Het_IT), colnames(Other_IDRisk))

df_combined <- rbind(
  RoH_Het_IT[common_cols],
  RoH_Het_FR[common_cols],
  Other_IDRisk[common_cols]
)

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/IDRisk_V2.pdf", height=6, width=6, useDingbats = F)
ggplot(df_combined, aes(x = Sample, y = IDRisk, fill = Origin)) +
  geom_col() +
  scale_fill_manual(values = c("Int-ITA" = "#6BCBDA", "Nat-ITA" = "#0673B3","Nat-FRA" = "#DD6A27", "Int-FRA" = "#EFE808")) +
  theme_bw() +
  labs(x = "ID Risk", y = "Value", title = "Barplot of ID Risk by Origin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))
dev.off()

###########################
# 3. PURGING  -----
###########################
# 3.1 Calculate Rxy for Italian-origin samples ----

# Load Italian Purging dataset - Relative frequencies and Jackknifing information
Ita_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Italian_RefSeq_freq.tsv", h=T)

# Calculate Lxy and Lyx ratios
Ita_Purging <- Ita_Purging %>%
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
Rxy_Italian <- Ita_Purging %>%
  dplyr::select(Rxy_High, Rxy_Moderate, Rxy_Low, Rxy_Modifier)

# Reshape data into long format for analysis
Rxy_Italian <- Rxy_Italian %>% 
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

calculate_minmax(Rxy_Italian)

# Visualize Italian-origin 
# Boxplot
Ita_RefSeq <- ggplot(Rxy_Italian, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#0673B3",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle =  "Italian",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.97, 1.04))


# 3.2 Calculate Rxy for French-origin samples ----
Fra_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_French_RefSeq_freq.tsv", h=T)

# Get Lxy and Lyx 
Fra_Purging <- Fra_Purging %>%
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
Rxy_French<- Fra_Purging %>%
  dplyr::select(Rxy_High,
                Rxy_Moderate,
                Rxy_Low,
                Rxy_Modifier)

# Transform data 
Rxy_French<- Rxy_French %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Visualize French-origin 
# BoxPlot 
Fra_RefSeq <- ggplot(Rxy_French, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#DD6A27",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle = "French",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.77, 1.04))

qq<- Ita_RefSeq + Fra_RefSeq
print(qq)  



# Confirm data as factors and check min max values after jackknifing
calculate_minmax(Rxy_French)



# 3.3 Calculate Rxy for each population ----

# Define populations
italian_pops <- c("BB", "DL", "NF", "SH", "VT", "WS")
french_pops <- c("BU", "WB", "WE")

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
    Population %in% italian_pops ~ "Italian",
    Population %in% french_pops ~ "French",
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
# Create Italian plot with green color
boxplot_Italian <- Rxy_all %>%
  filter(Group == "Italian") %>%
  ggplot(aes(x = Rxy, y = Impact)) +
  geom_boxplot(fill = "#0673B3",  # Green for Italian
               outlier.shape = NA, 
               alpha = 0.6) +
  facet_wrap(~ Population, nrow = 2) +
  theme_bw() +
  scale_x_continuous(limits = c(0.75, 1.1)) +
  labs(x = "Rxy", y = "Impact Category") +
  theme(legend.position = "none")

# Create French plot with orange color
boxplot_French <- Rxy_all %>%
  filter(Group == "French") %>%
  ggplot(aes(x = Rxy, y = Impact)) +
  geom_boxplot(fill = "#DD6A27",  # Orange for French
               outlier.shape = NA, 
               alpha = 0.6) +
  facet_wrap(~ Population, nrow = 1) +
  theme_bw() +
  scale_x_continuous(limits = c(0.75, 1.1)) +
  labs(x = "Rxy", y = "Impact Category") +
  theme(legend.position = "none")


qq<- (Ita_RefSeq + Fra_RefSeq) / boxplot_Italian / boxplot_French + plot_layout(heights = c(1, 1, 0.5))
pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Purging_V1.pdf", height=12, width=10, useDingbats = F)
print(qq)
dev.off()


#### Purging for Supplementary Material

# Load Italian Purging dataset - Relative frequencies and Jackknifing information
Ita_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Italian_Ensemble_freq.tsv", h=T)

# Calculate Lxy and Lyx ratios
Ita_Purging <- Ita_Purging %>%
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
Rxy_Italian <- Ita_Purging %>%
  dplyr::select(Rxy_High, Rxy_Moderate, Rxy_Low, Rxy_Modifier)

# Reshape data into long format for analysis
Rxy_Italian <- Rxy_Italian %>% 
  pivot_longer(cols = starts_with("Rxy_"),   # Select Rxy columns
               names_to = "Impact",          # Create column for categories
               values_to = "Rxy") %>%        # Assign values to Rxy
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean category names

calculate_minmax(Rxy_Italian)

# Visualize Italian-origin 
# Boxplot
Ita_Ensemble <- ggplot(Rxy_Italian, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#0673B3",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle =  "Italian",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.90, 1.04))

# Load Italian Purging dataset - Relative frequencies and Jackknifing information
Ita_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_Italian_Tiberius_freq.tsv", h=T)

# Calculate Lxy and Lyx ratios
Ita_Purging <- Ita_Purging %>%
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
Rxy_Italian <- Ita_Purging %>%
  dplyr::select(Rxy_High, Rxy_Moderate, Rxy_Low, Rxy_Modifier)

# Reshape data into long format for analysis
Rxy_Italian <- Rxy_Italian %>% 
  pivot_longer(cols = starts_with("Rxy_"),   # Select Rxy columns
               names_to = "Impact",          # Create column for categories
               values_to = "Rxy") %>%        # Assign values to Rxy
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean category names

calculate_minmax(Rxy_Italian)

# Visualize Italian-origin 
# Boxplot
Ita_Tiberius <- ggplot(Rxy_Italian, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#0673B3",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle =  "Italian",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.90, 1.04))

# 3.2 Calculate Rxy for French-origin samples ----
Fra_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_French_Ensemble_freq.tsv", h=T)

# Get Lxy and Lyx 
Fra_Purging <- Fra_Purging %>%
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
Rxy_French<- Fra_Purging %>%
  dplyr::select(Rxy_High,
                Rxy_Moderate,
                Rxy_Low,
                Rxy_Modifier)

# Transform data 
Rxy_French<- Rxy_French %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Visualize French-origin 
# BoxPlot 
Fra_Ensemble <- ggplot(Rxy_French, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#DD6A27",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle = "French",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.76, 1.04))

Fra_Purging <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Purging/All_French_Tiberius_freq.tsv", h=T)

# Get Lxy and Lyx 
Fra_Purging <- Fra_Purging %>%
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
Rxy_French<- Fra_Purging %>%
  dplyr::select(Rxy_High,
                Rxy_Moderate,
                Rxy_Low,
                Rxy_Modifier)

# Transform data 
Rxy_French<- Rxy_French %>% 
  pivot_longer(cols = starts_with("Rxy_"),  # Select all R_ columns
               names_to = "Impact",       # New column for impact categories
               values_to = "Rxy") %>%      # New column for Rxy values
  mutate(Impact = str_remove(Impact, "Rxy_")) # Clean up impact names

# Visualize French-origin 
# BoxPlot 
Fra_Tiberius <- ggplot(Rxy_French, aes(x = Rxy, y = Impact, fill = Impact)) +
  geom_boxplot(fill = "#DD6A27",  
               outlier.shape = NA, 
               alpha = 0.8) +
  labs(subtitle = "French",
       x = "Rxy",
       y = "Impact Category") +
  theme_bw() +
  scale_x_continuous(limits = c(0.76, 1.04))

pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/Purging_Supp.pdf", height=8, width=10, useDingbats = F)
(Ita_Ensemble + Fra_Ensemble) / (Ita_Tiberius + Fra_Tiberius)
dev.off()



#####################################################
### Plot the position of Recombination HotSpots
#####################################################

# ---- 1. Read chromosome sizes ----
chrom <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Chromosome_Length.txt", h=F)
colnames(chrom) <- c("chr", "start", "end")

# ---- 2. Read VCF (skip header lines starting with ##) ----
HotSpots <- read.table("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/RecombHotSpots/All_Merged_Hotspots.bed", header = T)

HotSpots <- HotSpots %>%
  mutate(chr = case_when(
    chr == "CM014743.1" ~ "NC_041312.1",
    chr == "CM014744.1" ~ "NC_041313.1",
    chr == "CM014745.1" ~ "NC_041314.1",
    chr == "CM014746.1" ~ "NC_041315.1",
    chr == "CM014747.1" ~ "NC_041316.1",
    chr == "CM014748.1" ~ "NC_041317.1",
    chr == "CM014749.1" ~ "NC_041318.1",
    chr == "CM014750.1" ~ "NC_041319.1",
    chr == "CM014751.1" ~ "NC_041320.1",
    chr == "CM014752.1" ~ "NC_041321.1",
    chr == "CM014753.1" ~ "NC_041322.1",
    chr == "CM014754.1" ~ "NC_041323.1",
    chr == "CM014755.1" ~ "NC_041324.1",
    chr == "CM014756.1" ~ "NC_041325.1",
    chr == "CM014757.1" ~ "NC_041326.1",
    chr == "CM014758.1" ~ "NC_041327.1",
    chr == "CM014759.1" ~ "NC_041328.1",
    chr == "CM014760.1" ~ "NC_041329.1",))

# ---- 3. Plot ----
pdf("C:/Users/feiner/Dropbox/MS_UK_wallies/Plots/RecombHotspots.pdf", height=6, width=16, useDingbats = F)
ggplot() + geom_segment(data = chrom, aes(x = start, xend = end, y = chr, yend = chr), size = 3, color = "grey70") + 
  geom_point(data = HotSpots, aes(x = start, y = chr, size = mean_rate),  alpha=0.1, color = "red") +
  theme_bw() + ggtitle("Hotspots (N=4113") + labs(x = "Genomic position", y = "Chromosome") + theme(panel.grid = element_blank())
dev.off()


##########################
### Position of HIGH impact variants
##########################


library(ggplot2)
library(dplyr)

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

### Now same with overlapping variants:

# ---- 2. Read VCF (skip header lines starting with ##) ----
vcf <- read.table("Y:/Projects/UKwallies/scripts/shared_sites_RefSeq_Tiberius_renamed.vcf.gz", comment.char = "#", header = FALSE)

# Extract relevant columns
colnames(vcf)[1:2] <- c("chr", "pos")
vcf <- vcf[, c("chr", "pos")]

# ---- 3. Plot ----
p_RefSeqTib <- ggplot() + geom_segment(data = chrom, aes(x = start, xend = end, y = chr, yend = chr), size = 3, color = "grey70") + 
  geom_point(data = vcf, aes(x = pos, y = chr), size = 2, alpha = 0.7, color = "red") +
  theme_bw() + ggtitle("RefSeq_Tiberius (N=127)") + labs(x = "Genomic position", y = "Chromosome") + theme(panel.grid = element_blank())


# ---- 2. Read VCF (skip header lines starting with ##) ----
vcf <- read.table("Y:/Projects/UKwallies/scripts/shared_sites_renamed.vcf.gz", comment.char = "#", header = FALSE)

# Extract relevant columns
colnames(vcf)[1:2] <- c("chr", "pos")
vcf <- vcf[, c("chr", "pos")]

# ---- 3. Plot ----
p_RefSeqEns <- ggplot() + geom_segment(data = chrom, aes(x = start, xend = end, y = chr, yend = chr), size = 3, color = "grey70") + 
  geom_point(data = vcf, aes(x = pos, y = chr), size = 2, alpha = 0.7, color = "red") +
  theme_bw() + ggtitle("RefSeq_Ensemble (N=100)") + labs(x = "Genomic position", y = "Chromosome") + theme(panel.grid = element_blank())

print(p_RefSeqTib+p_RefSeqEns)







########## Everything below is not needed!





## Genetic load (RefSeq High impact variants; N=382)


###
### Genetic load (i.e., High impact variants)
High_alleles <- read.table("Y:/Projects/UKwallies/Purging/Polarised_Impacts_NCBI/genotypes_high.tsv", h=T)
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
Lizards <- readWorkbook("C:/Users/feiner/Dropbox/MS_UK_wallies/Data/Samples_Santiago.xlsx", sheet = 2)
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

#############
#### NE estimation
#############

# ------------------------------
# 1. Read bcftools ROH file
# ------------------------------
roh <- read.table("Y:/Projects/UKwallies/Results_RoH_Het_FR/roh.pseudo.qual",
                  header = FALSE,
                  comment.char = "#",
                  stringsAsFactors = FALSE)

# Assign column names based on your file
colnames(roh) <- c("RG", "IND", "CHR", "START", "END",
                   "LENGTH_BP", "NSNP", "QUAL")

# Convert length to Mb
roh$MB <- roh$LENGTH_BP / 1e6

# Remove dummy sample
roh <- subset(roh, IND != "Pxx_Hz")
roh <- subset(roh, MB >= 0.1)


# ------------------------------
# 2. Define ROH bins (Mb)
# ------------------------------
bins <- c(0.1, 0.5, 1, 2, 4, 8, 16)

roh$bin <- cut(roh$MB, breaks = bins, include.lowest = TRUE, right = FALSE)

# ------------------------------
# 3. Summarise per individual + bin
# ------------------------------
roh_summary <- roh %>%
  group_by(IND, bin) %>%
  summarise(total_MB = sum(MB), .groups = "drop")

# ------------------------------
# 4. Genome length (adjust!)
# ------------------------------
genome_length <- 1400  # Mb (change for your species)

roh_summary$fROH <- roh_summary$total_MB / genome_length

# ------------------------------
# 5. Compute bin midpoints automatically
# ------------------------------
# Extract numeric bin limits
get_midpoint <- function(bin_label) {
  nums <- as.numeric(unlist(regmatches(bin_label, gregexpr("[0-9.]+", bin_label))))
  mean(nums)
}

roh_summary$mid_MB <- sapply(as.character(roh_summary$bin), get_midpoint)

# ------------------------------
# 6. Convert to generations
# ------------------------------
roh_summary$T_gen <- 100 / (2 * roh_summary$mid_MB)

# ------------------------------
# 7. Estimate Ne
# ------------------------------
roh_summary$Ne <- 1 / (2 * roh_summary$fROH)
roh_summary <- merge(roh_summary, LizardsData, by.x="IND", by.y="ID")

# ------------------------------
# 8. Plot (per individual)
# ------------------------------
ggplot(roh_summary, aes(x = T_gen, y = Ne, color = Group, group = IND)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values= c("#6BCBDA","#0673B3")) +
  scale_x_reverse() +
  labs(x = "Generations ago", y = "Estimated Ne") +
  theme_minimal()


