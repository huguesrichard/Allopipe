library(tidyverse)
library(glue)
library(ggforce)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(data.table)


setwd("../../output/indiv_vcf/joint_genotyping/hard-filtered/R5482")

theme_set(theme_bw())
graph_theme = theme(plot.title = element_text(face="bold",size = 12),
                    panel.border = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "grey50"),
                    axis.ticks = element_line(colour = "grey60", size = 0.2),
                    panel.grid.major = element_line(colour = "grey60", size = 0.2),
                    legend.justification="top")

dataset = read.csv("DP_unfiltered_donor_recipient_5482_pair.csv")

dp_recipient <- ggplot(dataset,
       aes(x=X1_y))+
  geom_histogram(binwidth = 10)+
  labs(title="Sequencing Depth counts in recipient of pair R5482 on all positions",
       x="Sequencing Depth",
       y="Counts")
dp_recipient
#ggsave("../../../../plots/plots_gdoc/methods/dp_recipient_R5482.png",
#       width = 751,
#       height = 468,
#       unit = "px",
#       dpi = 100)

dp_donor <- ggplot(dataset,
                       aes(x=X1_x))+
  geom_histogram(binwidth = 10)+
  labs(title="Sequencing Depth counts in donor of pair R5482 on all positions",
       x="Sequencing Depth",
       y="Counts")
dp_donor


mismatches = dataset %>% filter(X0_x != X0_y)
dp_donor_mismatches <- ggplot(mismatches,
                        aes(x=X1_x))+
  geom_histogram(binwidth = 10)+
  labs(title="Sequencing Depth counts in donor of pair R5482 on mismatched positions",
     x="Sequencing Depth",
     y="Counts")
dp_donor_mismatches

dp_recipient_mismatches <- ggplot(mismatches,
                              aes(x=X1_y))+
  geom_histogram(binwidth = 10)+
  labs(title="Sequencing Depth counts in recipient of pair R5482 on mismatched positions",
       x="Sequencing Depth",
       y="Counts")
dp_recipient_mismatches

dp_mismatches <- ggplot(mismatches,
                        aes(x=X1_x,
                            y=X1_y))+
  geom_point()+
  labs(title="Sequencing Depth counts of pair R5482 on mismatched positions",
       x="Sequencing Depth Donor",
       y="Sequencing Depth Recipient")
dp_mismatches

dp_mismatches_400 <- ggplot(mismatches %>% filter(X1_x<=400 & X1_y <=400),
                        aes(x=X1_x,
                            y=X1_y))+
  geom_point(shape="+",size=3)+
  labs(title="Sequencing Depth counts of pair R5482 on mismatched positions with DP<= 400",
       x="Sequencing Depth Donor",
       y="Sequencing Depth Recipient")
dp_mismatches_400


