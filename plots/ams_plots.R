library(tidyverse)
library(glue)
library(ggforce)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(data.table)
setwd("../../output/important_files")

theme_set(theme_bw())
boxplot_theme = theme(plot.title = element_text(face="bold",size = 12),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "grey50"),
                      axis.ticks = element_line(colour = "grey60", size = 0.2),
                      panel.grid.major = element_line(colour = "grey60", size = 0.2))


dataset_clin = read.csv("csh_clinical_ams.csv")
dataset_clin = read.csv("csh_clinical_ams_with_snps_and_dp_with_gnomad.csv")
dataset_clin$GVHc=factor(dataset_clin$GVHc,levels=c("yes","no","unknown"))
dataset_clin$GVHa=factor(dataset_clin$GVHa,levels=c("yes","no","unknown"))
dataset_clin$GVH=factor(dataset_clin$GVH,levels=c("yes","no","unknown"))
#acute and all
gvha_all <- ggplot(dataset_clin,
                   aes(x=GVHa,
                       y=AMS_stop,
                       fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AMS vs GVHa with HSC data",
       x="Acute GVH",
       y="AMS",
       fill="Acute GVH")
gvha_all
#chronic and all
gvhc_all <- ggplot(dataset_clin,
                   aes(x=GVHc,
                       y=AMS_stop,
                       fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AMS vs GVHc with HSC data",
       x="Chronic GVH",
       y="AMS",
       fill="Chronic GVH")
gvhc_all


# Acute GVH vs no GVH
gvha <- ggplot(dataset_clin %>% filter(GVH_status %in% c("acute_only","no_gvh")),
               aes(x=GVH_status,
                   y=AMS_gnomad,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AMS in patients segregated on acute GVH response for HSC data",
       x="GVH response to the graft",
       y="AMS",
       fill="GVH response")+
  boxplot_theme
gvha
ggsave("../plots/plots_gdoc/gvh_clinical/gvha_gnomad_ams_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)


# Chronic GVH vs no GVH
gvhc <- ggplot(dataset_clin %>% filter(GVH_status %in% c("chronic_only","no_gvh")),
               aes(x=GVH_status,
                   y=AMS_gnomad,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AMS in patients segregated on chronic GVH response for HSC data",
       x="GVH response to the graft",
       y="AMS",
       fill="GVH response")+
  boxplot_theme
gvhc
ggsave("../plots/plots_gdoc/gvh_clinical/gvhc_gnomad_ams_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)


dataset_clin$GVH_status=factor(dataset_clin$GVH_status,levels=c("acute_only","no_gvh","chronic_only"))
# chronic and acute GVH vs no GVH
gvh <- ggplot(dataset_clin %>% filter(GVH_status %in% c("acute_only","no_gvh","chronic_only")),
              aes(x=GVH_status,
                  y=AMS_gnomad,
                  fill=GVH_status,
                  label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AMS in patients segregated on GVH response for HSC data",
       x="GVH response to the graft",
       y="AMS",
       fill="GVH response")+
  boxplot_theme
gvh
ggsave("../plots/plots_gdoc/gvh_clinical/gvh_gnomad_ams_800_500.png",
       width = 800,
       height = 500,
       unit = "px",
       dpi = 100)


