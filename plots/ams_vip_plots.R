library(tidyverse)
library(glue)
library(ggforce)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(data.table)
setwd("~/Documents/AMS_workflow/output/important_files")


dataset_clin = read.csv("csh_clinical_ams_with_snps_and_dp_with_gnomad.csv")

dataset_clin$GVHc=factor(dataset_clin$GVHc,levels=c("yes","no","unknown"))
dataset_clin$GVHa=factor(dataset_clin$GVHa,levels=c("yes","no","unknown"))
dataset_clin$GVH=factor(dataset_clin$GVH,levels=c("yes","no","unknown"))

theme_set(theme_bw())
boxplot_theme = theme(plot.title = element_text(face="bold",size = 12),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "grey50"),
                      axis.ticks = element_line(colour = "grey60", size = 0.2),
                      panel.grid.major = element_line(colour = "grey60", size = 0.2))



#acute and all
gvha_all <- ggplot(dataset_clin,
                   aes(x=GVHa,
                       y=AMS_vip,
                       fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AMS in VIP regions vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")

#chronic and all
gvhc_all <- ggplot(dataset_clin,
                   aes(x=GVHc,
                       y=AMS_vip,
                       fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AMS in VIP regions vs GVHc with HSC data",
       x="Chronic GVH",
       fill="Chronic GVH")
# Acute GVH vs no GVH
gvha <- ggplot(dataset_clin %>% filter(GVH_status %in% c("acute_only","no_gvh")),
               aes(x=GVH_status,
                   y=AMS_vip,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AMS in VIP regions in patients segregated on acute GVH response for HSC data",
       x="GVH response to the graft",
       y="AMS in VIP regions",
       fill="GVH response")+
  boxplot_theme
gvha
ggsave("../../../plots/plots_gdoc/gvh_clinical/gvha_ams_vip_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)


# Chronic GVH vs no GVH
gvhc <- ggplot(dataset_clin %>% filter(GVH_status %in% c("chronic_only","no_gvh")),
               aes(x=GVH_status,
                   y=AMS_vip,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AMS in VIP regions in patients segregated on chronic GVH response for HSC data",
       x="GVH response to the graft",
       y="AMS in VIP regions",
       fill="GVH response")+
  boxplot_theme
gvhc
ggsave("../../../plots/plots_gdoc/gvh_clinical/gvhc_ams_vip_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)


gvh <- ggplot(dataset_clin,
              aes(x=GVH,
                  y=AMS_vip,
                  fill=GVH))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AMS in VIP regions vs GVH status with HSC data",
       x="GVH",
       fill="GVH")



#acute and all
gvha_all <- ggplot(dataset_clin,
                   aes(x=GVHa,
                       y=AAMS_vip,
                       fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AAMS in VIP regions vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")

#chronic and all
gvhc_all <- ggplot(dataset_clin,
                   aes(x=GVHc,
                       y=AAMS_vip,
                       fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AAMS in VIP regions vs GVHc with HSC data",
       x="Chronic GVH",
       fill="Chronic GVH")
# Acute GVH vs no GVH
gvha <- ggplot(dataset_clin %>% filter(GVH_status %in% c("acute_only","no_gvh")),
               aes(x=GVH_status,
                   y=AAMS_vip,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AAMS in VIP regions in patients segregated on acute GVH response for HSC data",
       x="GVH response to the graft",
       y="AAMS in VIP regions",
       fill="GVH response")+
  boxplot_theme
gvha
ggsave("../../../plots/plots_gdoc/gvh_clinical/gvha_aams_vip_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)



# Chronic GVH vs no GVH
gvhc <- ggplot(dataset_clin %>% filter(GVH_status %in% c("chronic_only","no_gvh")),
               aes(x=GVH_status,
                   y=AAMS_vip,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AAMS in VIP regions in patients segregated on chronic GVH response for HSC data",
       x="GVH response to the graft",
       y="AAMS in VIP regions",
       fill="GVH response")+
  boxplot_theme
gvhc
ggsave("../../../plots/plots_gdoc/gvh_clinical/gvhc_aams_vip_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)


gvh <- ggplot(dataset_clin,
              aes(x=GVH,
                  y=AAMS_vip,
                  fill=GVH))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AAMS in VIP regions vs GVH status with HSC data",
       x="GVH",
       fill="GVH")

