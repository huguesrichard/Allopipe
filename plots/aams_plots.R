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
                       y=AAMS_EL2,
                       fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AAMS vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")

#chronic and all
gvhc_all <- ggplot(dataset_clin,
                   aes(x=GVHc,
                       y=AAMS_EL2,
                       fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AAMS vs GVHc with HSC data",
       x="Chronic GVH",
       fill="Chronic GVH")
# Acute GVH vs no GVH AAMS_EL2
gvha <- ggplot(dataset_clin %>% filter(GVH_status %in% c("acute_only","no_gvh")),
               aes(x=GVH_status,
                   y=AAMS_EL2,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AAMS v1 in patients segregated on acute GVH response for HSC data",
       x="GVH response to the graft",
       y="AAMS for EL ranks under 2%",
       fill="GVH response")+
  boxplot_theme
gvha
ggsave("../plots/plots_gdoc/gvh_clinical/gvha_el2_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)


# Chronic GVH vs no GVH AAMS_EL2
gvhc <- ggplot(dataset_clin %>% filter(GVH_status %in% c("chronic_only","no_gvh")),
               aes(x=GVH_status,
                   y=AAMS_EL2,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AAMS v1 in patients segregated on chronic GVH response for HSC data",
       x="GVH response to the graft",
       y="AAMS for EL ranks under 2%",
       fill="GVH response")+
  boxplot_theme
gvhc
ggsave("../plots/plots_gdoc/gvh_clinical/gvhc_el2_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)

# Acute GVH vs no GVH AAMS_EL05
gvha <- ggplot(dataset_clin %>% filter(GVH_status %in% c("acute_only","no_gvh")),
               aes(x=GVH_status,
                   y=AAMS_EL05,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AAMS v2 in patients segregated on acute GVH response for HSC data",
       x="GVH response to the graft",
       y="AAMS for EL ranks under 0.5%",
       fill="GVH response")+
  boxplot_theme
gvha
ggsave("../plots/plots_gdoc/gvh_clinical/gvha_el05_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)


# Chronic GVH vs no GVH AAMS_EL05
gvhc <- ggplot(dataset_clin %>% filter(GVH_status %in% c("chronic_only","no_gvh")),
               aes(x=GVH_status,
                   y=AAMS_EL05,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AAMS v2 in patients segregated on chronic GVH response for HSC data",
       x="GVH response to the graft",
       y="AAMS for EL ranks under 0.5%",
       fill="GVH response")+
  boxplot_theme
gvhc
ggsave("../plots/plots_gdoc/gvh_clinical/gvhc_el05_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)


gvh <- ggplot(dataset_clin,
       aes(x=GVH,
           y=AAMS_EL2,
           fill=GVH))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AAMS vs GVH status with HSC data",
       x="GVH",
       fill="GVH")
gvh
# ggarrange(gvha, gvhc, gvha_all, gvhc_all, gvh + rremove("x.text"), 
#          labels = c("A", "B", "C","D","E"),
#          ncol = 5, nrow = 1)

acute <- dataset_clin %>% filter(GVH_status %in% c("acute_only"))
chronic <- dataset_clin %>% filter(GVH_status %in% c("chronic_only"))
nogvh <- dataset_clin %>% filter(GVH_status %in% c("no_gvh"))
t.test(acute$AAMS_EL2,nogvh$AAMS_EL2)
t.test(chronic$AAMS_EL2,nogvh$AAMS_EL2)

ac <- dataset_clin %>% filter(GVH_status %in% c("acute_only","no_gvh"))
mood.test(ac$AMS_stop~ac$GVH_status)
mood.test(ac$AAMS_EL2~ac$GVH_status)
mood.test(ac$AAMS_EL05~ac$GVH_status)
ch <- dataset_clin %>% filter(GVH_status %in% c("chronic_only","no_gvh"))
mood.test(ch$AMS_stop~ch$GVH_status)
mood.test(ch$AAMS_EL2~ch$GVH_status)
mood.test(ch$AAMS_EL05~ch$GVH_status)


ggplot(dataset_clin,
       aes(x=AAMS_EL2,
           y=AAMS_EL05,
           label=pair))+
  geom_point()+
  geom_text_repel(position=position_jitter(width=1,height=1),max.overlaps=Inf)+
  geom_smooth()+
  labs(title="Comparison of AAMS with EL=2 and EL=0.5",
       x="AAMS EL=2",
       y="AAMS EL=0.5")


#acute and all
ggplot(dataset_clin,
       aes(x=GVH,
           y=AAMS_EL2,
           fill=GVH)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AAMS vs GVH with HSC data",
       x="GVH",
       fill="GVH")



# AAMS for GVH
dataset_clin$GVH_status=factor(dataset_clin$GVH_status,levels=c("acute_only","no_gvh","chronic_only","both_gvh"))
gvh <- ggplot(dataset_clin %>% filter(GVH_status %in% c("acute_only","no_gvh","chronic_only")),
               aes(x=GVH_status,
                   y=AAMS_EL05,
                   fill=GVH_status,
                   label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="AAMS v2 in patients segregated on GVH response for HSC data",
       x="GVH response to the graft",
       y="AAMS for EL ranks under 0.5%",
       fill="GVH response")+
  boxplot_theme
gvh
ggsave("../plots/plots_gdoc/gvh_clinical/gvh_aams_v2_800_500.png",
       width = 800,
       height = 500,
       unit = "px",
       dpi = 100)


