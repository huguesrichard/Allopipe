library(tidyverse)
library(glue)
library(ggforce)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(plotly)
setwd("../../output/important_files")

theme_set(theme_bw())
graph_theme = theme(plot.title = element_text(face="bold",size = 12),
                    panel.border = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "grey50"),
                    axis.ticks = element_line(colour = "grey60", size = 0.2),
                    panel.grid.major = element_line(colour = "grey60", size = 0.2),
                    legend.justification="top")

# required number of pairs
# AMS
dataset <- read.csv("students_custom_pvalues.csv")
ggplot(dataset %>% filter(mismatch=="AMS"&condition=="acute_GVH"&test=="one-tail"&N<=100&forces==1),
       aes(x = N,
           y = log10(pval)))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red",size=0.7)+ 
  labs(title = "Estimation of required number of pairs for significant effect in acute GVH with AMS",
       x="N = alpha(pN1+pN2)")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  geom_point(aes(0,-2),color="red",size=0.1)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AMS_GVHa_pthres_0.01_width_900_height_500_diff.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)

ggplot(dataset %>% filter(mismatch=="AMS"&condition=="chronic GVH"&test=="one-tail"&N<=300&forces==1),
       aes(x = N,
           y = log10(pval)))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red",size=0.7)+ 
  labs(title = "Estimation of required number of pairs for significant effect in chronic GVH with AMS",
       x="N = alpha(pN1+pN2)")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  geom_point(aes(0,-2),color="red",size=0.1)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AMS_GVHc_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)

# AAMS

ggplot(dataset %>% filter(mismatch=="AAMS_EL2"&condition=="acute_GVH"&test=="one-tail"&N<=100&forces==1),
       aes(x = N,
           y = log10(pval)))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red",size=0.7)+ 
  labs(title = "Estimation of required number of pairs for significant effect in acute GVH with AAMS v1",
       x="N = alpha(pN1+pN2)")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  geom_point(aes(0,-2),color="red",size=0.1)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AAMS_v1_GVHa_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)

ggplot(dataset %>% filter(mismatch=="AAMS_EL2"&condition=="chronic GVH"&test=="one-tail"&N<=300&forces==1),
       aes(x = N,
           y = log10(pval)))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red",size=0.7)+ 
  labs(title = "Estimation of required number of pairs for significant effect in chronic GVH with AAMS v1",
       x="N = alpha(pN1+pN2)")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  geom_point(aes(0,-2),color="red",size=0.1)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AAMS_v1_GVHc_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)


ggplot(dataset %>% filter(mismatch=="AAMS_EL05"&condition=="acute_GVH"&test=="one-tail"&N<=150&forces==1),
       aes(x = N,
           y = log10(pval)))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red",size=0.7)+ 
  labs(title = "Estimation of required number of pairs for significant effect in acute GVH with AAMS v2",
       x="N = alpha(pN1+pN2)")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  geom_point(aes(0,-2),color="red",size=0.1)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AAMS_v2_GVHa_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)


ggplot(dataset %>% filter(mismatch=="AAMS_EL05"&condition=="chronic GVH"&test=="one-tail"&forces==1),
       aes(x = N,
           y = log10(pval)))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red",size=0.7)+ 
  labs(title = "Estimation of required number of pairs for significant effect in chronic GVH with AAMS v2",
       x="N = alpha(pN1+pN2)")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  geom_point(aes(0,-2),color="red",size=0.1)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AAMS_v2_GVHc_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)



# strength of effect

dataset <- read.csv("students_custom_pvalues.csv")
ggplot(dataset %>% filter(mismatch=="AMS"&condition=="acute_GVH"&test=="one-tail"&N<=100),
       aes(x = N,
           y = log10(pval),
           color = factor(forces)))+
  geom_line()+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red")+ 
  labs(title = "Effect of an increase of the strength of the effect on acute GVH data with AMS",
       x="N = alpha(pN1+pN2)",
       color="Strength factor")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AMS_strength_GVHa_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)

ggplot(dataset %>% filter(mismatch=="AMS"&condition=="chronic GVH"&test=="one-tail"&N<=200),
       aes(x = N,
           y = log10(pval),
           color = factor(forces)))+
  geom_line()+
  geom_hline(aes(yintercept=log10(0.01),linetype="pvalue = 0.01"),color="red")+ 
  labs(title = "Effect of an increase of the strength of the effect on chronic GVH data with AMS",
       x="N = alpha(pN1+pN2)",
       color = "Strength factor")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AMS_strength_GVHc_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)


ggplot(dataset %>% filter(mismatch=="AAMS_EL2"&condition=="acute_GVH"&test=="one-tail"&N<=200),
       aes(x = N,
           y = log10(pval),
           color = factor(forces)))+
  geom_line()+
  geom_hline(aes(yintercept=log10(0.01),linetype="pvalue = 0.01"),color="red")+ 
  labs(title = "Effect of an increase of the strength of the effect on acute GVH data with AAMS v1",
       x="N = alpha(pN1+pN2)",
       color = "Strength factor")+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  ylim(-20,0)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AAMS_strength_v1_GVHa_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)

ggplot(dataset %>% filter(mismatch=="AAMS_EL2"&condition=="chronic GVH"&test=="one-tail"&N<=200),
       aes(x = N,
           y = log10(pval),
           color = factor(forces)))+
  geom_line()+
  geom_hline(aes(yintercept=log10(0.01),linetype="pvalue = 0.01"),color="red")+ 
  labs(title = "Effect of an increase of the strength of the effect on chronic GVH data with AAMS v1",
       x="N = alpha(pN1+pN2)",
       color = "Strength factor")+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  ylim(-20,0)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AAMS_strength_v1_GVHc_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)

ggplot(dataset %>% filter(mismatch=="AAMS_EL05"&condition=="acute_GVH"&test=="one-tail"&N<=200),
       aes(x = N,
           y = log10(pval),
           color = factor(forces)))+
  geom_line()+
  geom_hline(aes(yintercept=log10(0.01),linetype="pvalue = 0.01"),color="red")+ 
  labs(title = "Effect of an increase of the strength of the effect on acute GVH data with AAMS v2",
       x="N = alpha(pN1+pN2)",
       color = "Strength factor")+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  ylim(-20,0)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AAMS_strength_v2_GVHa_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)

ggplot(dataset %>% filter(mismatch=="AAMS_EL05"&condition=="chronic GVH"&test=="one-tail"),
       aes(x = N,
           y = log10(pval),
           color = factor(forces)))+
  geom_line()+
  geom_hline(aes(yintercept=log10(0.01),linetype="pvalue = 0.01"),color="red")+ 
  labs(title = "Effect of an increase of the strength of the effect on chronic GVH data with AAMS v2",
       x="N = alpha(pN1+pN2)",
       color = "Strength factor")+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  ylim(-20,0)+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AAMS_strength_v2_GVHc_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)


# AMS and AAMS
ggplot(dataset %>% filter(mismatch %in% c("AMS","AAMS_EL2")&condition=="acute_GVH"&test=="one-tail"&N<=100&forces==1),
       aes(x = N,
           y = log10(pval),
           color=mismatch))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red",size=0.7)+ 
  labs(title = "Estimation of required number of pairs for significant effect in acute GVH with both scores",
       x="N = alpha(pN1+pN2)")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  geom_point(aes(0,-2),color="red",size=0.1)+
  scale_color_discrete(name="Mismatch type",labels=c("AAMS v1","AMS"))+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AMS_AAMS_v1_GVHa_pthres_0.01_width_900_height_500_diff.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)

ggplot(dataset %>% filter(mismatch %in% c("AMS","AAMS_EL2")&condition=="chronic GVH"&test=="one-tail"&N<=300&forces==1),
       aes(x = N,
           y = log10(pval),
           color=mismatch))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=log10(0.01), linetype = "pvalue = 0.01"),color="red",size=0.7)+ 
  labs(title = "Estimation of required number of pairs for significant effect in chronic GVH with both scores",
       x="N = alpha(pN1+pN2)")+
  scale_linetype_manual(name = "p-value threshold", values = c(2, 2), 
                        guide = guide_legend(override.aes = list(color = "red")))+
  scale_x_continuous(expand=c(0,0),breaks = scales::pretty_breaks(n = 10))+
  geom_point(aes(0,-2),color="red",size=0.1)+
  scale_color_discrete(name="Mismatch type",labels=c("AAMS v1","AMS"))+
  graph_theme
ggsave("../plots/plots_gdoc/custom_student/AMS_AAMS_v1_GVHc_pthres_0.01_width_900_height_500.png",
       width = 900,
       height = 500,
       unit = "px",
       dpi = 100)
