library(tidyverse)
library(glue)
library(ggforce)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(data.table)
library(boot)

setwd("../../output/important_files")

AMS_vip = read.csv("AMS_vip_df.csv")
test <- ggplot(AMS_vip,
       aes(x=AMS,
           y=AMS_vip,
           label=pair))+
  geom_point()+
  geom_text_repel(position=position_jitter(width=1,height=1),max.overlaps=Inf)+
  geom_smooth()+
  labs(title=glue("Comparison of AMS and filtered on VIP regions"))+
  xlab("Classic Allogenomics Mismatch Score")+
  ylab("AMS on VIP regions")
test + theme_bw()

theme_set(theme_bw())

graph_theme = theme(plot.title = element_text(face="bold",size = 12),
                    panel.border = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "grey50"),
                    axis.ticks = element_line(colour = "grey60", size = 0.2),
                    panel.grid.major = element_line(colour = "grey60", size = 0.2),
                    legend.justification="top")
test + graph_theme

data_depth = read.csv("unfiltered_depth.csv")

depth <- ggplot(data_depth,
       aes(x=DP))+
  geom_histogram(bins=100)+
  labs(title=glue("Histogram of Read Depth values counts in a VCF of one individual"),
       x = "Read Depth",
       y = "occurrences in the unfiltered VCF")

depth + graph_theme
ggsave("../plots/Depth_R383_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)
full_df = read.csv("csh_clinical_ams_with_snps_and_dp_with_gnomad.csv")

snps_fx <- ggplot(full_df,
                  aes(x=DP_donor,
                      y=SNPs_donor,
                      label=pair))+
  geom_point(aes(color = AMS_stop))+
  geom_text_repel(data=subset(full_df,SNPs_donor>1.1*mean(SNPs_donor)),position=position_jitter(width=1,height=1),max.overlaps=Inf,aes(label=pair))+
  geom_smooth()+
  labs(title=glue("Depth and SNPs count for donors"))+
  scale_color_continuous(name="AMS")
snps_fx+graph_theme

ggsave("../plots/plots_gdoc/ams_control/AMS_DP_SNPs_renamed_legend_751_468.png",
       width = 751,
       height = 468,
       unit = "px",
       dpi = 100)

ams <- subset(full_df,select = c("pair","AMS_stop","AMS_norm"))
ams <- ams %>% rename("AMS"="AMS_stop",
               "norm_AMS"="AMS_norm")
ams <- ams %>% pivot_longer(cols = c("AMS","norm_AMS"),
                            names_to = "type_AMS",
                            values_to = "AMS_value")
normalized <- ggplot(ams,
                     aes(x=type_AMS,
                         y=AMS_value))+
  geom_point()+
  geom_line(aes(group=pair))+
  labs(title("Effect of the normalization of the AMS on the order of the pairs"))
normalized


AMS_dataset = read.csv("AMS_df_fisher.csv")
fisher <- ggplot(AMS_dataset,
       aes(x=AMS_pre_fisher,
           y=AMS_post_fisher,
           label=pair)) +
  geom_point() +
  geom_text_repel(position=position_jitter(width=1,height=1),max.overlaps=Inf)+
  geom_smooth()+
  labs(title=glue("Comparison of AMS before and after Fisher filtering"))

fisher
ggplotly(fisher)

dataset_clin = read.csv("csh_clinical_ams.csv")
ams_intervals = read.csv("csh_2021_12_13_counts_chr_10000_20_400_5_0.8_0.csv")
clin_intervals = dataset_clin %>% inner_join(ams_intervals,by="pair")

clin_filt = clin_intervals %>% group_by(interval) %>% 
  mutate(size=n(),
         size.interval = glue('{sprintf("%02d", size)}pairs-chr{interval}'))

clin_filt$pair <- factor(clin_filt$pair)
common_intervals =   clin_intervals %>% group_by(interval) %>% 
  summarize(npairs=n(), pairs= pair %>% unique() %>% 
              sort() %>% str_c(collapse=","),
            AMS_mean = mean(AMS_giab) %>% signif(d=2),
            AMS_max = max(AMS_giab),
            AMS_min = min(AMS_giab),
            nbGVHa = uniqueN(pair[GVH_status=="acute_only"]),
            nbGVHc = uniqueN(pair[GVH_status=="chronic_only"]),
            nbGVH = uniqueN(pair[GVH_status == "both_gvh"]),
            nbnoGVH = uniqueN(pair[GVH_status == "no_gvh"]),
            pGVHa = mean(GVH_status=="acute_only"),
            pGVHc = mean(GVH_status=="chronic_only"),
            pGVH = mean(GVH_status=="both_gvh"),
            pnoGVH = mean(GVH_status=="no_gvh")) %>% 
  arrange(desc(npairs))
gvh_cases <- dataset_clin %>%
  group_by(GVH_status) %>%
  count()
gvh_cases$GVH_status[gvh_cases$GVH_status=="acute_only"] <- "nbGVHa"
gvh_cases$GVH_status[gvh_cases$GVH_status=="both_gvh"] <- "nbGVH"
gvh_cases$GVH_status[gvh_cases$GVH_status=="chronic_only"] <- "nbGVHc"
gvh_cases$GVH_status[gvh_cases$GVH_status=="no_gvh"] <- "nbnoGVH"
gvh_cases <- gvh_cases %>%
  rename(columns=GVH_status,
         counts = n)
sum_indiv <- sum(gvh_cases$counts)
gvh_cases <- gvh_cases %>% pivot_wider(names_from = columns,values_from = counts)
gvh_cases$total = sum_indiv

mymat = matrix(nrow = 17931,ncol=8)
for(i in 1:17931){
  tab1 = matrix(c(gvh_cases$total-gvh_cases$nbGVHa,gvh_cases$nbGVHa,common_intervals$npairs[i]-common_intervals$nbGVHa[i],common_intervals$nbGVHa[i]),nrow=2)
  a <- fisher.test(tab1)
  tab2 = matrix(c(gvh_cases$total-gvh_cases$nbGVHc,gvh_cases$nbGVHc,common_intervals$npairs[i]-common_intervals$nbGVHc[i],common_intervals$nbGVHc[i]),nrow=2)
  b <- fisher.test(tab2)
  tab3 = matrix(c(gvh_cases$total-gvh_cases$nbGVH,gvh_cases$nbGVH,common_intervals$npairs[i]-common_intervals$nbGVH[i],common_intervals$nbGVH[i]),nrow=2)
  c <- fisher.test(tab3)
  tab4 = matrix(c(gvh_cases$total-gvh_cases$nbnoGVH,gvh_cases$nbnoGVH,common_intervals$npairs[i]-common_intervals$nbnoGVH[i],common_intervals$nbnoGVH[i]),nrow=2)
  d <- fisher.test(tab4)
  mymat[i,1] = a$estimate
  mymat[i,2] = b$estimate
  mymat[i,3] = c$estimate
  mymat[i,4] = d$estimate
  mymat[i,5] = a$p.value
  mymat[i,6] = b$p.value
  mymat[i,7] = c$p.value
  mymat[i,8] = d$p.value
}
mydf <- as.data.frame(mymat)
common_intervals$odds_GVHa = mydf$V1
common_intervals$odds_GVHc = mydf$V2
common_intervals$odds_GVH = mydf$V3
common_intervals$odds_noGVH = mydf$V4
common_intervals$pval_GVHa = mydf$V5
common_intervals$pval_GVHc = mydf$V6
common_intervals$pval_GVH = mydf$V7
common_intervals$pval_noGVH = mydf$V8
common_intervals = common_intervals %>% filter((npairs>=15|(pval_GVHa<0.02|pval_GVHc<0.02|pval_GVH<0.02|pval_noGVH<0.02)))

common_intervals$adj_pval_GVHa = p.adjust(common_intervals$pval_GVHa,"fdr")
common_intervals$adj_pval_GVHc = p.adjust(common_intervals$pval_GVHc,"fdr")
common_intervals$adj_pval_GVH = p.adjust(common_intervals$pval_GVH,"fdr")
common_intervals$adj_pval_noGVH = p.adjust(common_intervals$pval_noGVH,"fdr")
write.csv(common_intervals,"table_pour_adele.csv",row.names = FALSE)

# BH/fdr ends up not being the right method to adjust p-values
# too many 1 values
plot(sort(common_intervals$pval_GVHa),log="y")
plot(sort(common_intervals$pval_GVHc),log="y")
plot(sort(common_intervals$pval_GVH),log="y")
plot(sort(common_intervals$pval_noGVH),log="y")


# replace counts_heatmap.csv with desired counts file in important files repo
counts_heatmap = read.csv("counts_heatmap.csv")
counts_heatmap[is.na(counts_heatmap)] = 0
counts_num = as.matrix(counts_heatmap[,2:length(counts_heatmap)])
rownames(counts_num) = sapply(counts_heatmap$pair,function(x)
  as.character(x))
pheatmap(counts_num,cluster_cols=FALSE ,main="Heatmap of AMS on windows of genome in pair, window = 100k, thresh = 5")
ggplot(counts_heatmap, aes(y= pair, x=size.pos, fill = AMS_count) ) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 90))+
  labs(title="Heatmap of AMS regions shared by 11 or more pairs",x="number of pairs per position")

# OUTDATED

dataset_kidney = read.csv("kidney.csv")
dataset_kidney_r = read.csv("kidney_relations.csv")
dataset_kidney_rel_rm = read.csv("kidney_relations_without_9_22_24.csv")
dataset_kidney_rm_pairs = subset(dataset_kidney,!(pair %in% c(9,22,24)))


snps_DP_marrow = read.csv("SNPs_DP_df_filt.csv")
snps_d_marrow = read.csv("SNPs_donors_reverse_ecdf.csv")
snps_r_marrow = read.csv("SNPs_recipients_reverse_ecdf.csv")
snps_d_marrow_giab = read.csv("SNPs_donors_reverse_ecdf_giab.csv")
snps_r_marrow_giab = read.csv("SNPs_recipients_reverse_ecdf_giab.csv")
marrow_AMS_comp = read.csv("marrow_AMS.csv")
norm_test = read.csv("normtest.csv")
ams_marrow = read.csv("marrow_ams_SNPs_DP.csv")
dataset_marrow = read.csv("marrow.csv")
dataset_clin = read.csv("csh_clin_giab_norm_dropna_plots_test.csv")
dataset_clin = read.csv("csh_clinical_ams.csv")

ggplotWithLabels <- function(dataset,col1,col2,point_colors){
  require(ggplot2)
  ggplot(dataset,
         aes_string(x=col1,
             y=col2,
             label="pair",
             col=point_colors)) +
    geom_point() +
    geom_text_repel(position=position_jitter(width=1,height=1),max.overlaps=Inf)+
    geom_smooth()+
    labs(title=glue("Comparison of {col1} and {col2}"))
}


ggplotWithLabels2 <- function(dataset,col1,col2,point_colors,xlim,ylim){
  require(ggplot2)
  ggplot(dataset,
         aes_string(x=col1,
                    y=col2,
                    col=point_colors)) +
    geom_line(size=0.3) +
    guides(col="none")+
    scale_x_continuous(col1,limits=xlim)+
    scale_y_continuous(col2,limits=ylim)+
    labs(title=glue("Comparison of {col1} and {col2}"))
}

ggplotWithLabels2(snps_d_marrow,"DP","r_ecdf","pair",c(10,50),c(65000,77000))
ggplotWithLabels2(snps_r_marrow,"DP","r_ecdf","pair",c(10,50),c(65000,77000))
ggplotWithLabels2(snps_d_marrow_giab,"DP","r_ecdf","pair",c(10,50),c(50000,60000))
ggplotWithLabels2(snps_r_marrow_giab,"DP","r_ecdf","colpair",c(10,50),c(50000,60000))


ggplot(dataset_clin,
       aes(x=GVHc %>% as.character(),
           y=AMS_semi_outer,
           fill=GVHc %>% as.character())) +
  geom_boxplot(notch = TRUE,notchwidth = 0.5) +
  geom_jitter() + 
  labs(title="Boxplot of raw AMS vs chronical GVH with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")

ggplot(dataset_clin,
       aes(x=GVHc %>% as.character(),
           y=AMS_norm)) +
  geom_boxplot(notch = TRUE,notchwidth = 0.5) +
  geom_jitter() + 
  labs(title="Boxplot of normalized AMS vs chronical GVH with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")

ggplot(dataset_clin,
       aes(x=GVHa %>% as.character(),
           y=AMS_semi_outer)) +
  geom_boxplot(notch = TRUE,notchwidth = 0.5) +
  geom_jitter() + 
  labs(title="Boxplot of raw AMS vs acute GVH with HSC data",
       x="Acute GVH",
       fill="Acute GVH")

ggplot(dataset_clin,
       aes(x=GVHa %>% as.character(),
           y=AMS_norm)) +
  geom_boxplot(notch = TRUE,notchwidth = 0.5) +
  geom_jitter() + 
  labs(title="Boxplot of normalized AMS vs acute GVH with HSC data")

ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " pval =",signif(summary(fit)$coef[2,4], 5)))+
    xlab("Ratio of Ref positions")+
    ylab("Classic AMS")
}
fit_lm <- lm(AMS_stop~ref_ratio,dataset_clin)
plotthis <- ggplotRegression(fit_lm)
plotthis + theme(plot.title = element_text(face="bold"),
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line = element_line(colour = "grey50"),
             aspect.ratio = 0.5)


dataset_clin$GVHc=factor(dataset_clin$GVHc,levels=c("yes","no","unknown"))
dataset_clin$GVHa=factor(dataset_clin$GVHa,levels=c("yes","no","unknown"))
dataset_clin$GVH=factor(dataset_clin$GVH,levels=c("yes","no","unknown"))

l <- ggplot(dataset_clin,
       aes(x=GVHc %>% as.character(),
           y=AMS_giab,
           fill=GVHc %>% as.character())) +
  geom_violin() +
  geom_sina(alpha=0.5) + 
  labs(title="Violin/Sinaplot of raw AMS vs GVHc with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")

l+ e

# without coloring by sex of the patients
# RAW AMS VS GHV_A
ggplot(dataset_clin,
       aes(x=GVHa,
           y=AMS_giab_gb,
           fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of raw AMS vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")

ggplot(dataset_clin,
       aes(x=GVHa,
           y=test,
           fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of normalized AMS vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")



# RAW AMS VS GVH_C
ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVHc,
           y=AMS_giab_gb,
           fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5) + 
  labs(title="Violin/Sinaplot of raw AMS vs GVHc with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")

# NORM AMS VS GVH_C
ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVHc,
           y=test,
           fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5) + 
  labs(title="Violin/Sinaplot of normalized AMS vs GVHc with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")

# NORM AMS VS GVH
ggplot(dataset_clin,
       aes(x=GVH,
           y=test,
           fill=GVH)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of normalized AMS vs GVH with HSC data",
       x="GVH",
       fill="GVH")

# coloring by sex of patients
# RAW AMS VS GVH_A
ggplot(dataset_clin,
       aes(x=GVHa,
           y=AMS_giab_gb,
           fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(aes(col=pair_sex,group=GVHa)) + 
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Violin/Sinaplot of raw AMS vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")
# NORM AMS VS GVH_A
ggplot(dataset_clin,
       aes(x=GVHa,
           y=AMS_giab_norm,
           fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(aes(col=pair_sex,group=GVHa)) + 
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Violin/Sinaplot of normalized AMS vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")
# RAW AMS VS GVH_C
ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVHc,
           y=AMS_giab_gb,
           fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_brewer(palette="Pastel1")+
  geom_sina(aes(col=pair_sex,group=GVHc)) + 
  labs(title="Violin/Sinaplot of raw AMS vs GVHc with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")


# NORM AMS VS GVH_C
ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVHc,
           y=AMS_giab_norm,
           fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_brewer(palette="Pastel1")+
  geom_sina(aes(col=pair_sex,group=GVHc)) + 
  labs(title="Violin/Sinaplot of normalized AMS vs GVHc with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")

ggplot(dataset_clin,
       aes(x=GVH,
           y=AMS_giab_norm,
           fill=GVH)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(aes(col=pair_sex,group=GVH)) + 
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Violin/Sinaplot of normalized AMS vs GVH with HSC data",
       x="GVH",
       fill="GVH")

#### 
ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVHc,
           y=AMS_giab_norm,
           col=pair_sex,
           fill=GVHc)) +
  geom_violin() +
  scale_fill_brewer(palette="Pastel1")+
  geom_sina(alpha=0.5) + 
  labs(title="Violin/Sinaplot of normalized AMS vs GVHc with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")

# RAW AMS SEX PAIR
ggplot(dataset_clin,
       aes(x=pair_sex,
           y=AMS_giab_gb,
           fill=pair_sex)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5) + 
  labs(title="Violin/Sinaplot of raw AMS vs sex mismatch with HSC data",
       x="Sex of the patients (D-R)",
       fill="Sex of the patients (D-R)")
# NORM AMS SEX PAIR
ggplot(dataset_clin,
       aes(x=pair_sex,
           y=AMS_giab_norm,
           fill=pair_sex)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5) + 
  labs(title="Violin/Sinaplot of normalized AMS vs sex mismatch with HSC data",
       x="Sex of the patients (D-R)",
       fill="Sex of the patients (D-R)")

# GVHa SEX MISMATCH
ggplot(dataset_clin,
       aes(x=GVHa,
           y=sex_mismatch,
           fill=GVHa)) +
  geom_jitter() +
  labs(title="Violin/Sinaplot of GVHa vs sex mismatch with HSC data",
       x="GVHa",
       fill="GVHa")

# GVHc SEX MISMATCH
ggplot(dataset_clin,
       aes(x=GVHc,
           y=sex_mismatch,
           fill=GVHc)) +
  geom_violin() +
  geom_sina(alpha=0.5) + 
  labs(title="Violin/Sinaplot of GVHc vs sex mismatch with HSC data",
       x="GVHc",
       fill="GVHc")


# RAW AMS VS GVH_A
ggplot(dataset_clin,
       aes(x=GVHa,
           y=AMS_giab_gb,
           fill=pair_sex)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Violin/Sinaplot of raw AMS vs GVHa with HSC data",
       x="Acute GVH",
       fill="pair_sex")
# NORM AMS VS GVH_A
ggplot(dataset_clin,
       aes(x=GVHa,
           y=AMS_giab_norm,
           fill=pair_sex)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Violin/Sinaplot of normalized AMS vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")


# RAW AMS VS GVH_C
ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVHc,
           y=AMS_giab_gb,
           fill=pair_sex)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Violin/Sinaplot of raw AMS vs GVHc with HSC data",
       x="Chronical GVH",
       fill="Sex pair")
# NORM AMS VS GVH_C
ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVHc,
           y=AMS_giab_norm,
           fill=pair_sex)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Violin/Sinaplot of normalized AMS vs GVHc with HSC data",
       x="Chronical GVH",
       fill="Sex pair")

ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVH,
           y=AMS_giab_norm,
           fill=pair_sex)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_brewer(palette="Pastel1") +
  labs(title="Violin/Sinaplot of normalized AMS vs GVH with HSC data",
       x="GVH",
       fill="Sex pair")


ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVH,
           y=AMS_giab_norm,
           col=pair_sex,
           fill=GVH)) +
  geom_violin(draw_quantiles=0.5) +
  scale_fill_brewer(palette="Pastel1")+
  geom_sina(alpha=0.5) + 
  labs(title="Violin/Sinaplot of normalized AMS vs GVH with HSC data",
       x="GVH",
       fill="GVH")

#acute and all
gvha_all <- ggplot(dataset_clin,
                   aes(x=GVHa,
                       y=AMS_giab_norm,
                       fill=GVHa)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of normalized AMS vs GVHa with HSC data",
       x="Acute GVH",
       fill="Acute GVH")

#chronic and all
gvhc_all <- ggplot(dataset_clin,
                   aes(x=GVHc,
                       y=AMS_giab_norm,
                       fill=GVHc)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of normalized AMS vs GVHc with HSC data",
       x="Chronic GVH",
       fill="Chronic GVH")
# Acute GVH vs no GVH
gvha <- ggplot(subset(dataset_clin,GVH_status=="acute"|GVH_status=="none"),
       aes(x=GVH_status,
           y=test,
           fill=GVH_status,
           label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of normalized AMS vs acute GVH status with HSC data",
       x="Acute GVH",
       fill="Acute GVH")
ggplotly(gvha)

# Chronic GVH vs no GVH
ggplot(dataset_clin %>% filter(GVH_status %in% c("chronic_only","no_gvh")),
       aes(x=GVH_status,
           y=AMS_stop,
           fill=GVH_status,
           label=pair))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of AMS vs chronical GVH status with HSC data",
       x="Chronical GVH",
       fill="Chronical GVH")


ggplot(subset(dataset_clin,pair!="R2655"),
       aes(x=GVH,
           y=test,
           fill=GVH))+
  geom_violin(draw_quantiles = 0.5)+
  geom_sina(alpha=0.5)+
  labs(title="Violin/Sinaplot of normalized AMS vs GVH status with HSC data",
       x="GVH",
       fill="GVH")
ggarrange(gvha, gvhc, gvha_all, gvhc_all, gvh + rremove("x.text"), 
          labels = c("A", "B", "C","D","E"),
          ncol = 5, nrow = 1)

#This is the custom Student's test part

ggplot(dataset_clin,
       aes(x=AMS_vip,
           y=AMS_stop))+
  geom_point()+
  geom_hline(yintercept = 0.05)


a = dataset_clin  %>% pivot_longer(c(AMS_norm, AMS_stop:AAMS_vip) , 
                                   names_to = "AMS_type", values_to = "value")

extreme_pairs = c("R603", "R6241", "R3215", "R3682")

ggplot( data = a %>% filter(! (pair  %in% extreme_pairs), 
                            GVH_status %in% c("acute_only", "no_gvh") ), 
        aes(x = GVH_status, y = value, fill = GVH_status)) + 
  geom_violin( draw_quantiles = c(0.5)) + geom_point() + 
  facet_wrap(~AMS_type, scales = "free")


ggplot( data = a %>% filter(! (pair  %in% extreme_pairs), 
                            GVH_status %in% c("chronic_only", "no_gvh") ), 
        aes(x = GVH_status, y = value, fill = GVH_status)) + 
  geom_violin( draw_quantiles = c(0.5)) + geom_point() + 
  facet_wrap(~AMS_type, scales = "free")


samp_mean <- function(x,index){
  mean(x[index])
}


nogvh <- dataset_clin %>% filter(GVH_status=="no_gvh")
acute <- dataset_clin %>% filter(GVH_status=="acute_only")
chronic <- dataset_clin %>% filter(GVH_status=="chronic_only")
set.seed(123)
boot_acute_AMS <- boot(data=acute$AMS_stop,statistic=samp_mean,R=1000)
boot_chronic_AMS <- boot(data=chronic$AMS_stop,statistic=samp_mean,R=1000)
boot_nogvh_AMS <- boot(data=nogvh$AMS_stop,statistic=samp_mean,R=1000)

boot_acute_AMS
boot_chronic_AMS
boot_nogvh_AMS


boot_acute_AAMS_v1 <- boot(data=acute$AAMS_EL2,statistic=samp_mean,R=1000)
boot_chronic_AAMS_v1 <- boot(data=chronic$AAMS_EL2,statistic=samp_mean,R=1000)
boot_nogvh_AAMS_v1 <- boot(data=nogvh$AAMS_EL2,statistic=samp_mean,R=1000)

boot_acute_AAMS_v1
boot_chronic_AAMS_v1
boot_nogvh_AAMS_v1

boot_acute_AAMS_v2 <- boot(data=acute$AAMS_EL05,statistic=samp_mean,R=1000)
boot_chronic_AAMS_v2 <- boot(data=chronic$AAMS_EL05,statistic=samp_mean,R=1000)
boot_nogvh_AAMS_v2 <- boot(data=nogvh$AAMS_EL05,statistic=samp_mean,R=1000)


boot_nogvh_AAMS_v1$t %>% hist(br = 50)
boot_acute_AAMS_v1$t %>% hist(br = 50)
bind_rows(boot_acute_AAMS_v1$t, boot_nogvh_AAMS_v1$t, .id = c('acute', 'nogvh'))
bind_rows(list(boot_acute_AAMS_v1$t, boot_nogvh_AAMS_v1$t), .id = c('acute', 'nogvh'))
bind_rows(list(acute=boot_acute_AAMS_v1$t, nogvh=boot_nogvh_AAMS_v1$t), .id = c('acute', 'nogvh'))
bind_rows(list(acute=tibble(val = boot_acute_AAMS_v1$t), nogvh=tibble(val = boot_nogvh_AAMS_v1$t)), .id = 'type')
bind_rows(list(acute=tibble(val = boot_acute_AAMS_v1$t[,1]), nogvh=tibble(val = boot_nogvh_AAMS_v1$t[,1])), .id = 'type')
toto = bind_rows(list(acute=tibble(val = boot_acute_AAMS_v1$t[,1]), nogvh=tibble(val = boot_nogvh_AAMS_v1$t[,1])), .id = 'type')
toto
ggplot(data = toto, aes(x = val, fill = type)) + geom_histogram(position = 'dodge')
ggplot(data = toto, aes(x = val, fill = type)) + geom_density()
ggplot(data = toto, aes(x = val, fill = type)) + geom_density( alpha = 0.5)
toto = bind_rows(list(acute=tibble(val = boot_acute_AAMS_v1$t[,1]), nogvh=tibble(val = boot_nogvh_AAMS_v1$t[,1]), chronic = tibble(val = boot_chronic_AAMS_v1$t[,1])), .id = 'type')
ggplot(data = toto, aes(x = val, fill = type)) + geom_density( alpha = 0.5)
boot_acute_AAMS_v2
boot_chronic_AAMS_v2
boot_nogvh_AAMS_v2

boot_acute_AMS_vip <- boot(data=acute$AMS_vip,statistic=samp_mean,R=1000)
boot_chronic_AMS_vip <- boot(data=chronic$AMS_vip,statistic=samp_mean,R=1000)
boot_nogvh_AMS_vip <- boot(data=nogvh$AMS_vip,statistic=samp_mean,R=1000)

boot_acute_AMS_vip
boot_chronic_AMS_vip
boot_nogvh_AMS_vip

boot_acute_AAMS_vip <- boot(data=acute$AAMS_vip,statistic=samp_mean,R=1000)
boot_chronic_AAMS_vip <- boot(data=chronic$AAMS_vip,statistic=samp_mean,R=1000)
boot_nogvh_AAMS_vip <- boot(data=nogvh$AAMS_vip,statistic=samp_mean,R=1000)

boot_acute_AAMS_vip
boot_chronic_AAMS_vip
boot_nogvh_AAMS_vip



setwd("~/Documents/ALLOGENOMICS/AMS_workflow/output/important_files/")
netmhcpan_df = read.csv("NetMHCpan_full_df_non_filtered.csv")

ggplot(netmhcpan_df %>% filter(indiv %in% c("R383-D0","R383-R0")&EL_Rank<10),
       aes(x = EL_Rank,
           y = EL.score))+
  geom_point(size=0.001)

ggplot(netmhcpan_df,
       aes(x=EL_Rank))+
  geom_histogram(bins=100)


merged_mismatches_1000 = read.csv("geno_csh_merged_df_0_1000_0_08_20_3.csv")

ggplot(merged_mismatches_1000,
       aes(x=DP_x,
           y=DP_y,
           color=pair))+
  geom_point()



merged_mismatches = read.csv("geno_csh_merged_df_20_400_5_08_20_3.csv")

ggplot(merged_mismatches,
       aes(x=DP_x,
           y=DP_y,
           color=pair))+
  geom_point()



setwd("~/Documents/ALLOGENOMICS/AMS_workflow/output/indiv_vcf/Full_Match")

data_fm = read.csv("clin_full_match_pairs.csv")
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " pval =",signif(summary(fit)$coef[2,4], 5)))+
    xlab("Creatinine")+
    ylab("AMS Full Match")
}
theme_set(theme_bw())
graph_theme = theme(plot.title = element_text(face="bold",size = 12),
                    panel.border = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "grey50"),
                    axis.ticks = element_line(colour = "grey60", size = 0.2),
                    panel.grid.major = element_line(colour = "grey60", size = 0.2),
                    legend.justification="top")

fit_lm <- lm(AMS~creat_12+creat_24+Creat_36+creat_48,data_fm)
plotthis <- ggplotRegression(fit_lm)
plotthis + graph_theme
