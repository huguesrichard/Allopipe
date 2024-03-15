library(tidyverse)
library(data.table)
library(glue)
library(ggforce)
library(plotly)
library(ggrepel)

setwd("../../output/netMHCpan_out")

aams_columns = read.csv("AMS_AAMS_df.csv")
aams = read.csv("AMS_AAMS_long_format.csv")

ggplot(aams,
       aes(x=type,
           y=value,
           label=pair)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_sina(alpha=0.5)+
  labs(title="AMS values vs AAMS values")
ggplotly(aams_violin_plot)

ams_comp <- ggplot(aams_columns,
                       aes(x=AMS,
                           y=AAMS,
                           label=pair))+
  geom_point() +
  geom_text_repel(position=position_jitter(width=1,height=1),max.overlaps=Inf)+
  geom_smooth()+
  labs(title=glue("Comparison of AMS and AAMS"))
ggplotly(ams_comp)


##### OUTDATED



fasta_r383 = read.table("../indiv_vcf/joint_genotyping/hard-filtered/R383/test_kmers.fa",
                        comment.char = ">",
                        header=FALSE)
r383_d_netmhc = read.csv("R383-D0-N1.csv")
r383_r_netmhc = read.csv("R383-R0-N1.csv")
fasta_r383$V1 %in% r383_d_netmhc$Peptide %>% sum()
fasta_r383$V1 %in% r383_r_netmhc$Peptide %>% sum()

fasta_r908 = read.table("../indiv_vcf/joint_genotyping/hard-filtered/R908/test_kmers.fa",
                        comment.char = ">",
                        header=FALSE)
r908_d_netmhc = read.csv("R908-D0-N1.csv")
r908_r_netmhc = read.csv("R908-R0-N1.csv")
fasta_r908$V1 %in% r908_d_netmhc$Peptide %>% sum()
fasta_r908$V1 %in% r908_r_netmhc$Peptide %>% sum()

fasta_r1048 = read.table("../indiv_vcf/joint_genotyping/hard-filtered/R1048/test_kmers.fa",
                         comment.char = ">",
                         header=FALSE)
r1048_d_netmhc = read.csv("R1048-D0-N1.csv")
r1048_r_netmhc = read.csv("R1048-R0-N1.csv")
fasta_r1048$V1 %in% r1048_d_netmhc$Peptide %>% sum()
fasta_r1048$V1 %in% r1048_r_netmhc$Peptide %>% sum()

fasta_r6764 = read.table("../indiv_vcf/joint_genotyping/hard-filtered/R6764/test_kmers.fa",
                         comment.char = ">",
                         header=FALSE)
r6764_d_netmhc = read.csv("R6764-D0-N2.csv")
r6764_r_netmhc = read.csv("R6764-R0-N1.csv")
fasta_r6764$V1 %in% r6764_d_netmhc$Peptide %>% sum()
fasta_r6764$V1 %in% r6764_r_netmhc$Peptide %>% sum()


ggplot(r383_d_netmhc,
       aes(x=EL_Rank))+
  geom_histogram(bins = 100)+
  labs(title=glue("Histogram of EL_Rank counts for R383 donor"))

ggplot(r383_r_netmhc,
       aes(x=EL_Rank))+
  geom_histogram(bins = 100)+
  labs(title=glue("Histogram of EL_Rank counts for R383 recipient"))

ggplot(r908_d_netmhc,
       aes(x=EL_Rank))+
  geom_histogram(bins = 100)+
  labs(title=glue("Histogram of EL_Rank counts for R908 donor"))

ggplot(r908_r_netmhc,
       aes(x=EL_Rank))+
  geom_histogram(bins = 100)+
  labs(title=glue("Histogram of EL_Rank counts for R908 recipient"))

ggplot(r1048_d_netmhc,
       aes(x=EL_Rank))+
  geom_histogram(bins = 100)+
  labs(title=glue("Histogram of EL_Rank counts for R1048 donor"))

ggplot(r1048_r_netmhc,
       aes(x=EL_Rank))+
  geom_histogram(bins = 100)+
  labs(title=glue("Histogram of EL_Rank counts for R1048 recipient"))

ggplot(r6764_d_netmhc,
       aes(x=EL_Rank))+
  geom_histogram(bins = 100)+
  labs(title=glue("Histogram of EL_Rank counts for R6764 donor"))

ggplot(r6764_r_netmhc,
       aes(x=EL_Rank))+
  geom_histogram(bins = 100)+
  labs(title=glue("Histogram of EL_Rank counts for R6764 recipient"))


merged_pep = read.csv("../indiv_vcf/joint_genotyping/hard-filtered/merged_peptide_dataframe.csv")
merged_test = merged_pep %>% group_by(peptide) %>% 
  mutate(size = n())
merged_test$indiv <- factor(merged_test$indiv)
nb_common =  merged_pep %>% group_by(peptide) %>% 
  summarize(npairs=n(), pairs= indiv %>% unique() %>% 
              sort() %>% str_c(collapse=","))

ggplot(filter(merged_test,size>=3), aes(y=indiv, x=peptide)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 90))+
  labs(title="Heatmap of peptides shared by pairs",
       x="peptides",
       y="pairs")
