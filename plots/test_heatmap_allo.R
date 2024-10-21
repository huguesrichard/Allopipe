library(tidyverse)
library(glue)

setwd("../../output/important_files")
a = read_csv("counts_heatmap_10000_20_400_5_0.8_0.csv")


b=pivot_longer(a, 2:last_col(), 
               names_to = "pos", values_to = "AMS_count") %>% 
  filter(!is.na(AMS_count))


b_filt = b %>% group_by(pos) %>% 
  mutate(size = n(), 
         size.pos = glue('{sprintf("%02d", size)}pairs-chr{pos}')) %>% 
  filter(size >= 11)
  
b_filt$pair <- factor(b_filt$pair)
nb_common =  b %>% group_by(pos) %>% 
  summarize(npairs=n(), pairs= pair %>% unique() %>% 
              sort() %>% str_c(collapse=","),
            AMS_mean = mean(AMS_count) %>% signif(d=2),
            AMS_max = max(AMS_count),
            AMS_min = min(AMS_count)) %>% 
  arrange(desc(npairs))

ggplot(b_filt, aes(y= pair, x=size.pos, fill = AMS_count) ) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 90))+
  labs(title="Heatmap of AMS regions shared by 11 or more pairs",x="number of pairs per position")
