library(curatedMetagenomicData)
library(phyloseq)
library(mia)
library(vegan)
library(dplyr)
library(microeco)
library(file2meco)
library(ggplot2)
library(Maaslin2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

setwd("e:/ma16/")

load(file = "country_sample.csv")

country_crt = country_md_cr_t %>%
  returnSamples("relative_abundance")

country_crt1 = country_md_cr_t %>%
  returnSamples("relative_abundance", counts = T)

con_cr = makePhyloseqFromTreeSummarizedExperiment(country_crt, abund_values = "relative_abundance")

#======================================================================
df_otu = as.data.frame(otu_table(con_cr))
rownames(df_otu) = sapply(rownames(df_otu),function(x) unlist(strsplit(x, split = "\\|"))[7])
dfmd = country_md_cr_t
dfmd = dfmd %>% mutate(age_fl = case_when(age <= 45 ~ "n1",
                                   age > 45 ~ "n2"))
dfmd = as.data.frame(dfmd)
rownames(dfmd) = dfmd$sample_id

fit = Maaslin2(input_data = df_otu, input_metadata = dfmd, output = "ma_output",
               normalization = "NONE", transform = "LOG",
               fixed_effects = c("country", "gender","age_fl"),
               random_effects = "study_name",
               reference = "country,CHN",
               min_prevalence = 0.1,
               min_abundance = 0.0001)
fr = fit$results

mr_c_chn = fr %>% filter(metadata == "country" & qval < 0.1)
mr_c_usa = fr %>% filter(metadata == "country" & qval < 0.1)
mr_c_ind = fr %>% filter(metadata == "country" & qval < 0.1)
mr_c_mdg = fr %>% filter(metadata == "country" & qval < 0.1)
mr_c_gbr = fr %>% filter(metadata == "country" & qval < 0.1)
mr_c_nld = fr %>% filter(metadata == "country" & qval < 0.1)


mr_c_chn$c = "CHN"
mr_c_usa$c = "USA"
mr_c_ind$c = "IND"
mr_c_mdg$c = "MDG"
mr_c_gbr$c = "GBR"
mr_c_nld$c = "NLD"

dfd = mr_c_u
taxdf = as.data.frame(tax_table(con_cr))
rownames(taxdf) = sapply(rownames(taxdf),function(x) unlist(strsplit(x, split = "\\|"))[7])
dfd = cbind(dfd, taxdf[dfd[,1], c(1,2,6,7)])
sort(table(dfd$Phylum))

mr_c_u = rbind(mr_c_chn, mr_c_usa, mr_c_ind, mr_c_mdg, mr_c_gbr, mr_c_nld)
sort(table(mr_c_u[ ,1]))
df = mr_c_u[mr_c_u[,1] %in% names(b[b>=15]),]

dadf = table(mr_c_gbr[,1])
dadf[dadf>=5]
df = mr_c_gbr[mr_c_gbr[,1] %in% names(dadf[dadf>=5]),]

dadf = table(mr_c_nld[,1])
dadf[dadf>=5]

corec = table(mr_c_chn[,1])
corec = corec[corec>=3]
df_corec = data.frame()


mr_c_u = Reduce(union, list(mr_c_chn[,1], mr_c_usa[,1], mr_c_ind[,1], mr_c_mdg[,1], mr_c_gbr[,1], mr_c_nld[,1]))
mr_c_u = sapply(mr_c_u, function(x) unlist(strsplit(x, split = "__"))[2])
mr_c_u = sapply(mr_c_u, function(x) unlist(strsplit(x, split = "_"))[1])
sort(table(mr_c_u))

mr_c = mr_c_chn
corec = table(mr_c[,1])
corec = corec[corec>=3]
df_corec = data.frame()

for (i in 1:length(corec)){
  dfc = mr_c[mr_c[ ,1] == names(corec)[i], c(1,3,4)]
  df_corec = rbind(df_corec, dfc)
}

df_tax_w = tidyr::spread(df_corec, key = "feature", value = "coef", fill = 0)
rownames(df_tax_w) = df_tax_w[ ,1]
df_tax_w = t(df_tax_w[ ,-1])
pheatmap::pheatmap(df_tax_w, 
                   display_numbers = matrix(ifelse(df_tax_w != 0, signif(df_tax_w,3), ""), nrow(df_tax_w)),
                   fontsize = 16,
                   #fontfamily=  "serif" 
                   )


usa_r = nrow(mr_c %>% filter(value == "USA" & coef > 0))
usa_d = nrow(mr_c %>% filter(value == "USA" & coef < 0))
df_phy1 = data.frame(sample = "CHN VS USA", y = c(usa_r,usa_d), group = c("r","d"))

ind_r = nrow(mr_c %>% filter(value == "IND" & coef > 0))
ind_d = nrow(mr_c %>% filter(value == "IND" & coef < 0))
df_phy2 = data.frame(sample = "CHN VS IND", y = c(ind_r,ind_d), group = c("r","d"))

mdg_r = nrow(mr_c %>% filter(value == "MDG" & coef > 0))
mdg_d = nrow(mr_c %>% filter(value == "MDG" & coef < 0))
df_phy3 = data.frame(sample = "CHN VS MDG", y = c(mdg_r,mdg_d), group = c("r","d"))

gbr_r = nrow(mr_c %>% filter(value == "GBR" & coef > 0))
gbr_d = nrow(mr_c %>% filter(value == "GBR" & coef < 0))
df_phy5 = data.frame(sample = "CHN VS GBR", y = c(gbr_r,gbr_d), group = c("r","d"))

nld_r = nrow(mr_c %>% filter(value == "NLD" & coef > 0))
nld_d = nrow(mr_c %>% filter(value == "NLD" & coef < 0))
df_phy6 = data.frame(sample = "CHN VS NLD", y = c(nld_r,nld_d), group = c("r","d"))
df_phy = rbind(df_phy1,df_phy2,df_phy3,df_phy5,df_phy6)

colnames(df_phy)[2] = "number"
ggplot(df_phy, aes(x = sample, y = number, fill = group)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_ipsum() + scale_fill_brewer(palette="BuPu") +  labs(x = "", y = "Differential abundance features") +
  scale_fill_discrete(labels = c("Deleted", "Enriched")) + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(family = "Times", size = 15),
        legend.text = element_text(family = "Times", size = 14)
  ) + labs(fill = "")


tax_t = unique(mr_c[,1])
tax_t = sapply(tax_t,function(x) unlist(strsplit(x, split = "\\."))[7])
tax_t = sapply(tax_t,function(x) unlist(strsplit(x, split = "__"))[2])
tax_t = sapply(tax_t,function(x) gsub("_", " ", x))

con_cr_t = subset_taxa(con_cr, Species %in% tax_t)
con_cr_t = subset_samples(con_cr_t, country == "CHN")

con_cr_tree = phy_tree(con_cr_t)
con_cr_tree$tip.label = sapply(con_cr_tree$tip.label,function(x) unlist(strsplit(x, split = "\\|"))[7])
con_cr_tree$tip.label = sapply(con_cr_tree$tip.label,function(x) unlist(strsplit(x, split = "__"))[2])
con_cr_tree$tip.label = sapply(con_cr_tree$tip.label,function(x) gsub("_", " ", x))
phy_tree(con_cr_t) = con_cr_tree

#p = ggtree(con_cr_tree, layout="fan", open.angle=10, size=0.5) 
#p1 = p +  geom_tiplab(size=3)

p = ggtree(con_cr_t, layout="fan", open.angle=10) + 
  geom_tippoint(mapping=aes(color=Phylum), 
                size=1.5,
                show.legend=FALSE) 

getPalette = colorRampPalette(brewer.pal(12, "Paired"))

taxd_ring = mr_c %>% mutate(type = ifelse(coef > 0, "Increased", "Decreased")) %>% select(feature, coef, type, value)
taxd_ring[ ,1] = sapply(taxd_ring[ ,1],function(x) gsub("\\.", "\\|", x))
colnames(taxd_ring)[4] = "Country"
taxd_ring = taxd_ring %>% mutate(alpha = abs(coef/max(abs(coef))))
taxd_ring$Type = paste(taxd_ring$Country, taxd_ring$type, sep=" ")
taxd_ring$alpha = 0.01
mcon = psmelt(con_cr_t) %>% select(OTU, val=Abundance)
mcon$val = mcon$val / 100
p1 = p + new_scale_fill() + 
  geom_fruit(
  data=taxd_ring,
  geom=geom_tile,
  mapping=aes(y=feature, x=Country, alpha = alpha, fill=Type),
  pwidth=0.15,color = "grey50", offset = 0.04,size = 0.02,
  axis.params=list(
    axis="x", # add axis text of the layer.
    text.angle=-90, # the text angle of x-axis.
    hjust=0,  # adjust the horizontal position of text of axis.
    text.size = 1.5
  )
)  +  scale_fill_manual(
  #values=c("#b22222", "#005500", "#0000be", "#9f1f9f", "#dfac03"),
  values = getPalette(10),
  guide=guide_legend(keywidth=0.5, keyheight=0.5, order=4)
) + scale_alpha_continuous(
  range=c(0, 0.4), # the range of alpha
  #guide=guide_legend(keywidth=0.5, keyheight=0.5, order=5)
) + guides(alpha = "none")

p2 <- p1 + new_scale_fill() + 
  geom_fruit(
    data=mcon,
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=val,
      group=label,
      fill=Phylum,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 2,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) +  theme(#legend.position=c(0.96, 0.5), # the position of legend.
    legend.title=element_text(size=6), # the title size of legend.
    legend.text=element_text(size=5), # the text size of legend.
    legend.spacing.y = unit(0.02, "cm")  # the distance of legends (y orientation).
  ) 
#======================================================================
con_cr = subset_samples(con_cr,gender %in% c("female","male"))
con_cr = subset_samples(con_cr,!is.na(age))

df_otu = as.data.frame(otu_table(con_cr))
rownames(df_otu) = sapply(rownames(df_otu),function(x) unlist(strsplit(x, split = "\\|"))[7])
dfmd = country_md_cr_t
dfmd = dfmd %>% mutate(age_fl = case_when(age <= 45 ~ "n1",
                                          age > 45 ~ "n2"))
dfmd = dfmd %>% mutate(country_chn = ifelse(country == "CHN", "CHN", "R"))

dfmd = as.data.frame(dfmd)
rownames(dfmd) = dfmd$sample_id

fit = Maaslin2(input_data = df_otu, input_metadata = dfmd, output = "ma_output",
               normalization = "NONE", transform = "LOG",
               fixed_effects = c("country_chn", "gender","age_fl"),
               random_effects = "study_name",
               reference = "country_chn,CHN",
               min_prevalence = 0.1,
               min_abundance = 0.0001)
fr = fit$results
mr_bm_chn = fr %>% filter(metadata == "country_chn" & qval < 0.1)
#======================================================================
con_cr = subset_samples(con_cr,gender %in% c("female","male"))
con_cr = subset_samples(con_cr,!is.na(age))

df_otu = as.data.frame(otu_table(con_cr))
rownames(df_otu) = sapply(rownames(df_otu),function(x) unlist(strsplit(x, split = "\\|"))[7])
dfmd = country_md_cr_t
dfmd = dfmd %>% mutate(age_fl = case_when(age <= 45 ~ "n1",
                                          age > 45 ~ "n2"))
dfmd = dfmd %>% mutate(country_chn = ifelse(country == "USA", "USA", "R"))

dfmd = as.data.frame(dfmd)
rownames(dfmd) = dfmd$sample_id

fit = Maaslin2(input_data = df_otu, input_metadata = dfmd, output = "ma_output",
               normalization = "NONE", transform = "LOG",
               fixed_effects = c("country_chn", "gender","age_fl"),
               random_effects = "study_name",
               reference = "country_chn,USA",
               min_prevalence = 0.1,
               min_abundance = 0.0001)
fr = fit$results
mr_bm_usa = fr %>% filter(metadata == "country_chn" & qval < 0.1)
#======================================================================
con_cr = subset_samples(con_cr,gender %in% c("female","male"))
con_cr = subset_samples(con_cr,!is.na(age))

df_otu = as.data.frame(otu_table(con_cr))
rownames(df_otu) = sapply(rownames(df_otu),function(x) unlist(strsplit(x, split = "\\|"))[7])
dfmd = country_md_cr_t
dfmd = dfmd %>% mutate(age_fl = case_when(age <= 45 ~ "n1",
                                          age > 45 ~ "n2"))
dfmd = dfmd %>% mutate(country_chn = ifelse(country == "IND", "IND", "R"))

dfmd = as.data.frame(dfmd)
rownames(dfmd) = dfmd$sample_id

fit = Maaslin2(input_data = df_otu, input_metadata = dfmd, output = "ma_output",
               normalization = "NONE", transform = "LOG",
               fixed_effects = c("country_chn", "gender","age_fl"),
               random_effects = "study_name",
               reference = "country_chn,IND",
               min_prevalence = 0.1,
               min_abundance = 0.0001)
fr = fit$results
mr_bm_ind = fr %>% filter(metadata == "country_chn" & qval < 0.1)
#======================================================================
con_cr = subset_samples(con_cr,gender %in% c("female","male"))
con_cr = subset_samples(con_cr,!is.na(age))

df_otu = as.data.frame(otu_table(con_cr))
rownames(df_otu) = sapply(rownames(df_otu),function(x) unlist(strsplit(x, split = "\\|"))[7])
dfmd = country_md_cr_t
dfmd = dfmd %>% mutate(age_fl = case_when(age <= 45 ~ "n1",
                                          age > 45 ~ "n2"))
dfmd = dfmd %>% mutate(country_chn = ifelse(country == "MDG", "MDG", "R"))

dfmd = as.data.frame(dfmd)
rownames(dfmd) = dfmd$sample_id

fit = Maaslin2(input_data = df_otu, input_metadata = dfmd, output = "ma_output",
               normalization = "NONE", transform = "LOG",
               fixed_effects = c("country_chn", "gender","age_fl"),
               random_effects = "study_name",
               reference = "country_chn,MDG",
               min_prevalence = 0.1,
               min_abundance = 0.0001)
fr = fit$results
mr_bm_mdg = fr %>% filter(metadata == "country_chn" & qval < 0.1)
#======================================================================
con_cr = subset_samples(con_cr,gender %in% c("female","male"))
con_cr = subset_samples(con_cr,!is.na(age))

df_otu = as.data.frame(otu_table(con_cr))
rownames(df_otu) = sapply(rownames(df_otu),function(x) unlist(strsplit(x, split = "\\|"))[7])
dfmd = country_md_cr_t
dfmd = dfmd %>% mutate(age_fl = case_when(age <= 45 ~ "n1",
                                          age > 45 ~ "n2"))
dfmd = dfmd %>% mutate(country_chn = ifelse(country == "GBR", "GBR", "R"))

dfmd = as.data.frame(dfmd)
rownames(dfmd) = dfmd$sample_id

fit = Maaslin2(input_data = df_otu, input_metadata = dfmd, output = "ma_output",
               normalization = "NONE", transform = "LOG",
               fixed_effects = c("country_chn", "gender","age_fl"),
               random_effects = "study_name",
               reference = "country_chn,GBR",
               min_prevalence = 0.1,
               min_abundance = 0.0001)
fr = fit$results
mr_bm_gbr = fr %>% filter(metadata == "country_chn" & qval < 0.1)


