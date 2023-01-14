library(curatedMetagenomicData)
library(phyloseq)
library(mia)
library(dplyr)
library(microeco)
library(file2meco)
library(RColorBrewer)
library(hrbrthemes)
library(ggplot2)
library(magrittr)
library(scater)
library(DirichletMultinomial)
library(Maaslin2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
#======================================================================
setwd("e:/ma16/")

load(file = "country_sample.csv")

country_crt = country_md_cr_t %>%
  returnSamples("relative_abundance")

country_crt1 = country_md_cr_t %>%
  returnSamples("relative_abundance", counts = T)

#======================================================================
con_cr = makePhyloseqFromTreeSummarizedExperiment(country_crt, abund_values = "relative_abundance")
m = phyloseq2meco(con_cr)
m1 = clone(m)
m1$sample_table = subset(m1$sample_table, country == "CHN")
m1$tidy_dataset()
m1$cal_abund()
t = trans_abund$new(dataset = m1, taxrank = "Phylum", groupmean = "country")
chn_p_df = t$abund_data
chn_p_df$Abundance = signif(chn_p_df$Abundance,3)

t = trans_abund$new(dataset = m1, taxrank = "Genus", groupmean = "country")
chn_g_df = t$abund_data
chn_g_df$Abundance = signif(chn_g_df$Abundance,3)

m = phyloseq2meco(con_cr)
m1 = clone(m)
m1$sample_table = subset(m1$sample_table, country == "CHN")
m1$tax_table = as.data.frame(tax_table(con_cr))
m1$tax_table$Species = paste("s__", m1$tax_table$Species, sep = "")
m1$cal_abund()
t = trans_abund$new(dataset = m1, taxrank = "Species", groupmean = "country")
chn_s_df = t$abund_data
chn_s_df$Abundance = signif(chn_s_df$Abundance,3)
#======================================================================
con_cr = subset_samples(con_cr,gender %in% c("female","male"))
con_cr = subset_samples(con_cr,!is.na(age))
con_cr_pc = tax_glom(con_cr, taxrank = rank_names(con_cr)[2])
tdf = as.data.frame(tax_table(con_cr_pc))
dfotu = as.data.frame(otu_table(con_cr_pc))
dfotu = as.data.frame(t(dfotu))
dfotu = cbind(dfotu, as.data.frame(sample_data(con_cr_pc))[ ,c("study_name", "country", "gender", "age")])
colnames(dfotu)[21] = "Country"

fit = Maaslin2(input_data = dfotu[,"bfr",drop = F], input_metadata = dfotu[ ,c("country", "gender", "age")], output = "ma_output",
               normalization = "NONE", transform = "NONE",
               fixed_effects = c("country", "gender","age"),
               random_effects = "study_name",
               reference = "country,CHN",
               min_prevalence = 0,
               min_abundance = 0)
fr = fit$results

ggplot(dfotu, aes(x = Country, y = bfr, fill = Country)) + geom_boxplot(alpha=0.3) +
  scale_fill_brewer(palette="Paired") + ylim(0,50) + ylab("B/F ratio") + theme_ipsum() + 
  theme(legend.position="none") + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15))
#======================================================================
con_cr = filter_taxa(con_cr, function(x) mean(x) > 0, TRUE)
con_cr_pc = tax_glom(con_cr, taxrank = rank_names(con_cr)[2])
dfotu = as.data.frame(otu_table(con_cr_pc))

con_cr_pc = tax_glom(con_cr, taxrank = rank_names(con_cr)[6])
con_cr_pc = subset_taxa(con_cr_pc, Phylum=="Euryarchaeota")
taxdf = as.data.frame(tax_table(con_cr_pc))
dfotu = as.data.frame(otu_table(con_cr_pc))

con_cr_pc = subset_taxa(con_cr, Phylum=="Euryarchaeota")
taxdf = as.data.frame(tax_table(con_cr_pc))
dfotu = as.data.frame(otu_table(con_cr_pc))[9,]
dfotu = t(dfotu)
colnames(dfotu) = "Methanobrevibacter_smithii"
dfotu = cbind(dfotu, as.data.frame(sample_data(con_cr_pc))[ ,"country"])
dfotu[,1] = ifelse(dfotu[,1]>0,1,0)
df = full_join(dfotu %>% group_by(country) %>% summarise(n1 = sum(Methanobrevibacter_smithii)),
               dfotu %>% group_by(country) %>% count())
df$n =  df$n-df$n1
pairwiseNominalIndependence(df,chisq = T, fisher = F, gtest = F)
#======================================================================
m = phyloseq2meco(con_cr)
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = 10, groupmean = "country")

t$plot_bar(others_color = "grey70", legend_text_italic = FALSE) + theme_ipsum() +
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
        ) +  theme(legend.title = element_text(size = 15, family = "Times"),
             legend.text = element_text(size = 14, family = "Times")
)

df_t_p = t$abund_data
df_t_p$Abundance = signif(df_t_p$Abundance,3)
write.csv(df_t_p, file = "tax_df_p.tsv", row.names=F, quote = F)

con_cr_p = tax_glom(con_cr, rank_names(con_cr)[2])
df1 = as.data.frame(otu_table(con_cr_p))
df1[df1 > 0] = 1
dfc = apply(df1,1,function(x) sum(x)/ncol(df1))
dfn = apply(df1,1,function(x) sum(x))
dfp = data.frame(dfc,dfn)
rownames(dfp) = rownames(df1) = sapply(rownames(df1),function(x) unlist(strsplit(x, split = "\\|"))[2])
df1$dfc = dfc
df1$dfn = dfn

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
m = phyloseq2meco(con_cr)
m$tax_table = as.data.frame(tax_table(con_cr))
m$tax_table$Genus = paste("g__", m$tax_table$Genus, sep = "")
m$cal_abund()
t1 = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = 20, groupmean = "country")

t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE, use_colors = getPalette(20)) + theme_ipsum() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
  ) + theme(legend.title = element_text(size = 15, family = "Times"),
            legend.text = element_text(size = 14, family = "Times")
)


df_t_g = t1$abund_data
df_t_g$Abundance = signif(df_t_g$Abundance,3)
write.csv(df_t_g, file = "tax_df_g.tsv", row.names=F, quote = F)

df1 = df_t_g %>% filter(Sample == "CHN") 
df2 = df_t_g %>% filter(Sample == "USA") 
df3 = df_t_g %>% filter(Sample == "GBR") 
df4 = df_t_g %>% filter(Sample == "NLD") 
df5 = df_t_g %>% filter(Sample == "IND")
df6 = df_t_g %>% filter(Sample == "MDG")

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
m = phyloseq2meco(con_cr)
m$tax_table = as.data.frame(tax_table(con_cr))
m$tax_table$Species = paste("s__", m$tax_table$Species, sep = "")
m$cal_abund()
t1 = trans_abund$new(dataset = m, taxrank = "Species", ntaxa = 30, groupmean = "country")

t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE, use_colors = getPalette(30)) + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
  )
+ theme(legend.title = element_text(size = 15, family = "Times"),
            legend.text = element_text(size = 14, family = "Times"))

df_t_s = t1$abund_data
df_t_s$Abundance = signif(df_t_s$Abundance,3)
write.csv(df_t_s, file = "tax_df_g.tsv", row.names=F, quote = F)
t1 = trans_abund$new(dataset = m, taxrank = "Species", ntaxa = 30)
t1$plot_box(group = "country")

con_crs = filter_taxa(con_cr, function(x) mean(x) > 0.1, TRUE)
m = phyloseq2meco(con_crs)
m$cal_abund()
t1 = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = 20, groupmean = "country")
df_t_s = t1$abund_data
df_w = tidyr::spread(df_t_s, key = "Taxonomy", value = "Abundance", fill = 0)
rownames(df_w) = df_w[,1]
df_w = df_w[,-1]
df_w = t(df_w)
df_w = df_w/100
pheatmap::pheatmap(df_w, cluster_rows = F,cellwidth = 15, border_color = NA,legend = T, 
         clustering_distance_rows = "correlation", show_rownames = F, fontsize = 18)

#======================================================================
library(forcats)
con_cr_pc = tax_glom(con_cr, taxrank = rank_names(con_cr)[2])

con_cr_chn = subset_samples(con_cr_pc, country == "CHN")
con_cr_chn_share90 = filter_taxa(con_cr_chn,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_chn), TRUE)
df_chn_s90 = as.data.frame(tax_table(con_cr_chn_share90))
con_p90 = subset_taxa(con_cr_pc, Phylum %in% df_chn_s90$Phylum)
m = phyloseq2meco(con_p90)
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = 5)
df1 = t$abund_data %>% filter(country == "CHN")
chn_cog = full_join(df1 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df1 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))

con_cr_usa = subset_samples(con_cr_pc, country == "USA")
con_cr_usa_share90 = filter_taxa(con_cr_usa,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_usa), TRUE)
df_usa_s90 = as.data.frame(tax_table(con_cr_usa_share90))
con_p90 = subset_taxa(con_cr_pc, Phylum %in% df_usa_s90$Phylum)
m = phyloseq2meco(con_p90)
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = 5)
df2 = t$abund_data %>% filter(country == "USA") 
usa_cog = full_join(df2 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df2 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))

con_cr_ind = subset_samples(con_cr_pc, country == "IND")
con_cr_ind_share90 = filter_taxa(con_cr_ind,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_ind), TRUE)
df_ind_s90 = as.data.frame(tax_table(con_cr_ind_share90))
con_p90 = subset_taxa(con_cr_pc, Phylum %in% df_ind_s90$Phylum)
m = phyloseq2meco(con_p90)
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = 5)
df3 = t$abund_data %>% filter(country == "IND") 
ind_cog = full_join(df3 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df3 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))

con_cr_mdg = subset_samples(con_cr_pc, country == "MDG")
con_cr_mdg_share90 = filter_taxa(con_cr_mdg,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_mdg), TRUE)
df_mdg_s90 = as.data.frame(tax_table(con_cr_mdg_share90))
con_p90 = subset_taxa(con_cr_pc, Phylum %in% df_mdg_s90$Phylum)
m = phyloseq2meco(con_p90)
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = 5)
df4 = t$abund_data %>% filter(country == "MDG") 
mdg_cog = full_join(df4 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df4 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
con_cr_gbr = subset_samples(con_cr_pc, country == "GBR")
con_cr_gbr_share90 = filter_taxa(con_cr_gbr,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_gbr), TRUE)
df_gbr_s90 = as.data.frame(tax_table(con_cr_gbr_share90))
con_p90 = subset_taxa(con_cr_pc, Phylum %in% df_gbr_s90$Phylum)
m = phyloseq2meco(con_p90)
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = 5)
df5 = t$abund_data %>% filter(country == "GBR") 
gbr_cog = full_join(df5 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df5 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
con_cr_nld = subset_samples(con_cr_pc, country == "NLD")
con_cr_nld_share90 = filter_taxa(con_cr_nld,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_nld), TRUE)
df_nld_s90 = as.data.frame(tax_table(con_cr_nld_share90))
con_p90 = subset_taxa(con_cr_pc, Phylum %in% df_nld_s90$Phylum)
m = phyloseq2meco(con_p90)
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = 5)
df6 = t$abund_data %>% filter(country == "NLD") 
nld_cog = full_join(df6 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df6 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
df = rbind(df1,df2,df3,df4,df5,df6)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
df %>% mutate(Sample = fct_reorder(Sample, Abundance)) %>%
  ggplot(aes(x = Sample, y = Abundance, color = Taxonomy, fill = Taxonomy)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_bar(stat="identity", position="stack") + facet_grid(.~country, scales = "free_x", space = "free_x") +
  xlab("") + scale_color_manual(values = brewer.pal(12, "Paired")) + scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
  )

#======================================================================
con_cr_pc = tax_glom(con_cr, taxrank = rank_names(con_cr)[6])

con_cr_chn = subset_samples(con_cr_pc, country == "CHN")
con_cr_chn_share90 = filter_taxa(con_cr_chn,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_chn), TRUE)
df_chn_s90 = as.data.frame(tax_table(con_cr_chn_share90))
con_p90 = subset_taxa(con_cr_pc, Genus %in% df_chn_s90$Genus)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Genus = paste("g__", m$tax_table$Genus, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = nrow(df_chn_s90))
df1 = t$abund_data %>% filter(country == "CHN")
chn_cog = full_join(df1 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
   df1 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                            max = max(Abundance)))

con_cr_usa = subset_samples(con_cr_pc, country == "USA")
con_cr_usa_share90 = filter_taxa(con_cr_usa,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_usa), TRUE)
df_usa_s90 = as.data.frame(tax_table(con_cr_usa_share90))
con_p90 = subset_taxa(con_cr_pc, Genus %in% df_usa_s90$Genus)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Genus = paste("g__", m$tax_table$Genus, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = nrow(df_usa_s90))
df2 = t$abund_data %>% filter(country == "USA") 
usa_cog = full_join(df2 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df2 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))

con_cr_ind = subset_samples(con_cr_pc, country == "IND")
con_cr_ind_share90 = filter_taxa(con_cr_ind,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_ind), TRUE)
df_ind_s90 = as.data.frame(tax_table(con_cr_ind_share90))
con_p90 = subset_taxa(con_cr_pc, Genus %in% df_ind_s90$Genus)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Genus = paste("g__", m$tax_table$Genus, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = nrow(df_ind_s90))
df3 = t$abund_data %>% filter(country == "IND") 
ind_cog = full_join(df3 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df3 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))

con_cr_mdg = subset_samples(con_cr_pc, country == "MDG")
con_cr_mdg_share90 = filter_taxa(con_cr_mdg,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_mdg), TRUE)
df_mdg_s90 = as.data.frame(tax_table(con_cr_mdg_share90))
con_p90 = subset_taxa(con_cr_pc, Genus %in% df_mdg_s90$Genus)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Genus = paste("g__", m$tax_table$Genus, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = nrow(df_mdg_s90))
df4 = t$abund_data %>% filter(country == "MDG") 
mdg_cog = full_join(df4 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df4 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
con_cr_gbr = subset_samples(con_cr_pc, country == "GBR")
con_cr_gbr_share90 = filter_taxa(con_cr_gbr,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_gbr), TRUE)
df_gbr_s90 = as.data.frame(tax_table(con_cr_gbr_share90))
con_p90 = subset_taxa(con_cr_pc, Genus %in% df_gbr_s90$Genus)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Genus = paste("g__", m$tax_table$Genus, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = nrow(df_gbr_s90))
df5 = t$abund_data %>% filter(country == "GBR") 
gbr_cog = full_join(df5 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df5 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
con_cr_nld = subset_samples(con_cr_pc, country == "NLD")
con_cr_nld_share90 = filter_taxa(con_cr_nld,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_nld), TRUE)
df_nld_s90 = as.data.frame(tax_table(con_cr_nld_share90))
con_p90 = subset_taxa(con_cr_pc, Genus %in% df_nld_s90$Genus)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Genus = paste("g__", m$tax_table$Genus, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = nrow(df_nld_s90))
df6 = t$abund_data %>% filter(country == "NLD") 
nld_cog = full_join(df6 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df6 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
df = rbind(df1,df2,df3,df4,df5,df6)
#con_cr_pc90 = filter_taxa(con_cr_pc,function(x) sum(x!=0) >= 0.90*nsamples(con_cr_pc), TRUE)
#df_s90 = as.data.frame(tax_table(con_cr_pc90))
#con_p90 = subset_taxa(con_cr_pc, Phylum %in% df_s90$Phylum)
#m = phyloseq2meco(con_p90)
#m$cal_abund()
#t = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = nrow(df_s90))
#df = t$abund_data
#p = t$plot_bar(others_color = "grey70", facet = "country", xtext_keep = FALSE, legend_text_italic = FALSE)
#p + geom_bar(aes(x = reorder(Sample, Abundance)))

#p = plot_bar(con_p90,fill="Phylum", facet_grid=~country) + theme(axis.text.x = element_blank())
getPalette = colorRampPalette(c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3")))
df %>% mutate(Sample = fct_reorder(Sample, Abundance)) %>%
 ggplot(aes(x = Sample, y = Abundance, color = Taxonomy, fill = Taxonomy)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_bar(stat="identity", position="stack") + facet_grid(.~country, scales = "free_x", space = "free_x") +
  xlab("") + scale_color_manual(values = getPalette(38)) + scale_fill_manual(values = getPalette(38)) + 
  theme(axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
  )


df_con = Reduce(union, list(df_usa_s90$Genus,df_ind_s90$Genus,df_mdg_s90$Genus,
                            df_gbr_s90$Genus,df_nld_s90$Genus))
df_u_chn = setdiff(df_chn_s90$Genus,df_con)

df_con_t = Reduce(union, list(df_chn_s90$Genus,df_usa_s90$Genus,df_ind_s90$Genus,
                              df_mdg_s90$Genus,df_gbr_s90$Genus,df_nld_s90$Genus))

chn_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_chn_s90$Genus, 1, 0)))
usa_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_usa_s90$Genus, 1, 0)))
ind_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_ind_s90$Genus, 1, 0)))
mdg_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_mdg_s90$Genus, 1, 0)))
gbr_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_gbr_s90$Genus, 1, 0)))
nld_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_nld_s90$Genus, 1, 0)))
core_df = data.frame(chn_ar,usa_ar,ind_ar,mdg_ar,gbr_ar,nld_ar)
colnames(core_df) = c("CHN","USA","IND","MDG","GBR","NLD")
pheatmap::pheatmap(core_df,cluster_rows = FALSE,cluster_cols = F, 
                   legend = FALSE, cellwidth = 15, color = c("#A6D854", "#FFD92F"),fontsize = 12)

#======================================================================
con_cr_pc = tax_glom(con_cr, taxrank = rank_names(con_cr)[7])

con_cr_chn = subset_samples(con_cr_pc, country == "CHN")
con_cr_chn_share90 = filter_taxa(con_cr_chn,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_chn), TRUE)
df_chn_s90 = as.data.frame(tax_table(con_cr_chn_share90))
con_p90 = subset_taxa(con_cr_pc, Species %in% df_chn_s90$Species)
#con_p90 = subset_samples(con_p90, country == "CHN")
#con_p90 = subset_taxa(con_p90, Genus != "unidentified")
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Species = paste("s__", m$tax_table$Species, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Species", ntaxa = nrow(df_s90))
df1 = t$abund_data %>% filter(country == "CHN")
chn_cog = full_join(df1 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df1 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))

con_cr_usa = subset_samples(con_cr_pc, country == "USA")
con_cr_usa_share90 = filter_taxa(con_cr_usa,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_usa), TRUE)
df_usa_s90 = as.data.frame(tax_table(con_cr_usa_share90))
con_p90 = subset_taxa(con_cr_pc, Species %in% df_usa_s90$Species)
m = phyloseq2meco(con_p90)
m$tidy_dataset()
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Species = paste("s__", m$tax_table$Species, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Species", ntaxa = nrow(df_s90))
df2 = t$abund_data %>% filter(country == "USA") 
usa_cog = full_join(df2 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df2 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))

con_cr_ind = subset_samples(con_cr_pc, country == "IND")
con_cr_ind_share90 = filter_taxa(con_cr_ind,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_ind), TRUE)
df_ind_s90 = as.data.frame(tax_table(con_cr_ind_share90))
con_p90 = subset_taxa(con_cr_pc, Species %in% df_ind_s90$Species)
m = phyloseq2meco(con_p90)
m$tidy_dataset()
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Species = paste("s__", m$tax_table$Species, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Species", ntaxa = nrow(df_s90))
df3 = t$abund_data %>% filter(country == "IND") 
ind_cog = full_join(df3 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df3 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))

con_cr_mdg = subset_samples(con_cr_pc, country == "MDG")
con_cr_mdg_share90 = filter_taxa(con_cr_mdg,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_mdg), TRUE)
df_mdg_s90 = as.data.frame(tax_table(con_cr_mdg_share90))
con_p90 = subset_taxa(con_cr_pc, Species %in% df_mdg_s90$Species)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Species = paste("s__", m$tax_table$Species, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Species", ntaxa = nrow(df_s90))
df4 = t$abund_data %>% filter(country == "MDG") 
mdg_cog = full_join(df4 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df4 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
con_cr_gbr = subset_samples(con_cr_pc, country == "GBR")
con_cr_gbr_share90 = filter_taxa(con_cr_gbr,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_gbr), TRUE)
df_gbr_s90 = as.data.frame(tax_table(con_cr_gbr_share90))
con_p90 = subset_taxa(con_cr_pc, Species %in% df_gbr_s90$Species)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Species = paste("s__", m$tax_table$Species, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Species", ntaxa = nrow(df_s90))
df5 = t$abund_data %>% filter(country == "GBR") 
gbr_cog = full_join(df5 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df5 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
con_cr_nld = subset_samples(con_cr_pc, country == "NLD")
con_cr_nld_share90 = filter_taxa(con_cr_nld,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_nld), TRUE)
df_nld_s90 = as.data.frame(tax_table(con_cr_nld_share90))
con_p90 = subset_taxa(con_cr_pc, Species %in% df_nld_s90$Species)
m = phyloseq2meco(con_p90)
m$tax_table = as.data.frame(tax_table(con_p90))
m$tax_table$Species = paste("s__", m$tax_table$Species, sep = "")
m$cal_abund()
t = trans_abund$new(dataset = m, taxrank = "Species", ntaxa = nrow(df_s90))
df6 = t$abund_data %>% filter(country == "NLD") 
nld_cog = full_join(df6 %>% group_by(Taxonomy) %>% summarise_each(function(x) min(x[x != 0]), min = Abundance),
                    df6 %>% group_by(Taxonomy) %>% summarise(mean = mean(Abundance), median = median(Abundance),
                                                             max = max(Abundance)))
df = rbind(df1,df2,df3,df4,df5,df6)
fs3d = df %>% mutate(Sample = fct_reorder(Sample, Abundance)) %>%
  ggplot(aes(x = Sample, y = Abundance, color = Taxonomy, fill = Taxonomy)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_bar(stat="identity", position="stack") + facet_grid(.~country, scales = "free_x", space = "free_x") +
  xlab("") + scale_color_manual(values = getPalette(38)) + scale_fill_manual(values = getPalette(38)) + 
  theme(axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
  )
#======================================================================
m = phyloseq2meco(con_cr1)
m1 = clone(m)
m1$sample_sums() %>% range
m1$rarefy_samples(sample.size = 701444)
m1$cal_alphadiv(measures = c("Chao1","Shannon","Simpson", "Observed"))
t1 = trans_alpha$new(dataset = m1, group = "country") 
#t1$plot_alpha(measure = "Shannon", group = "country", color_values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum()


zddf = m1$alpha_diversity
mddf = m1$sample_table

tot_df = cbind(zddf,mddf)

df_otu = tot_df[, 6, drop = F]
fit = Maaslin2(input_data = df_otu, input_metadata = tot_df, output = "ma_output",
               normalization = "NONE", transform = "NONE",
               fixed_effects = c("country", "gender","age"),
               random_effects = "study_name",
               reference = "country,CHN",
               min_prevalence = 0,
               min_abundance = 0)
fr = fit$results

gdf = mddf %>% group_by(country) %>% summarise(m = sum(gender == "male"), fm = sum(gender == "female"))

mtot_df = tot_df %>% group_by(country) %>% summarise(shn = mean(Shannon), chao = mean(Chao1), 
                                                     sim = mean(Simpson), obs = mean(Observed))
#plot alpha diversity
fs4a = ggplot(tot_df, aes(x = country, y = Shannon, color = country)) + geom_boxplot() + scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")) + labs(color = "Country") + theme_ipsum() +
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
  ) + theme(legend.title = element_text(family = "Times"),
            legend.text = element_text(family = "Times"))


ggplot(tot_df, aes(x = country, y = Chao1, color = country)) + geom_boxplot() + scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum()
ggplot(tot_df, aes(x = country, y = Simpson, color = country)) + geom_boxplot() + scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")) + labs(color = "Country") + theme_ipsum() 
ggplot(tot_df, aes(x = country, y = Observed, color = country)) + geom_boxplot() + scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum()

tot_df_f = tot_df %>% filter(gender == "female")
tot_df_m = tot_df %>% filter(gender == "male")

tot_df_f$country = factor(tot_df_f$country, levels = c("GBR","NLD","USA","CHN","MDG","IND"), ordered = T)
tot_df_m$country = factor(tot_df_m$country, levels = c("GBR","NLD","USA","CHN","MDG","IND"), ordered = T)


df1_chn = tot_df_f %>% filter(country == "CHN") 
df1_usa = tot_df_f %>% filter(country == "USA")
df1_ind = tot_df_f %>% filter(country == "IND")
df1_mdg = tot_df_f %>% filter(country == "MDG")
df1_nld = tot_df_f %>% filter(country == "NLD")
df1_gbr = tot_df_f %>% filter(country == "GBR")

ggplot(tot_df_f, aes(x = country, y = Simpson, color = country)) + geom_boxplot() + scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum()
ggplot(tot_df_m, aes(x = country, y = Chao1, color = country)) + geom_boxplot() + scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum()

m = phyloseq2meco(con_cr1)
m1 = clone(m)
m1$sample_table = subset(m1$sample_table, gender == "female")
m1$tidy_dataset()
m1$sample_sums() %>% range
m1$rarefy_samples(sample.size = 3050870)
m1$cal_alphadiv(measures = c("Chao1","Shannon","Simpson", "Observed"))
t1 = trans_alpha$new(dataset = m1, group = "country") 
t1$plot_alpha(measure = "Observed", group = "country",color_values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum()

m2 = clone(m)
m2$sample_table = subset(m2$sample_table, gender == "female")
m2$tidy_dataset()
m2$sample_sums() %>% range
m2$rarefy_samples(sample.size = 3050870)
m2$cal_alphadiv(measures = c("Chao1","Shannon","Simpson", "Observed"))
t2 = trans_alpha$new(dataset = m1, group = "country") 
t2$plot_alpha(measure = "Observed", group = "country",color_values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum()


dfz = m1$alpha_diversity$Shannon
dfd = m1$sample_table$age
ks.test(dfz,"pnorm")
ks.test(dfd,"pnorm")
cor.test(dfz,dfd, method = "spearman")

dfg = m1$sample_table$gender


ggplot(dfcor, aes(x = age, y = Shannon, color = country)) + geom_point(alpha = 0.5) + scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
  )

dfg = data.frame(g = m1$sample_table$gender, 
                 y1 = m1$alpha_diversity$Shannon,
                 c = m1$sample_table$country)

gm = dfg %>% filter(g == "male")
gf = dfg %>% filter(g == "female")

ggplot(dfcor, aes(x = age, y = country, fill = country)) +
  geom_density_ridges(adjust=1.5, alpha=.5) +
  theme_ridges() +
  theme(legend.position = "none")

m = phyloseq2meco(con_cr)
m$sample_table = subset(m$sample_table, gender %in% c("female","male"))
m$sample_table <- m$sample_table[!is.na(m$sample_table$age),]
m$tidy_dataset()
m$cal_betadiv(unifrac = F)
t1 = trans_beta$new(dataset = m, group = "country", measure = "bray")
t1$cal_ordination(ordination = "PCoA")

f3a = t1$plot_ordination(plot_color = "country", plot_point_alpha = 0.5) + theme_ipsum() + labs(color = "Country") +
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13),
        legend.title = element_text(family = "Times", size = 14),
        legend.text = element_text(family = "Times", size = 13)
  ) + theme(legend.title = element_text(family = "Times", size = 14),
            legend.text = element_text(family = "Times", size = 13))
save(f3a, file = "f3a.rda")
load(file = "f3a.rda")

t1$cal_manova(cal_manova_set = "country + age + gender + study_name",permutations = 999)
dfp = t1$res_manova

df_ord = t1$res_ordination$scores
mp = clone(m)
mp = mp$merge_taxa(taxa = "Phylum")
dfotu = mp$otu_table
p = as.numeric(dfotu[1,])
cor.test(p,df_ord$PCo1, method = "spearman")
cor.test(p,df_ord$PCo2, method = "spearman")

mp = clone(m)
mp = mp$merge_taxa(taxa = "Genus")
dfotu = mp$otu_table
dftax = mp$tax_table
p = as.numeric(dfotu[9,])
cor.test(p,df_ord$PCo1,method = "spearman")
cor.test(p,df_ord$PCo2,method = "spearman")
dfp = data.frame(t1 = as.numeric(dfotu[1,]),
                 t2 = as.numeric(dfotu[2,]),
                 t3 = as.numeric(dfotu[3,]),
                 t4 = as.numeric(dfotu[4,]),
                 t5 = as.numeric(dfotu[5,]),
                 t6 = as.numeric(dfotu[6,]),
                 t7 = as.numeric(dfotu[7,]),
                 t8 = as.numeric(dfotu[8,]),
                 t9 = as.numeric(dfotu[9,]),
                 Pco1 = df_ord$PCo1, Pco2 = df_ord$PCo2)

pc11 = ggplot(dfp, aes(x=Pco1, y=t1)) +
  geom_point(alpha = 1/10) + geom_smooth(method=lm) + ylab(dftax[1,6]) + theme_ipsum() +
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13)
  )

pc13 = ggplot(dfp, aes(x=Pco1, y=t3)) +
  geom_point(alpha = 1/10) + geom_smooth(method=lm) + ylab(dftax[3,6]) + theme_ipsum() +
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13)
  )

pc19 = ggplot(dfp, aes(x=Pco1, y=t9)) +
  geom_point(alpha = 1/10) + geom_smooth(method=lm) + ylab(dftax[9,6]) + theme_ipsum() + 
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13)
  )
pc21 = ggplot(dfp, aes(x=Pco2, y=t1)) +
  geom_point(alpha = 1/10) + geom_smooth(method=lm) + ylab(dftax[1,6]) + theme_ipsum() +
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13)
  )
pc22 = ggplot(dfp, aes(x=Pco2, y=t2)) +
  geom_point(alpha = 1/10) + geom_smooth(method=lm) + ylab(dftax[2,6]) + theme_ipsum() + 
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13)
  )
pc27 = ggplot(dfp, aes(x=Pco2, y=t7)) +
  geom_point(alpha = 1/10) + geom_smooth(method=lm) + ylab(dftax[7,6]) + theme_ipsum() + 
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13)
  )
library(patchwork)
pc11+pc13+pc19+pc21+pc22+pc27


df_ord = t1$res_ordination$scores
ggplot(df_ord, aes(x=country, y=PCo1, fill=country)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Paired") + coord_flip()+theme_ipsum() + 
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
  ) + xlab("Country") + labs(fill = "Country")
+ theme(legend.title = element_text(family = "Times", size = 14),
            legend.text = element_text(family = "Times", size = 13))


ggplot(df_ord, aes(x=country, y=PCo2, fill=country)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Paired") + coord_flip() + theme_ipsum() + 
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
  ) + xlab("Country") + labs(fill = "Country")
+ theme(legend.title = element_text(family = "Times", size = 14),
            legend.text = element_text(family = "Times", size = 13))


#===============================================================================
#mmdg
country_mdg = filter(sampleMetadata, disease=="healthy") %>% 
  filter(body_site =="stool")  %>%
  filter(country %in% c("MDG", "IND"))
country_mmdg = country_mdg %>%
  returnSamples("relative_abundance")
country_mmdg1 = country_mdg %>%
  returnSamples("relative_abundance", counts = T)

pmdg = makePhyloseqFromTreeSummarizedExperiment(country_mmdg, abund_values = "relative_abundance")
pmdg1 = makePhyloseqFromTreeSummarizedExperiment(country_mmdg1, abund_values = "relative_abundance")

mdgr = subset_samples(pmdg, country == "MDG")
mdg = subset_samples(pmdg1, country == "MDG")

summary(sample_sums(mdg))
set.seed(1)
ps.rarefied = rarefy_even_depth(mdg, rngseed=1, sample.size=min(sample_sums(mdg)), replace=F)
global_m = as.matrix(as.data.frame(otu_table(mdg)))
m_rare = rrarefy(t(global_m), min(colSums(global_m)))
shannon_div = diversity(m_rare, index = "shannon")
plot_richness(ps.rarefied, x="country", color="country", measures=c("Observed")) + geom_boxplot()
plot_richness(ps.rarefied, x="country", color="country", measures=c("Shannon")) + geom_boxplot()

global_mr = as.matrix(as.data.frame(otu_table(mdgr)))
bcd = vegdist(t(global_mr), method = "bray")
plot(as.vector(bcd))
dist = phyloseq::distance(mdgr, method="bray")

m = phyloseq2meco(mdgr)
m$tidy_dataset()
m$cal_betadiv(unifrac = F)
t1 = trans_beta$new(dataset = m, group = "country", measure = "bray")
t1$cal_group_distance()
t1$plot_group_distance(color_values = RColorBrewer::brewer.pal(12, "Paired"))
#===============================================================================
t1$cal_group_distance()
bdf = t1$res_group_distance
fs5a = t1$plot_group_distance(color_values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum() + 
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
save(fs5a, file = "fs5a.rda")
load(file = "fs5a.rda")

mbdf = bdf %>% group_by(country) %>% summarise(mb = mean(value))
mt_df = full_join(mtot_df, mbdf, by = "country")

df1 = df %>% filter(country %in% c("CHN vs MDG","CHN vs NLD","GBR vs CHN","IND vs CHN","USA vs CHN"))
t1$res_group_distance = df1
fs5b = t1$plot_group_distance(color_values = RColorBrewer::brewer.pal(12, "Paired")) + theme_ipsum() + labs(color = "Country") + 
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15))
save(fs5b, file = "fs5b.rda")
load(file = "fs5b.rda")

df_c =  df %>% filter(country %in% c("GBR vs CHN","USA vs CHN"))
wilcox.test(value~country, df_c)

m1 = clone(m)
m1$sample_table <- subset(m1$sample_table, country == "NLD")
m1$tidy_dataset()

m1 = clone(m)
m1$sample_table <- subset(m1$sample_table, country %in% c("CHN","IND"))
m1$sample_table <- subset(m1$sample_table, gender %in% c("female","male"))
m1$sample_table <- m1$sample_table[!is.na(m1$sample_table$age),]
m1$tidy_dataset()
m1$cal_betadiv(unifrac = F)
t11 = trans_beta$new(dataset = m1, group = "country", measure = "bray")
t11$cal_manova(cal_manova_set = "country + age + gender + study_name",permutations = 999)
t11$res_manova

m2 = clone(m)
m2$sample_table <- subset(m2$sample_table, country %in% c("CHN","MDG"))
m2$sample_table <- subset(m2$sample_table, gender %in% c("female","male"))
m2$sample_table <- m2$sample_table[!is.na(m2$sample_table$age),]
m2$tidy_dataset()
m2$cal_betadiv(unifrac = F)
t2 = trans_beta$new(dataset = m2, group = "country", measure = "bray")
t2$cal_manova(cal_manova_set = "country + age + gender + study_name",permutations = 999)
t2$res_manova

m3 = clone(m)
m3$sample_table <- subset(m3$sample_table, country %in% c("CHN","NLD"))
m3$sample_table <- subset(m3$sample_table, gender %in% c("female","male"))
m3$sample_table <- m3$sample_table[!is.na(m3$sample_table$age),]
m3$tidy_dataset()
m3$cal_betadiv(unifrac = F)
t3 = trans_beta$new(dataset = m3, group = "country", measure = "bray")
t3$cal_manova(cal_manova_set = "country + age + gender + study_name",permutations = 999)
t3$res_manova

m4 = clone(m)
m4$sample_table <- subset(m4$sample_table, country %in% c("CHN","GBR"))
m4$sample_table <- subset(m4$sample_table, gender %in% c("female","male"))
m4$sample_table <- m4$sample_table[!is.na(m4$sample_table$age),]
m4$tidy_dataset()
m4$cal_betadiv(unifrac = F)
t4 = trans_beta$new(dataset = m4, group = "country", measure = "bray")
t4$cal_manova(cal_manova_set = "country + age + gender + study_name",permutations = 999)
t4$res_manova

m5 = clone(m)
m5$sample_table <- subset(m5$sample_table, country %in% c("CHN","USA"))
m5$sample_table <- subset(m5$sample_table, gender %in% c("female","male"))
m5$sample_table <- m5$sample_table[!is.na(m5$sample_table$age),]
m5$tidy_dataset()
m5$cal_betadiv(unifrac = F)
t5 = trans_beta$new(dataset = m5, group = "country", measure = "bray")
t5$cal_manova(cal_manova_set = "country + age + gender + study_name",permutations = 999)
t5$res_manova

df_p = data.frame(country = c("USA","IND","MDG","NLD","GBR"), f = c(18, 43, 58, 221, 69))
df_p %>% mutate(country = reorder(country, f)) %>% ggplot(aes(x = country, y = f, fill = country)) + geom_bar(stat = "identity", width = 0.5) +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position="none") + theme_ipsum() + xlab("") +  ylab("Pseudo F-stat") + labs(fill = "Country") + 
  theme(axis.title.x = element_text(size = 16, family = "Times"), 
        axis.title.y = element_text(size = 16, family = "Times"),
        axis.text.x = element_text(family = "Times", size = 13),
        axis.text.y = element_text(family = "Times", size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)
  ) 
+ theme(legend.title = element_text(family = "Times", size = 14),
            legend.text = element_text(family = "Times", size = 13))


con_cr_pc = tax_glom(con_cr, taxrank = rank_names(con_cr)[6])
con_cr_pc = filter_taxa(con_cr_pc,function(x) mean(x) > 0, TRUE)
con_cr_pc
rank_names(con_cr)[6]

con_cr = filter_taxa(con_cr,function(x) mean(x) > 0, TRUE)

dfotu = as.data.frame(otu_table(con_cr))
dfotu[dfotu > 0] = 1
sc = apply(dfotu,1,sum)*100/ncol(dfotu)
length(sc[sc<=10])

hist(sc,ylim = c(0,500),
     col = 'orange',xlab = "Species prevalence across 2035 samples (%)",
     ylab = "Number of species",breaks=100)

con_cr_share100 = filter_taxa(con_cr,function(x) sum(x!=0) == nsamples(con_cr), TRUE)

con_cr_share90 = filter_taxa(con_cr,function(x) sum(x!=0) >= 0.90*nsamples(con_cr), TRUE)
df_s90 = as.data.frame(tax_table(con_cr_share90))

con_cr_share85 = filter_taxa(con_cr,function(x) sum(x!=0) >= 0.85*nsamples(con_cr), TRUE)
df_s85 = as.data.frame(tax_table(con_cr_share85))

con_cr_share80 = filter_taxa(con_cr,function(x) sum(x!=0) >= 0.8*nsamples(con_cr), TRUE)
df_s85 = as.data.frame(tax_table(con_cr_share80))

dfotu = as.data.frame(otu_table(con_cr))
tax_s90 = sapply(rownames(df_s90), function(x) which(x == rownames(dfotu)))
tax_s90_otu = dfotu[tax_s90,]
r_s90 = sum(tax_s90_otu)/ncol(tax_s90_otu)

con_cr_chn = subset_samples(con_cr, country == "CHN")
con_cr_chn_share90 = filter_taxa(con_cr_chn,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_chn), TRUE)
df_chn_s90 = as.data.frame(tax_table(con_cr_chn_share90))

con_cr_usa = subset_samples(con_cr, country == "USA")
con_cr_usa_share90 = filter_taxa(con_cr_usa,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_usa), TRUE)
df_usa_s90 = as.data.frame(tax_table(con_cr_usa_share90))

con_cr_ind = subset_samples(con_cr, country == "IND")
con_cr_ind_share90 = filter_taxa(con_cr_ind,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_ind), TRUE)
df_ind_s90 = as.data.frame(tax_table(con_cr_ind_share90))

con_cr_mdg = subset_samples(con_cr, country == "MDG")
con_cr_mdg_share90 = filter_taxa(con_cr_mdg,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_mdg), TRUE)
df_mdg_s90 = as.data.frame(tax_table(con_cr_mdg_share90))

con_cr_gbr = subset_samples(con_cr, country == "GBR")
con_cr_gbr_share90 = filter_taxa(con_cr_gbr,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_gbr), TRUE)
df_gbr_s90 = as.data.frame(tax_table(con_cr_gbr_share90))

con_cr_nld = subset_samples(con_cr, country == "NLD")
con_cr_nld_share90 = filter_taxa(con_cr_nld,function(x) sum(x!=0) >= 0.9*nsamples(con_cr_nld), TRUE)
df_nld_s90 = as.data.frame(tax_table(con_cr_nld_share90))

df_con = Reduce(union, list(df_usa_s90$Species,df_ind_s90$Species,df_mdg_s90$Species,
                            df_gbr_s90$Species,df_nld_s90$Species))
df_u_chn = setdiff(df_chn_s90$Species,df_con)

df_con_t = Reduce(union, list(df_chn_s90$Species,df_usa_s90$Species,df_ind_s90$Species,
                              df_mdg_s90$Species,df_gbr_s90$Species,df_nld_s90$Species))

chn_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_chn_s90$Species, 1, 0)))
usa_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_usa_s90$Species, 1, 0)))
ind_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_ind_s90$Species, 1, 0)))
mdg_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_mdg_s90$Species, 1, 0)))
gbr_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_gbr_s90$Species, 1, 0)))
nld_ar = sapply(df_con_t, function(x) return(ifelse(x %in% df_nld_s90$Species, 1, 0)))
core_df = data.frame(chn_ar,usa_ar,ind_ar,mdg_ar,gbr_ar,nld_ar)
colnames(core_df) = c("CHN","USA","IND","MDG","GBR","NLD")
fs3f = pheatmap::pheatmap(core_df,cluster_rows = FALSE,cluster_cols = F, 
                          legend = FALSE, cellwidth = 15, color = c("#A6D854", "#FFD92F"),fontsize = 12)
save(fs3f, file = "fs3f.rda")
load(file = "fs3f.rda")

#==============================================================================
dfotu = as.data.frame(otu_table(con_cr))
dfotu[dfotu > 0] = 1
sc = apply(dfotu,1,sum)*100/ncol(dfotu)

hist(sc,ylim = c(0,500),
     col = 'orange',xlab = "Species prevalence across 3063 samples (%)",
     ylab = "Number of species",breaks=100)
#==============================================================================

m = phyloseq2meco(con_cr)
m$cal_abund()
t1 = trans_abund$new(dataset = m, taxrank = "Phylum", ntaxa = 10, groupmean = "country")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)

df_ab = t1$abund_data 
df_ab_chn = df_ab %>% filter(Sample == "CHN")


getPalette = colorRampPalette(brewer.pal(12, "Paired"))
t1 = trans_abund$new(dataset = m, taxrank = "Genus", ntaxa = 20, groupmean = "country")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE, use_colors = getPalette(20)) 
#======================================================================
#con_crg = tax_glom(con_cr, taxrank = rank_names(con_cr)[2])
#con_cr_chn = subset_samples(con_cr, country == "CHN")

m = phyloseq2meco(con_cr)
m$tidy_dataset()
mc = clone(m)
mc$sample_table = subset(mc$sample_table, country == "IND")
mc$tidy_dataset()
#mc_m = mc$merge_taxa(taxa = "Genus")
dfotu = mc_m$otu_table
dfotu[dfotu > 0] = 1
rownames(dfotu) = sapply(rownames(dfotu),function(x) unlist(strsplit(x, split = "\\|"))[2])
dfc = apply(dfotu,1,function(x) sum(x)/ncol(dfotu))
dfotu$dfc = dfc
dfotu = dfotu[,"dfc", drop = F]
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
