library(curatedMetagenomicData)
library(phyloseq)
library(mia)
library(dplyr)
library(doParallel)
library(caret)
library(UpSetR)
library(iml)
library(pROC)

setwd("e:/ma16/")

con_cr = filter_taxa(con_cr, function(x) mean(x) > 1e-3, TRUE)
#=============================================================
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

dfotu = t(as.data.frame(otu_table(con_cr)))
dfotu = as.data.frame(dfotu)

dfotu$seq_p = sample_data(con_cr)$sequencing_platform
dfotu$ek = sample_data(con_cr)$DNA_extraction_kit

dfotu$country_chn = factor(ifelse(sample_data(con_cr)$country == "CHN", "CHN", "R"))
dfotu = na.omit(dfotu)

rf_id_chn = createDataPartition(dfotu$country_chn, p = 0.8, list = F)
dfotu_chn_train = dfotu[rf_id_chn, ]
dfotu_chn_test = dfotu[-rf_id_chn, ]

set.seed(1)
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = T,
  allowParallel = T,
  summaryFunction=twoClassSummary
)
rftune <- expand.grid(mtry = c(15,20,50,80,100,150,195))
rf_chn_con = caret::train(country_chn~., data = dfotu_chn_train, method = "rf", 
                          metric = "ROC", ntree = 500,tuneGrid = rftune, 
                          trControl = ctrl)
save(rf_chn_con, file = "rfchn.rdata")
save(rf_chn_con, file = "rfchn1.rdata")
load(file = "rfchn.rdata")
chnr = rf_chn_con$results
chnr$country = "CHN"

imptax_chn = varImp(rf_chn_con)$importance
imptax_chn = arrange(imptax_chn, desc(Overall))
rn =  rownames(imptax_chn)
imptax_chn = as.data.frame(imptax_chn[grepl("k_", rn), ])
rownames(imptax_chn) = rn[grepl("k_", rn)]
rownames(imptax_chn) = sapply(rownames(imptax_chn),function(x) unlist(strsplit(x, split = "\\|"))[7])
rownames(imptax_chn) = sapply(rownames(imptax_chn),function(x) unlist(strsplit(x, split = "\\`"))[1])
taxchn = rownames(imptax_chn)[1:30]

#taxchn1 = rownames(imptax_chn)[1:50]
bm_chn = intersect(mr_bm_chn[,1], taxchn)
rownames(mr_bm_chn) = mr_bm_chn[,1]
df_bm = mr_bm_chn[bm_chn,]
colnames(df_bm)[1] = "species"
ggplot(df_bm, aes(x = species, y = -coef)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_ipsum() + xlab("Species") + ylab("Coef") + scale_fill_brewer(palette="BuPu") + coord_flip() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(family = "Times", size = 15),
        legend.text = element_text(family = "Times", size = 14)
  ) 


result.predicted.prob <- predict(rf_chn_con, dfotu_chn_test, type="prob")
result1 <- predict(rf_chn_con, dfotu_chn_test)
cm = confusionMatrix(result1, reference = dfotu_chn_test$country_chn)
roc_chn<- roc(dfotu_chn_test$country_chn, result.predicted.prob[,1])
ggroc(roc_chn, legacy.axes = TRUE, size = 1.2)
#=============================================================
dfotu = t(as.data.frame(otu_table(con_cr)))
dfotu = as.data.frame(dfotu)
dfotu$seq_p = sample_data(con_cr)$sequencing_platform
dfotu$ek = sample_data(con_cr)$DNA_extraction_kit

dfotu$country_usa = factor(ifelse(sample_data(con_cr)$country == "USA", "USA", "R"))
dfotu = na.omit(dfotu)

rf_id_usa = createDataPartition(dfotu$country_usa, p = 0.8, list = F)
save(rf_id_usa, file = "rf_id_usa.rda")
dfotu_usa_train = dfotu[rf_id_usa, ]
dfotu_usa_test = dfotu[-rf_id_usa, ]

rftune <- expand.grid(mtry = c(15,20,50,80,100,150,195))
rf_usa_con = caret::train(country_usa~., data = dfotu_usa_train, method = "rf", 
                          metric = "ROC", ntree = 500,tuneGrid = rftune, trControl = ctrl)
usar = rf_usa_con$results
usar$country = "USA"

save(rf_usa_con, file = "rfusa.rdata")
load(file = "rfusa.rdata")

imptax_usa = varImp(rf_usa_con)$importance
imptax_usa = arrange(imptax_usa, desc(Overall))
rn =  rownames(imptax_usa)
imptax_usa = as.data.frame(imptax_usa[grepl("k_", rn), ])
rownames(imptax_usa) = rn[grepl("k_", rn)]
rownames(imptax_usa) = sapply(rownames(imptax_usa),function(x) unlist(strsplit(x, split = "\\|"))[7])
rownames(imptax_usa) = sapply(rownames(imptax_usa),function(x) unlist(strsplit(x, split = "\\`"))[1])
taxusa = rownames(imptax_usa)[1:30]

bm_usa = intersect(mr_bm_usa[,1], taxusa)
rownames(mr_bm_usa) = mr_bm_usa[,1]
df_bm = mr_bm_usa[bm_usa,]
colnames(df_bm)[1] = "species"
ggplot(df_bm, aes(x = species, y = -coef)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_ipsum() + xlab("Species") + ylab("Coef") + scale_fill_brewer(palette="BuPu") + coord_flip() +
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(family = "Times", size = 15),
        legend.text = element_text(family = "Times", size = 14)
  ) 


result.predicted.prob <- predict(rf_usa_con, dfotu_usa_test, type="prob")
result1 <- predict(rf_usa_con, dfotu_usa_test)
cm = confusionMatrix(result1, reference = dfotu_usa_test$country_usa)
roc_usa <- roc(dfotu_usa_test$country_usa, result.predicted.prob[,1])
ggroc(roc_usa, legacy.axes = TRUE, size = 1.2)
#=============================================================
dfotu = t(as.data.frame(otu_table(con_cr)))
dfotu = as.data.frame(dfotu)
dfotu$seq_p = sample_data(con_cr)$sequencing_platform
dfotu$ek = sample_data(con_cr)$DNA_extraction_kit

dfotu$country_ind = factor(ifelse(sample_data(con_cr)$country == "IND", "IND", "R"))
dfotu = na.omit(dfotu)

rf_id_ind = createDataPartition(dfotu$country_ind, p = 0.8, list = F)
save(rf_id_ind, file = "rf_id_ind.rda")
dfotu_ind_train = dfotu[rf_id_ind, ]
dfotu_ind_test = dfotu[-rf_id_ind, ]

rftune <- expand.grid(mtry = c(15,20,50,80,100,150,195))
rf_ind_con = caret::train(country_ind~., data = dfotu_ind_train, method = "rf", ntree = 1000,tuneGrid = rftune, 
                   metric = "ROC", trControl = ctrl)
indr = rf_ind_con$results
indr$country = "IND"

save(rf_ind_con, file = "rfind.rdata")
load(file = "rfind.rdata")

imptax_ind = varImp(rf_ind_con)$importance
imptax_ind = arrange(imptax_ind, desc(Overall))
rn =  rownames(imptax_ind)
imptax_ind = as.data.frame(imptax_ind[grepl("k_", rn), ])
rownames(imptax_ind) = rn[grepl("k_", rn)]
rownames(imptax_ind) = sapply(rownames(imptax_ind),function(x) unlist(strsplit(x, split = "\\|"))[7])
rownames(imptax_ind) = sapply(rownames(imptax_ind),function(x) unlist(strsplit(x, split = "\\`"))[1])
taxind = rownames(imptax_ind)[1:30]

bm_ind = intersect(mr_bm_ind[,1], taxind)
rownames(mr_bm_ind) = mr_bm_ind[,1]
df_bm = mr_bm_ind[bm_ind,]
colnames(df_bm)[1] = "species"
ggplot(df_bm, aes(x = species, y = -coef)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_ipsum() + xlab("Species") + ylab("Coef") + scale_fill_brewer(palette="BuPu") + coord_flip() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(family = "Times", size = 15),
        legend.text = element_text(family = "Times", size = 14)
  ) 


result.predicted.prob <- predict(rf_ind_con, dfotu_ind_test, type="prob")
result1 <- predict(rf_ind_con, dfotu_ind_test)
cm = confusionMatrix(result1, reference = dfotu_ind_test$country_ind)
roc_ind <- roc(dfotu_ind_test$country_ind, result.predicted.prob[,1])
ggroc(roc_ind, legacy.axes = TRUE, size = 1.2)
#=============================================================
dfotu = t(as.data.frame(otu_table(con_cr)))
dfotu = as.data.frame(dfotu)
dfotu$seq_p = sample_data(con_cr)$sequencing_platform
dfotu$ek = sample_data(con_cr)$DNA_extraction_kit
dfotu$country_mdg = factor(ifelse(sample_data(con_cr)$country == "MDG", "MDG", "R"))
dfotu = na.omit(dfotu)

rf_id_mdg = createDataPartition(dfotu$country_mdg, p = 0.8, list = F)
save(rf_id_mdg, file = "rf_id_mdg.rda")
dfotu_mdg_train = dfotu[rf_id_mdg, ]
dfotu_mdg_test = dfotu[-rf_id_mdg, ]

rftune <- expand.grid(mtry = c(15,20,50,80,100,150,195))
rf_mdg_con = caret::train(country_mdg~., data = dfotu_mdg_train, method = "rf", ntree = 500,tuneGrid = rftune, 
                   metric = "ROC", trControl = ctrl)
mdgr = rf_mdg_con$results
mdgr$country = "MDG"

save(rf_mdg_con, file = "rfmdg.rdata")
load(file = "rfmdg.rdata")

imptax_mdg = varImp(rf_mdg_con)$importance
imptax_mdg = arrange(imptax_mdg, desc(Overall))
rn =  rownames(imptax_mdg)
imptax_mdg = as.data.frame(imptax_mdg[grepl("k_", rn), ])
rownames(imptax_mdg) = rn[grepl("k_", rn)]
rownames(imptax_mdg) = sapply(rownames(imptax_mdg),function(x) unlist(strsplit(x, split = "\\|"))[7])
rownames(imptax_mdg) = sapply(rownames(imptax_mdg),function(x) unlist(strsplit(x, split = "\\`"))[1])
taxmdg = rownames(imptax_mdg)[1:30]

bm_mdg = intersect(mr_bm_mdg[,1], taxmdg)
rownames(mr_bm_mdg) = mr_bm_mdg[,1]
df_bm = mr_bm_mdg[bm_mdg,]
colnames(df_bm)[1] = "species"
ggplot(df_bm, aes(x = species, y = -coef)) + geom_bar(stat = "identity", position = "dodge", width = 0.1) + 
  theme_ipsum() + xlab("Species") + ylab("Coef") + scale_fill_brewer(palette="BuPu") + coord_flip() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(family = "Times", size = 15),
        legend.text = element_text(family = "Times", size = 14)
  ) 

result.predicted.prob <- predict(rf_mdg_con, dfotu_mdg_test, type="prob")
result1 <- predict(rf_mdg_con, dfotu_mdg_test)
cm = confusionMatrix(result1, reference = dfotu_mdg_test$country_mdg)
roc_mdg <- roc(dfotu_mdg_test$country_mdg, result.predicted.prob[,1])
ggroc(roc_mdg, legacy.axes = TRUE, size = 1.2)
#=============================================================
dfotu = t(as.data.frame(otu_table(con_cr)))
dfotu = as.data.frame(dfotu)
dfotu$seq_p = sample_data(con_cr)$sequencing_platform
dfotu$ek = sample_data(con_cr)$DNA_extraction_kit
dfotu$country_gbr = factor(ifelse(sample_data(con_cr)$country == "GBR", "GBR", "R"))
dfotu = na.omit(dfotu)

rf_id_gbr = createDataPartition(dfotu$country_gbr, p = 0.8, list = F)
save(rf_id_gbr, file = "rf_id_gbr.rda")
dfotu_gbr_train = dfotu[rf_id_gbr, ]
dfotu_gbr_test = dfotu[-rf_id_gbr, ]

rftune <- expand.grid(mtry = c(15,20,50,80,100,150,195))
rf_gbr_con = caret::train(country_gbr~., data = dfotu_gbr_train, method = "rf", ntree = 500,tuneGrid = rftune, 
                          metric = "ROC", trControl = ctrl)
gbrr = rf_gbr_con$results
gbrr$country = "GBR"

save(rf_gbr_con, file = "rfgbr.rdata")
load(file = "rfgbr.rdata")

imptax_gbr = varImp(rf_gbr_con)$importance
imptax_gbr = arrange(imptax_gbr, desc(Overall))
rn =  rownames(imptax_gbr)
imptax_gbr = as.data.frame(imptax_gbr[grepl("k_", rn), ])
rownames(imptax_gbr) = rn[grepl("k_", rn)]
rownames(imptax_gbr) = sapply(rownames(imptax_gbr),function(x) unlist(strsplit(x, split = "\\|"))[7])
rownames(imptax_gbr) = sapply(rownames(imptax_gbr),function(x) unlist(strsplit(x, split = "\\`"))[1])
taxgbr = rownames(imptax_gbr)[1:30]

bm_gbr = intersect(mr_bm_gbr[,1], taxgbr)
rownames(mr_bm_gbr) = mr_bm_gbr[,1]
df_bm = mr_bm_gbr[bm_gbr,]
colnames(df_bm)[1] = "species"
f5b_gbr = ggplot(df_bm, aes(x = species, y = -coef)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_ipsum() + xlab("Species") + ylab("Coef") + scale_fill_brewer(palette="BuPu") + coord_flip() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(family = "Times", size = 15),
        legend.text = element_text(family = "Times", size = 14)
  ) 
save(f5b_gbr, file = "f5b1_gbr.rda")

result.predicted.prob <- predict(rf_gbr_con, dfotu_gbr_test, type="prob")
result1 <- predict(rf_gbr_con, dfotu_gbr_test)
cm = confusionMatrix(result1, reference = dfotu_gbr_test$country_gbr)
roc_gbr <- roc(dfotu_gbr_test$country_gbr, result.predicted.prob[,1])
ggroc(roc_gbr, legacy.axes = TRUE, size = 1.2)
#=============================================================
dfotu = t(as.data.frame(otu_table(con_cr)))
dfotu = as.data.frame(dfotu)
dfotu$seq_p = sample_data(con_cr)$sequencing_platform
dfotu$ek = sample_data(con_cr)$DNA_extraction_kit
dfotu$country_nld = factor(ifelse(sample_data(con_cr)$country == "NLD", "NLD", "R"))
dfotu = na.omit(dfotu)

rf_id_nld = createDataPartition(dfotu$country_nld, p = 0.8, list = F)
save(rf_id_nld, file = "rf_id_nld.rda")
dfotu_nld_train = dfotu[rf_id_nld, ]
dfotu_nld_test = dfotu[-rf_id_nld, ]

rftune <- expand.grid(mtry = c(15,20,50,80,100,150,195))
rf_nld_con = caret::train(country_nld~., data = dfotu_nld_train, method = "rf", ntree = 500,tuneGrid = rftune, 
                   metric = "ROC", trControl = ctrl)
nldr = rf_nld_con$results
nldr$country = "NLD"

save(rf_nld_con, file = "rfnld.rdata")
load(file = "rfnld.rdata")

imptax_nld = varImp(rf_nld_con)$importance
imptax_nld = arrange(imptax_nld, desc(Overall))
rn =  rownames(imptax_nld)
imptax_nld = as.data.frame(imptax_nld[grepl("k_", rn), ])
rownames(imptax_nld) = rn[grepl("k_", rn)]
rownames(imptax_nld) = sapply(rownames(imptax_nld),function(x) unlist(strsplit(x, split = "\\|"))[7])
rownames(imptax_nld) = sapply(rownames(imptax_nld),function(x) unlist(strsplit(x, split = "\\`"))[1])
taxnld = rownames(imptax_nld)[1:30]

bm_nld = intersect(mr_bm_nld[,1], taxnld)
rownames(mr_bm_nld) = mr_bm_nld[,1]
df_bm = mr_bm_nld[bm_nld,]
colnames(df_bm)[1] = "species"
f5b_nld = ggplot(df_bm, aes(x = species, y = -coef)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_ipsum() + xlab("Species") + ylab("Coef") + scale_fill_brewer(palette="BuPu") + coord_flip() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(family = "Times", size = 15),
        legend.text = element_text(family = "Times", size = 14)
  ) 
save(f5b_nld, file = "f5b1_nld.rda")  

result.predicted.prob <- predict(rf_nld_con, dfotu_nld_test, type="prob")
result1 <- predict(rf_nld_con, dfotu_nld_test)
cm = confusionMatrix(result1, reference = dfotu_nld_test$country_nld)
roc_nld <- roc(dfotu_nld_test$country_nld, result.predicted.prob[,1])
ggroc(roc_nld, legacy.axes = TRUE, size = 1.2)
ggroc(list(roc_chn, roc_usa, roc_ind, roc_mdg, roc_gbr, roc_nld), 
      legacy.axes = TRUE, size = 1.2, alpha = 0.5) + theme_minimal() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)
  ) + guides(color = guide_legend(title = NULL)) + 
  scale_color_discrete(
    labels = c("CHN", "USA", "IND", "MDG", "GBR", "NLD")
  )

   
listtax = list(CHN = taxchn, USA = taxusa, IND = taxind, MDG = taxmdg, GBR = taxgbr, 
               NLD = taxnld) 
upset(fromList(listtax), nsets = 6, text.scale = 3)


cr = rbind(chnr, usar, indr, mdgr, gbrr, nldr)
cr$mtry = factor(cr$mtry)
cr$country = factor(cr$country)
ggplot(cr, aes(x = mtry, y = ROC, color = country, group = country)) + geom_point(size = 3) + geom_line()  + theme_ipsum() + 
  theme(axis.title.x = element_text(size = 18, family = "Times", face = "bold"), 
        axis.title.y = element_text(size = 18, family = "Times", face = "bold"),
        axis.text.x = element_text(family = "Times", size = 15),
        axis.text.y = element_text(family = "Times", size = 15),
        legend.title = element_text(family = "Times", size = 15),
        legend.text = element_text(family = "Times", size = 14)
  ) + labs(color = "Country") 

stopCluster(cl)
