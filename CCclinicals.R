####
#
# -- short report and meta analysis to the journal "Science Progress" -----------
#
#

setwd("/Users/seidmuhie/Documents/06032023_CNN_Application_MedicalImage_Classification/10072023_Cervical Cancer for Science Progress/10102023_Science Progress")

# read the cervical cancer risk factors (table from caracas venzuala)

cc_risk_factors_caracas = read.csv("risk_factors_cervical_cancer_Caracas_Venezuela.csv")
View(cc_risk_factors_caracas)
str(cc_risk_factors_caracas)

risk_factors_uganda = read.csv("41588_2020_673_MOESM3_ESM_uganda_ng.csv")
str(risk_factors_uganda)
View(risk_factors_uganda)
names(risk_factors_uganda)

library(xtable)
test_table = xtabs(Stage ~ HPV.clade + HPV.type + Grade + HIV.status + Final.histology, data = risk_factors_uganda)
ftable(test_table)


test_table = xtabs(~ Stage + HPV.clade  +  HIV.status + Final.histology , data = risk_factors_uganda)
ftable(test_table)
summary(test_table)


risk_count <- xtabs(Stage ~ HPV.clade + HPV.type + HIV.status + Final.histology, data = risk_factors_uganda)

test_table2 = table(risk_factors_uganda[,c(5,6,8,10,14)])
ftable(test_table2)
test_bable3 = as.data.frame(test_table2)
View(test_bable3)
write.csv(test_bable3, "testfrequency table for uganda data.csv")
ftable(test_bable3[,c(1,2,3,6)])

test_table4 = aggregate(Freq ~ Stage + HPV.clade, 
          data=test_bable3, 
          FUN=sum)

View(test_table4)


install.packages('epiDisplay')
library(epiDisplay)
#tab1(mtcars$cyl, sort.group = "decreasing", cum.percent = TRUE)
tab1(risk_factors_uganda$Stage, sort.group = "decreasing", cum.percent = TRUE)



# CrossTable in R
install.packages('gmodels')
library(gmodels)
CrossTable(risk_factors_uganda$Stage,  risk_factors_uganda$HIV.status, prop.t=TRUE, prop.r=TRUE, prop.c=TRUE)


# reading TCGA data sets and merging them
tcga_dataset1 = read.csv("TCGA Supplemental Table 1-CESC_Data_Summary_Table_Supplement_V2.csv")
tcga_dataset2 = read.csv("TCGA_Supplemental Table 2-CESC_Extended_Set_Data_Summary_Table_V2.csv")
names(tcga_dataset1)
tcga_dataset1

tab1(tcga_dataset1$CLIN.clinStage, sort.group = "decreasing", cum.percent = TRUE)

tab1(tcga_dataset2$CLIN.clinStage, sort.group = "decreasing", cum.percent = TRUE)

tab1(cc_risk_factors_caracas$Number.of.sexual.partners, sort.group = "decreasing", cum.percent = TRUE)

tab1(risk_factors_uganda$Stage_Pure, sort.group = "decreasing", cum.percent = TRUE)
tab1(risk_factors_uganda$Stage, sort.group = "decreasing", cum.percent = TRUE)


# read the cell genomics data from china
cell_genomics_china = read.csv("cell geneomics Supplemental_Tables_china.csv")
tab1(cell_genomics_china$Stage_2009_adj, sort.group = "decreasing", cum.percent = TRUE)


# read the dataset from Uganda nature genetics
risk_factors_uganda = read.csv("41588_2020_673_MOESM3_ESM_uganda_ng.csv")
tab1(risk_factors_uganda$Stage2, sort.group = "decreasing", cum.percent = TRUE)


### correlation matrix
names(cc_risk_factors_caracas)
names(risk_factors_uganda)
names(cell_genomics_china)
names(tcga_dataset1)




# read the cervical cancer risk factors (table from caracas venzuala)

cc_risk_factors_caracas = read.csv("risk_factors_cervical_cancer_Caracas_Venezuela.csv")
View(cc_risk_factors_caracas)
str(cc_risk_factors_caracas)
names(cc_risk_factors_caracas)


res <- cor(cc_risk_factors_caracas)

install.packages("Hmisc")
library("Hmisc")
res2 <- rcorr(as.matrix(cc_risk_factors_caracas))
res2


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res3 = flattenCorrMatrix(res2$r, res2$P)

View(res3)

install.packages("corrplot")
library(corrplot)
corrplot(res3, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

write.csv(res3, "uganda data correlation between the different variables.csv", row.names = F)



#################################################################################################################
#
#
#-------------- two groups of combined files for correlation and stage related factors identification ---------
#
#

#
#-----the variable importances from Regularized Random Forest (RRF) algorithm--------
#
# Train an RRF model and compute variable importance.
setwd("/Users/seidmuhie/Documents/06032023_CNN_Application_MedicalImage_Classification/10072023_Cervical Cancer for Science Progress/10102023_Science Progress")
rm (list = ls())
# read the combination of the cervical cancer data data sets from the three cohorts
tcga_mocc_htmcp = read.csv("For correlation TCGA HTMCP-Uganda MOCC-China.csv")
View(tcga_mocc_htmcp)

library("Hmisc")


tcga_mocc_htmcp.cor = cor(tcga_mocc_htmcp)



# read the dataset from Caracas (cervical cancer data set)
cc_caracas = read.csv("For correlation caracase Venezuela.csv")
names(cc_caracas)
str(cc_caracas)

par(mar = c(4, 4, 2, 2)) # Set the margin on all sides to 2
cc_caracas.cor = cor(cc_caracas[,-c(14,18)], use='pairwise.complete.obs')
View(cc_caracas.cor)
#install.packages("corrplot")
library(corrplot)

corrplot(cc_caracas.cor, type="upper", order="hclust")

library(RColorBrewer)
corrplot(cc_caracas.cor, method = "number", type="upper", order="hclust", 
         col=brewer.pal(n=8, name="RdBu"))


# read the one hot encoded cervical cancer data
tcga_mocc_htmcp_ohe = read.csv("For correlation TCGA HTMCP-Uganda MOCC-China_onehotendoded.csv")
View(tcga_mocc_htmcp_ohe)
str(tcga_mocc_htmcp_ohe)
names(tcga_mocc_htmcp_ohe)

tcga_mocc_htmcp_ohe.cor = cor(tcga_mocc_htmcp_ohe[,-c(1:2,31,33:34)], use='pairwise.complete.obs')
View(tcga_mocc_htmcp_ohe.cor)
#install.packages("corrplot")
library(corrplot)

corrplot(tcga_mocc_htmcp_ohe.cor)


tcga_mocc_htmcp_ohe.cor = cor(tcga_mocc_htmcp_ohe[,-c(1:2,31,33:34)], use='complete.obs')
View(tcga_mocc_htmcp_ohe.cor)
#install.packages("corrplot")
library(corrplot)

corrplot(tcga_mocc_htmcp_ohe.cor)

library(epiDisplay)
# read the combination of the cervical cancer data data sets from the three cohorts
tcga_mocc_htmcp = read.csv("For correlation TCGA HTMCP-Uganda MOCC-China.csv")
View(tcga_mocc_htmcp)
names(tcga_mocc_htmcp)
tab1(tcga_mocc_htmcp$Stage_Pure, sort.group = "decreasing", cum.percent = TRUE)

df = tcga_mocc_htmcp[,c(3,15)]



tableStack(as.data.frame(tcga_mocc_htmcp, tcga_mocc_htmcp$Stage))
tabpct(tcga_mocc_htmcp$Cohort, tcga_mocc_htmcp$Stage_Pure)

names(df)
myfrequency = df %>%
  group_by(Stage_Pure, Region) %>%
  count() # %>% tab1()

View(myfrequency)
names(myfrequency)

myfrequency2 = myfrequency[1:10,]
View(myfrequency2)

library(ggplot2)
library(viridis)
library(hrbrthemes)

# Graph
ggplot(myfrequency2, aes(fill=Stage_Pure, x=n, y=Stage_Pure)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = T, option = "E") +
  ggtitle("frequency of cervical cancer stages from three regions") +
  facet_wrap(~ Region) +
  xlab("number of cervical cancer patients") +
  ylab("cervial cancer stages")+
 
  theme_ipsum(plot_title_margin = 15, plot_title_size = 16, plot_title_face = "plain",axis_title_size = 13,axis_title_just = "middle") +
  #theme(panel.background = element_blank(), axis.line = element_line(colour = "grey"))+
  theme(legend.position="none") +
  geom_text(aes(label = sprintf("%.0f", n), x= n),  hjust = -0.1)+
  xlab("number of cervical cancer patients") +
  ylab("cervial cancer stages")+
  
  #theme(#panel.grid.minor = element_blank(), #panel.grid.major = element_blank(),
       # panel.background = element_blank(), axis.line = element_line(colour = "grey"))+
  #scale_x_discrete(expand = expansion(add = .6)) +
  scale_x_continuous(limits = c(0, 150), expand = c(0.1, 10)) 
  #theme(axis.text=element_text(size=12), #change font size of axis text
        #axis.title=element_text(size=13), #change font size of axis titles
        #plot.title=element_text(size=15), #change font size of plot title
        #legend.text=element_text(size=20), #change font size of legend text
        #legend.title=element_text(size=20)) #change font size of legend title  

  #ylab("cervial cancer stages")




################################################################
install.packages("RRF") # Regularized Random Forest
library(RRF)

set.seed(100)

cc_dx = as.factor(cc_caracas$Dx.Cancer)
cc_caracas_selected = cc_caracas[,-c(2,14,16:18)]

# first imputing NA values
cc_caracas_impute = rrfImpute(cc_caracas_selected, cc_dx)
rrfMod <- RRF(cc_dx ~ ., cc_caracas_impute, method="RRF")

varImpPlot(rrfMod, sort=TRUE, n.var=min(30, nrow(rrfMod$importance)),
           type=NULL, class=NULL, scale=TRUE, 
           main="Cervical cancer predictors\n ranked using regularized random forest") 



rrfImp <- varImp(rrfMod, 39, scale=F)

rrfImp[["importance"]] = importance; rrfImp
importance = read.csv("01142020_SelectedClassifiers/DiscoveryMale_RRFvarImportance_modified.csv", row.names = 1)


plot(rrfImp, top = 50, main='DiscoveryMale RRF Variable Importance')
rrf_importance = rrfImp$importance
View(rrf_importance)

plot(rrfImp, top = 20, main='DiscoveryMale RRF Variable Importance')

str(rrfImp)

nrow(varImp(rrfMod)$importance) #1305 variables extracted




tcga_mocc_htmcp_ohe = read.csv("For correlation TCGA HTMCP-Uganda MOCC-China_onehotendoded_v2.csv", row.names = 1)
View(tcga_mocc_htmcp_ohe)
str(tcga_mocc_htmcp_ohe)
names(tcga_mocc_htmcp_ohe)
is.na(tcga_mocc_htmcp_ohe)



# for ranking of risk factors using regularized random forest
#install.packages("RRF")
library(RRF)

set.seed(100)

cc_stage1 = as.factor(tcga_mocc_htmcp_ohe$Stage_I)

cc_tcga_mocc_htmcp_ohe = tcga_mocc_htmcp_ohe[,-1]
View(cc_tcga_mocc_htmcp_ohe)
# first imputing NA values
rrfMod_stage <- RRF(cc_tcga_mocc_htmcp_ohe, cc_stage1)
View(rrfMod_stage$importance)
varImpPlot(rrfMod_stage, sort=TRUE, n.var=min(30, nrow(rrfMod_stage$importance)),
           type=NULL, class=NULL, scale=TRUE, 
           main="Cervical cancer stage 1 predictors\n ranked using regularized random forest") 


names(tcga_mocc_htmcp_ohe)
tcga_mocc_htmcp_ohe_short = tcga_mocc_htmcp_ohe[,-c(10:12,25:27)]
tcga_mocc_htmcp_ohe.cor = cor(tcga_mocc_htmcp_ohe_short, use='pairwise.complete.obs')
View(tcga_mocc_htmcp_ohe.cor)
#install.packages("corrplot")
library(corrplot)
tcga_mocc_htmcp_ohe.cor_ordered = corrMatOrder(tcga_mocc_htmcp_ohe.cor)

corrplot(tcga_mocc_htmcp_ohe.cor, type="upper", order="hclust")

library(RColorBrewer)
corrplot(tcga_mocc_htmcp_ohe.cor, type="upper", order="hclust", 
         col=brewer.pal(n=8, name="RdBu"))




library(caret)
library(tidyverse)
varImp(rrfMod)$importance %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(Overall) %>%
  mutate(rowname = forcats::fct_inorder(rowname )) %>%
  ggplot()+
  geom_col(aes(x = rowname, y = Overall))+
  coord_flip()+
  theme_bw()

DiscoveryMale_RRFvarImportance = varImp(rrfMod)$importance %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(Overall) %>%
  mutate(rowname = forcats::fct_inorder(rowname ))
View(DiscoveryMale_RRFvarImportance)
write.csv(DiscoveryMale_RRFvarImportance, "01142020_SelectedClassifiers/DiscoveryMale_RRFvarImportance.csv")


#================================================================================
#
#------Least Absolute Shrinkage and Selection Operator (LASSO) Regression--------
#
#      is a type of regularization method that penalizes with L1-norm.
#

library(glmnet)


x <- as.matrix(biomarkermaledata[,-1]) # all X vars
y <- as.double(as.matrix(ifelse(biomarkermaledata[, 1]=='Negative', 0, 1))) # Only Class

# Fit the LASSO model (Lasso: Alpha = 1)
set.seed(100)
cv.lasso <- cv.glmnet(x, y, family='binomial', alpha=1, standardize=TRUE, type.measure='auc')

# Results
plot(cv.lasso)

# plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
cat('Min Lambda: ', cv.lasso$lambda.min, '\n 1Sd Lambda: ', cv.lasso$lambda.1se)
df_coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)

# See all contributing variables
df_coef[df_coef[, 1] != 0, ]
View(df_coef)
write.csv(df_coef, "01142020_SelectedClassifiers/DiscoveryMale_CV.LASSO_selected.csv")


#######################################################################################
#
#------------- statistical evaluation of the three cohorts--------------------
#

# read the file combining the cohorts from Uganda, USA and China
cc_clinical_3cohorts = read.csv("For correlation TCGA HTMCP-Uganda MOCC-China.csv", row.names = 1)
names(cc_clinical_3cohorts)
length(cc_clinical_3cohorts$Stage)

library(purrr)
#cc_clinical_3cohorts %>% map(table)

cc_stat = summary(cc_clinical_3cohorts$Age)

library(plyr)
cc_data_freq = cc_clinical_3cohorts[,-c(1,11)] %>% map(count)

n.obs <- sapply(cc_data_freq, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(cc_data_freq, "[", i = seq.max))

View(mat)


write.csv(mat, "cc_frequency_table.csv")

# subset for each stage (I, II and III_IV)
cc_clinical_3cohorts = read.csv("For correlation TCGA HTMCP-Uganda MOCC-China.csv", row.names = 1)
names(cc_clinical_3cohorts)
length(cc_clinical_3cohorts$Stage)


cc_clinical_3cohorts_stage1 = cc_clinical_3cohorts

library(plyr)

# subset and count features corresponding for the stage I patients
cc_data_freq_stage1 = cc_clinical_3cohorts[cc_clinical_3cohorts$Stage_III_IV=='I',-c(1,2,12)] %>% map(count)

n.obs_stage1 <- sapply(cc_data_freq_stage1, length)
seq.max_stage1 <- seq_len(max(n.obs_stage1))
mat_stage1 <- t(sapply(cc_data_freq_stage1, "[", i = seq.max_stage1))

View(mat_stage1)


# subset and count features corresponding for the stage II patients
cc_data_freq_stage2 = cc_clinical_3cohorts[cc_clinical_3cohorts$Stage_III_IV=='II',-c(1,2,12)] %>% map(count)

n.obs_stage2 <- sapply(cc_data_freq_stage2, length)
seq.max_stage2 <- seq_len(max(n.obs_stage2))
mat_stage2 <- t(sapply(cc_data_freq_stage2, "[", i = seq.max_stage2))

View(mat_stage2)

# subset and count features corresponding for the stage III_IV patients
cc_data_freq_stage3_4 = cc_clinical_3cohorts[cc_clinical_3cohorts$Stage_III_IV=='III_IV',-c(1,2,12)] %>% map(count)

n.obs_stage3_4 <- sapply(cc_data_freq_stage3_4, length)
seq.max_stage3_4 <- seq_len(max(n.obs_stage3_4))
mat_stage3_4 <- t(sapply(cc_data_freq_stage3_4, "[", i = seq.max_stage3_4))

View(mat_stage3_4)

# row bind the three outputs
cc_stages_freq = rbind(mat_stage1, mat_stage2, mat_stage3_4)
View(cc_stages_freq)
write.csv(cc_stages_freq, "Table 3 Frequencies of features related to each of the cervical cancer stages.csv")


# col bind the three outputs
cc_stages_freq_cbind = cbind(mat_stage1, mat_stage2, mat_stage3_4)
View(cc_stages_freq_cbind)
cc_stages_freq_cbind <- apply(cc_stages_freq_cbind,2,as.character)
write.csv(cc_stages_freq_cbind, "Table 3 Frequencies of features related to each of the cervical cancer stages_cbind.csv")

plot(cc_stages_freq_cbind)



###############################
#
#------------ plot the frequencies of features of cc ----
#

cc_freq = read.csv("Table 3 Frequencies of features of cc.csv")
View(cc_freq)
names(cc_freq)

cc_freq_I = cc_freq[cc_freq$Stage=="I",]


library(ggplot2)
library(viridis)
library(hrbrthemes)

# Graph
ggplot(cc_freq, aes(fill=Feature_Subtype, x=Freq, y=reorder(Feature_Subtype, reorder(Feature,Freq)))) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = T, option = "E") +
  #ggtitle("frequency of cervical cancer stages from three regions") +
  facet_wrap( ~ Stage ) +
  xlab("number of cervical cancer patients") +
  ylab("cervial cancer stages")+
  
  theme_ipsum(plot_title_margin = 15, plot_title_size = 16, plot_title_face = "plain",axis_title_size = 13,axis_title_just = "middle") +
  #theme(panel.background = element_blank(), axis.line = element_line(colour = "grey"))+
  theme(legend.position="none") +
  geom_text(aes(label = sprintf("%.0f", Freq), x= Freq),  hjust = -0.1)+
  xlab("number of cervical cancer patients") +
  ylab("clinical and demographic features")+
  
  #theme(#panel.grid.minor = element_blank(), #panel.grid.major = element_blank(),
  # panel.background = element_blank(), axis.line = element_line(colour = "grey"))+
  #scale_x_discrete(expand = expansion(add = .6)) +
  scale_x_continuous(limits = c(0, 230), expand = c(0.1, 10)) 
#theme(axis.text=element_text(size=12), #change font size of axis text
#axis.title=element_text(size=13), #change font size of axis titles
#plot.title=element_text(size=15), #change font size of plot title
#legend.text=element_text(size=20), #change font size of legend text
#legend.title=element_text(size=20)) #change font size of legend title  

#ylab("cervial cancer stages")





