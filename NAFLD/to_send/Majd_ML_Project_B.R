###Machine Learning Project
###Seyyed Amirreza Mousavi Majd
###610799001
###Question B:

###Initialization
###we employ the same codes for importing the data.
set.seed(142241)
library(ggplot2)
library(randomForest)
library(dplyr)
library(caret)
library(e1071)
meta_data=read.csv("meta_data.csv")
meta_data=na.omit(meta_data)
rownames(meta_data) = meta_data[,1]
counts_data=read.csv("Normal.counts.voom.txt",header = T,sep="\t")
counts_data=na.omit(counts_data)
counts_data= counts_data %>% rename("DLDR_0117"="DLDR_1017")
genes = counts_data[,1]
a = as.data.frame(t(counts_data[,-1]))
colnames(a) = genes
counts_data = a
rm(a,genes)
counts_data$Patient_ID = rownames(counts_data)
feature_vec = na.omit(right_join(counts_data,meta_data,by="Patient_ID"))

feature_vec = feature_vec %>% mutate(Run=NULL,SEX=as.factor(SEX),
                                     Diabet = as.factor(Diabet),
                                     Simplified_class = as.factor(Simplified_class))

View(feature_vec)

###from question A we had:
SD_rf_rfe = c("ENSG00000173917","ENSG00000164114","ENSG00000188517",
"ENSG00000177951","ENSG00000112139","ENSG00000078549","ENSG00000163017",
"ENSG00000164588","ENSG00000185585","ENSG00000006704")

rm(meta_data,counts_data)


###let's build a model that uses SD_rf_rfe as well as BMI, SEX, Age, 
###and Diabet as input

###svm radial and random forest were the best methods with SD_rf_rfe




rf_SD_rf_rfe = train( Simplified_class~., data = feature_vec[,c(SD_rf_rfe,"Simplified_class")]
                      , method = "rf",
                      metric = "Accuracy",
                      trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                      preProcess = c("center","scale"))

svm_radial_SD_rf_rfe =  train( Simplified_class~., data = feature_vec[,c(SD_rf_rfe,"Simplified_class")]
                               , method = "svmRadial",
                               metric = "Accuracy",
                               trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                               preProcess = c("center","scale"))





svm_radial_SD_rf_rfe_extended =  train( Simplified_class~., 
                               data = feature_vec[,c(SD_rf_rfe,
                                                     "SEX","BMI_surg",
                                                     "Age","Diabet",
                                                     "Simplified_class")],
                               method = "svmRadial",
                               metric = "Accuracy",
                               trControl = trainControl(method="repeatedcv", 
                                                        number=10, repeats=3), 
                               preProcess = c("center","scale"))

rf_SD_rf_rfe_extended = train( Simplified_class~., 
                      data = feature_vec[,c(SD_rf_rfe,"SEX","BMI_surg",
                                            "Age","Diabet","Simplified_class")],
                      method = "rf",
                      metric = "Accuracy",
                      trControl = trainControl(method="repeatedcv",
                                               number=10, repeats=3), 
                      preProcess = c("center","scale"))

###Perform random forest with recursive feature elimination; Extended version
rf_rfe_extended <- rfe(x = feature_vec %>% select(all_of(SD_rf_rfe),SEX,
                                         BMI_surg,Age,Diabet),
              y = (feature_vec %>% select(Simplified_class))[,1], 
              size=14, 
              rfeControl=rfeControl(functions=rfFuncs, method="cv", number=10))

###We will observe that the extended machines improve a little; But this is
###because they have 14 features in contrast to 10. Thus their better performance
###is not surprising.
print(svm_radial_SD_rf_rfe)
print(svm_radial_SD_rf_rfe_extended)
print(rf_SD_rf_rfe)
print(rf_SD_rf_rfe_extended)




print(rf_rfe_extended)


print(predictors(rf_rfe_extended))
###Notice that SEX,Age,BMI_surge,Diabet always reside in the last features 
###to be selected, thus they are useless. THe transcription dataset is 
###sufficient.
