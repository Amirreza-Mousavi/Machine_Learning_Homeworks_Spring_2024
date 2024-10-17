###Machine Learning Project
###S.Amirreza Mousavi Majd
###St.No. : 610799001
###Currently, I do not possess the required systems to perform large scale
###sequence alignment for the requested tasks(linux virtual machine, CPU, etc),
###nor the internet speed rate is sufficient for the download of the 
###sequencing data.(I'm not residing in Tehran)
###However, I have previously done a somewhat similar pipeline for the analysis
###of transcriptomics data of S. cerevisiae.
###if you are interested in the linux implementations of the pipeline, Please Refer to 
###https://github.com/Amirreza-Mousavi/Rodriguez2017_RNA_Seq_Analysis_Replication/blob/main/script.sh


###Pipeline: sratoolkit-fastqc-Trimmomatic-hisat2-featureCounts-DEseq2


###Link to the complete project:
###https://github.com/Amirreza-Mousavi/Rodriguez2017_RNA_Seq_Analysis_Replication/tree/main


###Here we continue from the meta_data.csv and Normal.counts.voom.txt data

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################


###initialization
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

###someone has messed with the Normal.counts.voom.txt data, because DLDR0117 
###does not exist, and DLDR1017 exists, which is probably a typo. 
###(The range should be from DLDR0001 to DLDR0192)
###we should change its name to the correct one, manually.

counts_data= counts_data %>% rename("DLDR_0117"="DLDR_1017")
###A somewhat dirty transposition
genes = counts_data[,1]
a = as.data.frame(t(counts_data[,-1]))
colnames(a) = genes
counts_data = a
rm(a,genes)


###data wrangling
counts_data$Patient_ID = rownames(counts_data)

View(counts_data)
View(meta_data)

feature_vec = na.omit(right_join(counts_data,meta_data,by="Patient_ID"))

feature_vec = feature_vec %>% mutate(Run=NULL,SEX=as.factor(SEX),
                       Diabet = as.factor(Diabet),
                       Simplified_class = as.factor(Simplified_class))

View(feature_vec)

#picks from counts data: Unknown (to be determined)
#Let's use some random feature vector of size = 8 and see what happens
#SD_Rand = select data A

SD_Rand = sample(colnames(counts_data[-17397]),size = 10)

N = nrow(feature_vec) # = 192


###The selection of features must be improved


###Feature Selection based on Principle Component Analysis

pca_obj=prcomp(counts_data[,-17397]) #17397 = Patient_ID = not numeric

ggplot()+geom_point(aes(x = pca_obj$x[,1] , y=pca_obj$x[,2],colour = (feature_vec$Simplified_class)))
ggplot()+geom_point(aes(x = pca_obj$x[,1] , y=pca_obj$x[,3],colour = (feature_vec$Simplified_class)))
ggplot()+geom_point(aes(x = pca_obj$x[,2] , y=pca_obj$x[,3],colour = (feature_vec$Simplified_class)))

###In a somewhat naive way for feature selection we choose the ones with high stdev.

print(pca_obj$sdev)
###we pick the first 10 PCAs all of them have stdev>10
pca_score = pca_obj$rotation[,1:10]
pca_score = abs(pca_score) #absolute value is imp

l = list()
for (j in 1:10){
l =  c(l,  names((sort(pca_score[,j],decreasing = T))[1:10])) #pick top ten from each pca
}

SD_PCA=unique(x = unlist(l)) #37 genes were duplicates in top ones



############
####Feature Selection by Filtering out highly correlated gene expressions
cor_table=cor(counts_data[,-17397]) #takes about a minute
fc=findCorrelation(cor_table,cutoff = 0.7,names = F) #to_be removed
counts_data_filtered = counts_data[,-fc]
cat(paste("In the first step, we reduced our",dim(counts_data)[2],"genes to"
          ,dim(counts_data_filtered)[2],"genes.\n"))



####remove useless intermediate data for better run
rm(cor_table,l,pca_obj,pca_score,fc,j,N)


####Implementing Recursive Feature Elimination(RFE) on the
####filtered datasets using a generic random forest model.

control <- rfeControl(functions=rfFuncs, method="cv", number=10)

####The code below takes about 3 minutes to be completed
rf_rfe <- rfe(x = counts_data_filtered[,-7753],
               y = feature_vec[,"Simplified_class"] , 
               sizes=c(1:15), rfeControl=control) #7753 = patient id

print(rf_rfe) # summarize the results

predictors(rf_rfe) # list the chosen features

plot(rf_rfe, type=c("g", "o")) # plot the results


if (length(predictors(rf_rfe)>10))  {
    SD_rf_rfe=predictors(rf_rfe)[1:10]
} else {
    SD_rf_rfe=predictors(rf_rfe)
}

rm(control,counts_data_filtered)
##ANOTHER IDEA FOR FEATURE SELECTION:
##If the un-normalized data was given, I could have implemented DESeq2 and 
##use the Differentially Expressed Genes as the REDUCED features.


##########################TESTING THE MACHINES
 

###KNN
###Running KNN with the whole feature counts data leads to this error:
###Error: protect(): protection stack overflow
###Although we can change the settings in R-studio,
###I chose to use the selected genes in the previous sections as the reduced features
knn_SD_Rand = train(Simplified_class ~., data = feature_vec[,c(SD_Rand,"Simplified_class")]
            , method = "knn",
             metric="Accuracy",
            preProcess = c("center","scale"),
             trControl = trainControl(method="repeatedcv", number=10, repeats=3)
             )

knn_SD_PCA = train(Simplified_class ~., data = feature_vec[,c(SD_PCA,"Simplified_class")]
                , method = "knn",
                metric="Accuracy",
                preProcess = c("center","scale"),
                trControl = trainControl(method="repeatedcv", number=10, repeats=3)
                )
  
knn_SD_rf_rfe  =  train(Simplified_class ~., data = feature_vec[,c(SD_rf_rfe,"Simplified_class")]
                        , method = "knn",
                        metric="Accuracy",
                        preProcess = c("center","scale"),
                        trControl = trainControl(method="repeatedcv", number=10, repeats=3)
                )


svm_radial_SD_Rand = train( Simplified_class~., data = feature_vec[,c(SD_Rand,"Simplified_class")]
                            , method = "svmRadial", 
                            metric = "Accuracy",
                            trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                            preProcess = c("center","scale"))

svm_linear_SD_Rand = train( Simplified_class~., data = feature_vec[,c(SD_Rand,"Simplified_class")]
             , method = "svmLinear",
             metric = "Accuracy",
              trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
             preProcess = c("center","scale"))


svm_radial_SD_PCA =  train( Simplified_class~., data = feature_vec[,c(SD_PCA,"Simplified_class")]
                            , method = "svmRadial",
                            metric = "Accuracy",
                            trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                            preProcess = c("center","scale"))

svm_linear_SD_PCA =  train( Simplified_class~., data = feature_vec[,c(SD_PCA,"Simplified_class")]
                            , method = "svmLinear",
                            metric = "Accuracy",
                            trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                            preProcess = c("center","scale"))

svm_radial_SD_rf_rfe =  train( Simplified_class~., data = feature_vec[,c(SD_rf_rfe,"Simplified_class")]
                               , method = "svmRadial",
                               metric = "Accuracy",
                               trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                               preProcess = c("center","scale"))

svm_linear_SD_rf_rfe =  train( Simplified_class~., data = feature_vec[,c(SD_rf_rfe,"Simplified_class")]
                               , method = "svmLinear",
                               metric = "Accuracy",
                               trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                               preProcess = c("center","scale"))


rf_SD_Rand = train( Simplified_class~., data = feature_vec[,c(SD_Rand,"Simplified_class")]
       , method = "rf",
       metric = "Accuracy",
       trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
       preProcess = c("center","scale"))

rf_SD_PCA = train( Simplified_class~., data = feature_vec[,c(SD_PCA,"Simplified_class")]
                  , method = "rf",
                  metric = "Accuracy",
                  trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                  preProcess = c("center","scale"))

rf_SD_rf_rfe = train( Simplified_class~., data = feature_vec[,c(SD_rf_rfe,"Simplified_class")]
       , method = "rf",
       metric = "Accuracy",
       trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
       preProcess = c("center","scale"))

###List of all the machines made: (13)
###knn_SD_Rand, knn_SD_PCA, knn_SD_rf_rfe,
###svm_radial_SD_Rand, svm_linear_SD_Rand, svm_radial_SD_PCA,
###svm_linear_SD_PCA, svm_radial_SD_rf_rfe, svm_linear_SD_rf_rfe,
###rf_SD_Rand, rf_SD_PCA, rf_SD_rf_rfe , rf_rfe



###knn machines report
print(knn_SD_Rand)
print(knn_SD_PCA)
print(knn_SD_rf_rfe)

###svm machines report
print(svm_radial_SD_Rand)
print(svm_linear_SD_Rand)
print(svm_radial_SD_PCA)
print(svm_linear_SD_PCA)
print(svm_radial_SD_rf_rfe)
print(svm_linear_SD_rf_rfe)

###random forest machines report
print(rf_SD_Rand)
print(rf_SD_PCA)
print(rf_SD_rf_rfe)

###random forest with rfe report
print (rf_rfe)

###I have reported their performances in the Explanation_A.pdf



###on question b
