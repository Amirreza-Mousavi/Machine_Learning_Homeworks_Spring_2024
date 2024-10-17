###Machine Learning Project
###Seyyed Amirreza Mousavi Majd
###610799001

###Machine Learning Project
###Seyyed Amirreza Mousavi Majd
###610799001
###Question C:

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

feature_vec = feature_vec %>% mutate(Patient_ID=NULL)
View(feature_vec)

###Because in question A, filtering out the correlated ones, and performing
###Recursive Feature Elimination yielded the best results, we employ this 
###method again; with a little modifications.


####Feature Selection by Filtering out highly correlated gene expressions
cor_table=cor(feature_vec[,1:17396]) #takes about a minute
fc=findCorrelation(cor_table,cutoff = 0.7,names = F) #to_be removed

feature_vec_filtered = feature_vec[,-fc]
cat(paste("In the first step, we reduced our",dim(feature_vec)[2],"features to"
          ,dim(feature_vec_filtered)[2],"features.","Features are genes & the 4 metadata."))


rm(cor_table,counts_data,meta_data,fc)
###

class_original=(feature_vec_filtered %>% select(Simplified_class))[,1]
feature_vec_filtered = feature_vec_filtered %>% mutate(Simplified_class = NULL)
View(feature_vec_filtered)

decision_layer_one = function(cls_orig){ ###we will use it just for class_original vector
  cls_orig=as.character(cls_orig)
  for (i in 1:length(cls_orig)){
    if(cls_orig[i]=="Non_advanced_Fibrosis"||
       cls_orig[i]=="Advanced_fibrosis"){
      cls_orig[i]="Fibrosis"
    }
  }
  return(as.factor(cls_orig))
}

#print(decision_layer_one(class_original))
class_layer_one = decision_layer_one(class_original)

###takes about 3 minutes
rf_rfe_C_layer1 = rfe(x = feature_vec_filtered,
              y = class_layer_one, 
              size=10, rfeControl= rfeControl(functions=rfFuncs, 
                                              method="cv", number=10))
               
plot(rf_rfe_C_layer1, type=c("g", "o"))
print(rf_rfe_C_layer1)
               

SD_rf_rfe_C_layer1 = predictors(rf_rfe_C_layer1)[1:10] #pick the top ten predictors.

###we will continue the work with SD_rf_rfe_C
###we will build a radial svm and random forest with SD_rf_rfe_C because these
###methods yielded the best results in question A

svm_radial_SD_rf_rfe_C_layer1 =  train(x = feature_vec_filtered%>%select(all_of(SD_rf_rfe_C_layer1)),
                                y = class_layer_one,
                               method = "svmRadial",
                               metric = "Accuracy",
                               trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                               preProcess = c("center","scale"))


rf_SD_rf_rfe_C_layer1 = train(x = feature_vec_filtered%>%select(all_of(SD_rf_rfe_C_layer1)),
                       y = class_layer_one,
                       method = "rf",
                       metric = "Accuracy",
                       trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                       preProcess = c("center","scale"))
                        

print(svm_radial_SD_rf_rfe_C_layer1)
print(rf_SD_rf_rfe_C_layer1)


choose_svm_layer1 = (max(svm_radial_SD_rf_rfe_C_layer1$results[,"Accuracy"]) > max(rf_SD_rf_rfe_C_layer1$results[,"Accuracy"]))
###choose_svm_layer1's value may be dependent on the random seed

rm(rf_rfe_C_layer1)

###Next step: classify the sickness into the advanced and non-advanced.
feature_vec_filtered_layer2=feature_vec_filtered %>% 
                               mutate(class_original) %>% 
                                 filter(class_original != "Normal")
class_layer_two = (feature_vec_filtered_layer2 %>% select(class_original))[,1]
feature_vec_filtered_layer2 = feature_vec_filtered_layer2 %>% mutate(
                            class_original = NULL)
View(feature_vec_filtered_layer2)
### 118 person have remained.
print(class_layer_two)

###let's train the second layer. We will use RFE again

feature_vec_filtered_layer2 = feature_vec_filtered_layer2 %>%
                              mutate(SEX = as.factor(SEX),
                              Diabet = as.factor(Diabet))


class_layer_two = droplevels.factor(class_layer_two)

###takes about 3 minutes
rf_rfe_C_layer2 = rfe(x = feature_vec_filtered_layer2,
                      y = class_layer_two, 
               size=10, rfeControl= rfeControl(functions=rfFuncs, 
                                               method="cv", number=10))


print(rf_rfe_C_layer2)
SD_rf_rfe_C_layer2 = predictors(rf_rfe_C_layer2)[1:10]

plot(rf_rfe_C_layer2, type=c("g", "o"))

rm(rf_rfe_C_layer2)

###Now let's pass the SD_rf_rfe_C_layer2 to 2 different ML machines
###SVM radial and random forest

svm_radial_SD_rf_rfe_C_layer2 =  train(x = feature_vec_filtered_layer2%>%select(all_of(SD_rf_rfe_C_layer2)),
                                       y = class_layer_two,
                                       method = "svmRadial",
                                       metric = "Accuracy",
                                       trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                                       preProcess = c("center","scale"))


rf_SD_rf_rfe_C_layer2 = train(x = feature_vec_filtered_layer2%>%select(all_of(SD_rf_rfe_C_layer2)),
                              y = class_layer_two,
                              method = "rf",
                              metric = "Accuracy",
                              trControl = trainControl(method="repeatedcv", number=10, repeats=3), 
                              preProcess = c("center","scale"))

choose_svm_layer2 = (max(svm_radial_SD_rf_rfe_C_layer2$results[,"Accuracy"]) > max(rf_SD_rf_rfe_C_layer2$results[,"Accuracy"]))



###Now we wrap things up.
#########################################################################
#########################################################################
###Performing Nested_Machine_Learning_Classifier
test_index=sample(1:192,size = 39) #39 = 20% of 192
test_data = feature_vec[test_index,]

if(choose_svm_layer1){
  pred_lay1 = predict(svm_radial_SD_rf_rfe_C_layer1,test_data[,-17401] %>% select(all_of(SD_rf_rfe_C_layer1)))
  
} else{
  pred_lay1 = predict(rf_SD_rf_rfe_C_layer1,test_data[,-17401] %>% select(all_of(SD_rf_rfe_C_layer1)))
  
}

t1 = table(pred_lay1,decision_layer_one(test_data[,"Simplified_class"])) #first layer classification
print(t1)
cat(paste("Accuracy is",sum(diag(t1))/39,"in the first layer of classification."))

test_data=test_data[which(pred_lay1 == "Fibrosis"),] #going for the second layer
test_data$Simplified_class=droplevels.factor(test_data$Simplified_class)


if(choose_svm_layer2){
  pred_lay2 = predict(svm_radial_SD_rf_rfe_C_layer2,test_data[,-17401] %>% select(all_of(SD_rf_rfe_C_layer2)))
} else{
  pred_lay2 = predict(rf_SD_rf_rfe_C_layer2,test_data[,-17401] %>% select(all_of(SD_rf_rfe_C_layer2)))
}

t2 = table(pred_lay2,test_data[,"Simplified_class"])
print(t2)
cat(paste("Accuracy is",sum(diag(t2))/sum(colSums(t2)),"in the second layer of classification"    ))
