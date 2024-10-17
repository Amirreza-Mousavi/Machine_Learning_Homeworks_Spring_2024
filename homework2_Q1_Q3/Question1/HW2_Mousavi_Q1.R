#HW2 Seyyed Amirreza Mousavi Majd
#Student Number = 610799001

####################

###Initialization: load libraries, import datasets, declare functions

library(ggplot2)
library(pheatmap)
library(factoextra)
library(fossil)

pathbased = read.table("pathbased.txt")
spiral = read.table("spiral.txt")
jain = read.table("jain.txt")
flame = read.table("flame.txt")
s1 = read.table("s1.txt")
s4 = read.table("s4.txt")


purity_index_calc = function(predicted_cluster_assignments,real_cluster_assignments){
  
  pred = tabulate(predicted_cluster_assignments)
  real = tabulate(real_cluster_assignments)
  conf_matrix = t( cbind(pred,real))
  purity_index = (max(pred)+max(real))/sum(conf_matrix)
  return(purity_index)
}


pdf("graphical_results.pdf")

#Q1


##for dataset = pathbased k = 3
##for dataset = spiral k = 3
##for dataset = jain k = 2
##for dataset = flame k = 2
##for dataset = s1 k = 15. purity index not applicable.
##for dataset = s4 k = 15. purity index not applicable




do_my_homework_Q1 = function(df,df_K,Is_purity_index_applicable,tag){

print(paste("Analysis on",tag))
rand_index_classifactions_df = as.data.frame(matrix(0,nrow = NROW(df),ncol = 4))
colnames(rand_index_classifactions_df)=c("knn","hc_average","hc_single","hc_complete")

###a&b&c

####i
k1=kmeans(df, centers = df_K, nstart = 25)

if(Is_purity_index_applicable){
print(fviz_cluster(k1,data=df[,-3],main = paste("knn clustering ",tag)))
}
else{
print(fviz_cluster(k1,data=df,main = paste("knn clustering",tag)))#s1 and s4 do not have column #3 i.e. class assignment
}

rand_index_classifactions_df$knn = k1$cluster

if(Is_purity_index_applicable){
print(paste("Purity index is",purity_index_calc(k1$cluster,df$V3)," for knn vs real"))
}

####ii
hc_average_df = hclust(d = dist(df,method = "euclidean"), method = "average")
print(plot(hc_average_df,main = paste("average hc for ",tag),cex = 0.1)) #cex is font size. small to avoid crowded plots.
dummy_storage = rect.hclust(hc_average_df,k=df_K,border = 1:df_K)
#better plot
print(ggplot(df)+
geom_point(aes(x = V1, y = V2, color = as.factor(cutree(hc_average_df,k = df_K))))
)

if(Is_purity_index_applicable){
print(paste("Purity index is",purity_index_calc(cutree(hc_average_df,k = df_K),df$V3)," for hc_average vs real"))
}

rand_index_classifactions_df$hc_average = cutree(hc_average_df,k = df_K)


####iii
hc_single_df = hclust(d = dist(df,method = "euclidean"), method = "single")
print(plot(hc_single_df,main = paste("single hc for",tag),cex = 0.1))
dummy_storage = rect.hclust(hc_single_df,k=df_K,border = 1:df_K)

#better plot
print(ggplot(df)+
geom_point(aes(x = V1, y = V2, color = as.factor( cutree(hc_single_df,k = df_K))))
)

if(Is_purity_index_applicable){
print(paste("Purity index is",purity_index_calc(cutree(hc_single_df,k = df_K),df$V3)," for hc_single vs real"))
}

rand_index_classifactions_df$hc_single = cutree(hc_single_df,k = df_K)


####iv
hc_complete_df = hclust(d = dist(df,method = "euclidean"), method = "complete")
print(plot(hc_complete_df,main = paste("complete hc for",tag), cex = 0.1))
dummy_storage = rect.hclust(hc_complete_df,k=df_K,border = 1:df_K)

#better plot
print(ggplot(df)+
geom_point(aes(x = V1, y = V2, color = as.factor( cutree(hc_complete_df,k = df_K))))
)

if(Is_purity_index_applicable){
print(paste("Purity index is",purity_index_calc(cutree(hc_complete_df,k = df_K),df$V3)," for hc_complete vs real"))
}


rand_index_classifactions_df$hc_complete = cutree(hc_complete_df,k = df_K)


########Rand index summing up
ridf = rand_index_classifactions_df #a name for cleaner codes
print("classification results point by point: head records")
print(head(ridf)) #for clarification

heat_matrix = matrix(0,nrow = 4,ncol = 4)
for(i in 1:4) {
  for(j in 1:4){
    heat_matrix[i,j] = rand.index(ridf[,i],ridf[,j])
  }
}

heat_matrix = as.data.frame(heat_matrix,row.names = colnames(ridf))
colnames(heat_matrix) = rownames(heat_matrix)

print("Rand index matrix is :")
print(heat_matrix)
pheatmap(heat_matrix,main = paste(tag," rand index heatmap"))

}
###############################################################################

################################################################################

########################################################################################

###############################################################################

##############################################################################



##for dataset = pathbased k = 3
##for dataset = spiral k = 3
##for dataset = jain k = 2
##for dataset = flame k = 2
##for dataset = s1 k = 15. purity index not applicable.
##for dataset = s4 k = 15. purity index not applicable


do_my_homework_Q1(df = pathbased,df_K = 3,Is_purity_index_applicable = T,tag = "pathbased")

do_my_homework_Q1(df = spiral,df_K = 3,Is_purity_index_applicable = T,tag = "spiral")

do_my_homework_Q1(df = jain,df_K = 2,Is_purity_index_applicable = T, tag = "jain")

do_my_homework_Q1(df = flame,df_K = 2,Is_purity_index_applicable = T, tag = "flame")

do_my_homework_Q1(df = s1,df_K = 15,Is_purity_index_applicable = F, tag = "s1")

do_my_homework_Q1(df = s4,df_K = 15,Is_purity_index_applicable = F, tag = "s4")

dev.off()
