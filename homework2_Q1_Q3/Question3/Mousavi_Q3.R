#Homework2 Question3
#Seyyed Amirreza Mousavi Majd
#St. No. -> 610799001

##Initialization
library(factoextra)
library(ggplot2)
anno_df = read.csv("Leaf_Root_annotation.csv")
anno_df[,1]=NULL

xp_raw_df = read.csv("Leaf_Root_raw_data.csv")
xp_raw_df[,1]=NULL

xp_norm_df = read.csv("Leaf_Root_normalized_data.csv")
xp_norm_df[,1]=NULL

get_tissue_from_SRR=function(srr){
  return(anno_df[which(anno_df$Sample==srr),4])
}

get_project_from_SRR=function(srr){
  return(anno_df[which(anno_df$Sample==srr),1])
}


##clustering is performed using kmeans clustering

calc_label_df = function(df, k_of_cluster){
  k_labels=kmeans(t(df),centers = k_of_cluster)$cluster
  
  res_df = data.frame(kmeans_label = k_labels, 
                      project_label = get_project_from_SRR(names(k_labels)),
                      tissue_label = get_tissue_from_SRR(names(k_labels))
  )
  return(res_df)
}

pca_study = function (df,tag){
  
  pca_df=prcomp(t(df))
  
  pc1 = pca_df$x[,"PC1"]
  pc2 = pca_df$x[,"PC2"]
  
  print(ggplot()+
    geom_point(aes(x=pc1,y=pc2,color = get_project_from_SRR(names(pc1)))) + 
    ggtitle(paste("PCA colored by project for",tag)))
  
  print(ggplot()+
    geom_point(aes(x=pc1,y=pc2,color = get_tissue_from_SRR(names(pc1)))) +
    ggtitle(paste("PCA colored by tissue for",tag)))
  
}

##a

print(calc_label_df(xp_raw_df,k_of_cluster = 2))

print(calc_label_df(xp_raw_df,k_of_cluster = 3))

print(calc_label_df(xp_norm_df,k_of_cluster = 2))

print(calc_label_df(xp_norm_df,k_of_cluster = 3))



##b
pdf("graphical_results_q3_b.pdf")
pca_study(xp_raw_df,"raw dataset")

pca_study(xp_norm_df,"normalized dataset")

dev.off()

