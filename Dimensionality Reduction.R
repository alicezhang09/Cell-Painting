# load libraries
library(readr)
library(dplyr)
library(gplots)
library(pdist)
library (ggplot2)
library(purrr)
library(lineup)
library(ppcor)
library(caret)

#function to import dataframes and remove rows and columns with mostly empty cells
#Parameter: list of file names
#Returns: dataframe combining all files from list, excluding columns and rows more than half empty

import<-function(x){
  x<-x%>%map_dfr(read.csv)
  x<-x[rowSums(is.na(x))<ncol(x)*0.9,colSums(is.na(x))<nrow(x)*0.9]
  return (x)
}


#normalize function
#Parameter: dataframe to normalize
#Returns: normalized dataframe obtained by subtracting the mean of the variable and dividing the result 
#by the variable's standard deviation, first five columns are ignored because they contain information
#about the metadata of the features, and are not features themselves

normalize<-function(x){
  x<- subset(x, select=-c(Metadata_PositionX,Metadata_Series, Metadata_Site))
  x <- x[ - as.numeric(which(apply(x, 2, var) == 0))]
  x<-x[,colSums(is.na(x))<nrow(x)]
  x<-x[rowSums(is.na(x))<ncol(x)*0.5,]
  x1<-x[,c(1:5)]
  x2<- as.data.frame(scale(x[,c(6:ncol(x))]))
  x<-cbind(x1, x2)
  return(x)
}

#Parameter: 2 dataframes to compare
#Returns shorter dataframe of two parameters

get_shorter<-function(x,y){
  if(nrow(x)>nrow(y)){
    return(y)
  }
  else {
    return(x)
  }
}

#Parameter: 2 dataframes to compare
#Returns longer dataframe of two parameters

get_longer<-function(x,y){
  if(nrow(x)>nrow(y)){
    return(x)
  }
  else{
    return(y)
  }
}

metadata<-function(df,main, x, y){
  df<-main%>% 
    filter(Metadata_PositionY==x|Metadata_PositionY==y)
  return(df)
}

#function that adds rows with the column's average so the shorter dataframe matches the dimensions of the 
#longer dataframe
#Parameter: two dataframes with the same number of columns
#Returns: the originally shorter dataframe with the additional rows
average<-function(longer,shorter){
  average<-as.data.frame.list(colMeans(shorter))
  average<-average%>% slice(rep(1, each=nrow(longer)-nrow(shorter)))
  shorter2<-rbind(average, shorter)
  return(shorter2)
}


#function that produces heatmap of correlation matrix between two dataframes 
#Parameters: two datframes with equal dimensions
#Returns: heatmap of correlation matrix (Method: Pearson)
heatmap_function<-function(x,y){
  corr<-cor(x,y)
  corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
  heatmap1<<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
                        col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
                        xlab=x, ylab=y, main= "Correlation matrix heatmap", margins = c(12,12))
}

heatmap1<-function(x){
  heatmap<<-heatmap.2(x, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
                      col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
                      xlab=x, ylab=x, main= "Correlation matrix heatmap", margins = c(12,12))
}

#import feature dataframes and remove rows and columns with more than half empty cells

setwd("C:/Users/xinli/Desktop/Work/Spreadsheets/Features")
cells <- list.files(pattern="*Cells.csv")
cells <- import(cells)
cells<- normalize(cells)
cytoplasm <- list.files(pattern="*Cytoplasm.csv")
cytoplasm <-import(cytoplasm)
cytoplasm<- normalize(cytoplasm)
nuclei <- list.files(pattern="*Nuclei.csv")
nuclei <- import(nuclei)
nuclei<- normalize(nuclei)


#import image dataframes

setwd("C:/Users/xinli/Desktop/Work/Spreadsheets/Image Data")

images<-list.files(pattern="*Image.csv")
images<-images%>% map_dfr(read.csv)

#separate dataframes into experimental conditions 

control_cells<- metadata(control_cells, cells, 2, 11)
control_nuclei<- metadata(control_nuclei,nuclei, 2, 11)
control_cyto<-metadata(control_cyto, cytoplasm, 2, 11)

DMSO_cells<- metadata(DMSO_cells, cells, 3, 10)
DMSO_nuclei <- metadata (DMSO_nuclei, nuclei, 3,10)
DMSO_cytoplasm <- metadata(DMSO_cytoplasm, cytoplasm, 3, 10)

CHIR_cells<-metadata(CHIR_cells, cells, 5, 8)
CHIR_nuclei<-metadata(CHIR_nuclei, nuclei, 5,8)
CHIR_cytoplasm<-metadata (CHIR_cytoplasm, cytoplasm, 5, 8)

Rock_cells<-metadata(Rock_cells, cells, 4, 9)
Rock_nuclei<-metadata(Rock_nuclei, nuclei, 4, 9)
Rock_cytoplasm<-metadata(Rock_cytoplasm, cytoplasm, 4,9)

IWP_cells<-metadata(IWP_cells, cells, 6, 7)
IWP_nuclei<-metadata(IWP_nuclei, nuclei, 6, 7)
IWP_cytoplasm<-metadata(IWP_cytoplasm, cytoplasm, 6, 7)


# to create a heatmap of a correlation matrix: 
y<- average(get_longer(DMSO_cells, CHIR_cells), get_shorter(DMSO_cells, CHIR_cells))
corr<-cor(DMSO_cells, y)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
              xlab="DMSO", ylab="CHIR", main= "Correlation matrix heatmap", margins = c(12,12))

y<- average(Rock_cells, DMSO_cells)
corr<-cor(Rock_cells, y)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
              xlab="DMSO", ylab="IWP", main= "Correlation matrix heatmap", margins = c(12,12))
write.csv(corr,"C:\\Users\\xinli\\Desktop\\RockvsDMSO.csv", row.names = TRUE)

#separate dataframes based on features

indx <- grepl('Area', colnames(cells))
cell_shape_features<-cells[indx]

indx<- grepl('Texture',colnames(cells))
cell_texture_features<-cells[indx]

indx<-grepl('Intensity', colnames(cells))
indx2<-grepl('Radial',colnames(cells))
cell_intensity_features<-cbind(cells[indx2],cells[indx])

indx<-grepl('Neighbor', colnames(cells))
cell_microenvironment_features<-cells[indx]

#construct correlation matrices and heatmaps to compare features within specific categories 
corr<-cor(cell_shape_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
              xlab="Cell shape features", ylab="Cell shape features", main= "Shape correlation heatmap", margins = c(12,12))

corr<-cor(cell_intensity_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
              xlab="Cell Intensity features", ylab="Cell Intensity features", main= "Intensity correlation heatmap", 
              margins = c(14,14))

corr<-cor(cell_texture_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
              xlab="Cell texture features", ylab="Cell texture features", main= "Texture correlation heatmap", 
              margins = c(12,12))

corr<-cor(cell_microenvironment_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none",
              main= "Microenvironment correlation heatmap", margins = c(12,12))



#separate features based on wavelength
indx<- grepl('mito',colnames(cells))
cell_mito_features<-cells[indx]
indx<-grepl('DNA', colnames(cells))
cell_DNA_features<- cells[indx]
indx<- grepl('Actin',colnames(cells))
cell_actin_features<-cells[indx]
indx<-grepl('Nucleoli', colnames(cells))
cell_nucleoli_features<- cells[indx]
indx<- grepl('ER',colnames(cells))
cell_ER_features<-cells[indx]

#construct correlation matrices and heatmaps to compare different wavelengths
corr<-cor(cell_mito_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
              xlab="Mitochondria features", ylab="Mitochondria features", main= "Mitochondria correlation heatmap", 
              margins = c(15,15))

corr<-cor(cell_ER_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
              xlab="ER features", ylab="ER features", main= "ER correlation heatmap", 
              margins = c(15,15))

#separate image data based on metadata

Count<-images[,c(1:9,15:17,126)]
Count<- aggregate(.~Metadata_PositionY, FUN=mean, data=Count)
Count<-Count%>% replace(Count=="8", "CHIR")
Count<-Count%>% replace(Count=="9", "Rock Inhibitor")
Count<-Count%>% replace(Count=="10", "DMSO")
Count<-Count%>% replace(Count=="7", "IWP")
Count<-Count%>% replace(Count=="11", "control")

#plot bar graphs to compare average of features extracted from image data

ggplot(Count, aes(x=Metadata_PositionY, y=AreaOccupied_AreaOccupied_Cells)) + geom_bar(stat="identity") + 
  labs(x="Experimental Condition", y="Area Occupied by Cells", main="Area Occupied vs Experimental Condition")
ggplot(Count, aes(x=Metadata_PositionY, y=AreaOccupied_AreaOccupied_Nuclei)) + geom_bar(stat="identity") + 
  labs(x="Experimental Condition", y="Area Occupied by Nuclei")
ggplot(Count, aes(x=Metadata_PositionY, y=Count_Cells)) + geom_bar(stat="identity") + 
  labs(x="Experimental Condition", y="Cell count")

#export dataframe 
write.csv(corr,"C:\\Users\\xinli\\Desktop\\MyData.csv", row.names = FALSE)

#construct distance matrix
cells_matrix<-matrix(as.numeric(unlist(cells[,-1])),nrow=nrow(cells))
nuclei_matrix<-matrix(as.numeric(unlist(nuclei)),nrow=nrow(nuclei))
cytoplasm_matrix<-matrix(as.numeric(unlist(cytoplasm)),nrow=nrow(cytoplasm))


DMSO_cells_matrix<-matrix(as.numeric(unlist(DMSO_cells)),nrow=nrow(DMSO_cells))
control_cells_matrix<-matrix(as.numeric(unlist(control_cells)),nrow=nrow(control_cells))
DMSO_nuclei_matrix<-matrix(as.numeric(unlist(DMSO_nuclei)),nrow=nrow(DMSO_nuclei))
control_nuclei_matrix<-matrix(as.numeric(unlist(control_nuclei)),nrow=nrow(control_nuclei))

dists <- pdist(t(DMSO_cells), t(control_cells))
dists<- as.matrix(dists)
heatmap.2(dists, Colv=NA, Rowv=NA, col=cm.colors(256), scale="row", margins=c(3,3), xlab= "DMSO", ylab="Control",
          main= "Feature-Feature Correlation Heatmap",xlab="DMSO", ylab="control")

