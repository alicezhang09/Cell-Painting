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
library(philentropy)
library(scales)
library(spade)
library(devtools)
library(Rclusterpp)

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

#metadata function
#Parameters: dataframe to reutrn, main dataframe (cells, nuclei, cytoplasm), columns corresponding to df's condition
#Returns: dataframe with only rows of specified condition

metadata<-function(df,main, x, y){
  df<-main%>% 
    filter(Metadata_PositionY==x|Metadata_PositionY==y)
  return(df)
}


#Shapiro- Wilks function
#Parameters: dataframe, threshold
#Returns: dataframe with only columns with W<threshold
shapiro<-function(df, threshold){
  df<- subset(df, select=-c(Metadata_PositionX,Metadata_Series, Metadata_Site))
  df<- df[ ,- as.numeric(which(apply(df, 2, var) == 0))]
  lshap <- lapply(df, shapiro.test)
  lres <- sapply(lshap, `[`, c("statistic","p.value"))
  lres<- as.data.frame(t(lres))
  lres<-lres[lres[,"statistic"]<threshold,]
  return (lres)
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


#separate dataframes based on site

cells_site10<-cells%>%
  filter(ImageNumber==1)
cells_site11<-cells%>%
  filter(ImageNumber==2)
cells_site12<-cells%>%
  filter(ImageNumber==3)
cells_site13<-cells%>%
  filter(ImageNumber==4)
cells_site14<-cells%>%
  filter(ImageNumber==5)
cells_site15<-cells%>%
  filter(ImageNumber==6)
cells_site16<-cells%>%
  filter(ImageNumber==7)
cells_site9<-cells%>%
  filter(ImageNumber==8)

# compute mean, median, variance for each cell feature

cells_mean<- aggregate(.~ImageNumber, FUN=mean, data=cells)
cells_median<- aggregate(.~ImageNumber, FUN=median, data=cells)
cells_var<- aggregate(.~ImageNumber, FUN=var, data=cells)
mean_var<-as.data.frame(sapply(cells_mean, var))

#scale median dataframe to 0-10, create median profile heatmap

ces <- data.frame(lapply(cells_median[,-1], function(x) (x-min(x))/(max(x) - min(x)) * 10))
ces<-cbind(cells_median[,1], ces)
rownames(ces) <- ces[,1]
heatmap.2(as.matrix(ces), scale = "none", dendrogram='none', Rowv=FALSE, Colv=FALSE, 
          col = colorpanel(1000, "white", "blue"), trace = "none", density.info = "none", 
          xlab="median profile", ylab="condition", main= "Median Profiles", margins = c(12,5))

#create histograms of distributions for a few features
Cells_text<-cells$Texture_Correlation_Actin_10_00
hist(Cells_text, breaks=100)
cells_solidity<-cells$AreaShape_Solidity
hist(cells_solidity, breaks=100)
nuclei_neighbors<-nuclei$Neighbors_NumberOfNeighbors_5
hist(nuclei_neighbors, breaks=10)
nuclei_intensity<-nuclei$Intensity_MinIntensity_Actin
hist(nuclei_intensity, breaks=200)

#create distance matrix and heatmap of cell feature means
dists <- dist(cells_mean)
dists<- as.matrix(dists)
dists<-round(dists, 2)
heatmap.2(dists, Rowv=FALSE, Colv=FALSE, col = colorpanel(1000, "red", "white", "blue"),  margins=c(3,3), cellnote=dists,
          notecex=1.0, notecol="black",na.color=par("bg"), main= "Site Distance Heatmap")

#create correlation matrix and heatmap of cell feature means
corr<-cor(t(cells_mean))
corr<-round(as.matrix(corr),6)
heatmap.2(corr, Rowv=FALSE, Colv=FALSE, col = colorpanel(1000, "red", "white", "blue"),  margins=c(3,3), cellnote=corr,         
          notecex=1.0, notecol="black",na.color=par("bg"), main= "Site Correlation Heatmap")


#running Shapiro-Wilks test
Shapiro_control_nuclei<-as.matrix(shapiro(control_nuclei,0.1))
shapiro_control_cells<-as.matrix(shapiro(control_cells, 0.1))
shapiro_control_cytoplasm<-as.matrix(shapiro(control_cyto, 0.1))

shapiro_cells<-as.matrix(shapiro(cells,0.2))
shapiro_nuclei<-as.matrix(shapiro(nuclei, 0.2))
shapiro_cytoplasm<- shapiro(cytoplasm, 0.2)

write.csv(shapiro_control_cytoplasm,"C:\\Users\\xinli\\Desktop\\shapiro_control_cytoplasm.csv", row.names = TRUE)

# to create a heatmap of a correlation matrix: 
y<- average(Rock_cells, DMSO_cells)
corr<-cor(Rock_cells, y)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "blue"), trace = "none", density.info = "none", 
              xlab="DMSO", ylab="IWP", main= "Correlation matrix heatmap", margins = c(12,12))

#export correlation matrix
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
              col = colorpanel(1000, "red", "white", "blue"), trace = "none", density.info = "none", 
              xlab="Cell shape features", ylab="Cell shape features", main= "Shape correlation heatmap", margins = c(12,12))

corr<-cor(cell_intensity_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "blue"), trace = "none", density.info = "none", 
              xlab="Cell Intensity features", ylab="Cell Intensity features", main= "Intensity correlation heatmap", 
              margins = c(14,14))

corr<-cor(cell_texture_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "blue"), trace = "none", density.info = "none", 
              xlab="Cell texture features", ylab="Cell texture features", main= "Texture correlation heatmap", 
              margins = c(12,12))

corr<-cor(cell_microenvironment_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "blue"), trace = "none", density.info = "none",
              main= "Microenvironment correlation heatmap", margins = c(12,12))



#separate cell features based on wavelength
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

#separate nuclei features based on wavelength
indx<- grepl('mito',colnames(nuclei))
nuc_mito_features<-nuclei[indx]
indx<-grepl('DNA', colnames(nuclei))
nuc_DNA_features<- nuclei[indx]
indx<- grepl('Actin',colnames(nuclei))
nuc_actin_features<-nuclei[indx]
indx<-grepl('Nucleoli', colnames(nuclei))
nuc_nucleoli_features<- nuclei[indx]
indx<- grepl('ER',colnames(nuclei))
nuc_ER_features<-nuclei[indx]

#construct correlation matrices and heatmaps to compare different wavelengths
corr<-cor(cell_mito_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "blue"), trace = "none", density.info = "none", 
              xlab="Mitochondria features", ylab="Mitochondria features", main= "Mitochondria correlation heatmap", 
              margins = c(15,15))

corr<-cor(cell_ER_features)
corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
p<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "blue"), trace = "none", density.info = "none", 
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

ex<-as.data.frame(colnames(cells))
write.csv(ex,"C:\\Users\\xinli\\Desktop\\Features.csv", row.names = FALSE)



dists <- pdist(t(DMSO_cells), t(control_cells))
dists<- as.matrix(dists)
heatmap.2(dists, Colv=NA, Rowv=NA, col=cm.colors(256), scale="row", margins=c(3,3), xlab= "DMSO", ylab="Control",
          main= "Feature-Feature Correlation Heatmap")
          
          