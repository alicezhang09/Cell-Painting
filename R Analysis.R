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
library(devtools)

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
  x<- subset(x, select=-c(Metadata_PositionX,Metadata_Series))
  x <- x[ - as.numeric(which(apply(x, 2, var) == 0))]
  x<-x[,colSums(is.na(x))<nrow(x)]
  x<-x[rowSums(is.na(x))<ncol(x)*0.5,]
  #x1<-x[,c(1:5)]
  #x2<- as.data.frame(scale(x[,c(6:ncol(x))]))
  #x<-cbind(x1, x2)
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
  df<- subset(df, select=-c(Metadata_Site))
  df<- df[ ,- as.numeric(which(apply(df, 2, var) == 0))]
  lshap <- lapply(df, shapiro.test)
  lres <- sapply(lshap, `[`, c("statistic","p.value"))
  lres<- as.data.frame(t(lres))
  lres<-lres[lres[,"p.value"]<threshold,]
  return (lres)
}

#Shapiro- Wilks function
#Parameters: dataframe, threshold
#Returns: dataframe with only columns with W>threshold
shapiro2<-function(df, threshold){
  df<- subset(df, select=-c(Metadata_Site))
  df<- df[ ,- as.numeric(which(apply(df, 2, var) == 0))]
  lshap <- lapply(df, shapiro.test)
  lres <- sapply(lshap, `[`, c("statistic","p.value"))
  lres<- as.data.frame(t(lres))
  lres<-lres[lres[,"p.value"]>threshold,]
  return (lres)
}

ks2<-function(df){
  df<- subset(df, select=-c(Metadata_Site))
  df<- df[ ,- as.numeric(which(apply(df, 2, var) == 0))]
  lshap <- lapply(df, function(t,d1){ks.test(df,"pnorm")})
  lres <- sapply(lshap, `[`, c("statistic","p.value"))
  return (lres)
}

#diagnol correlation function
#parameters: longer dataframe, shorter dataframe, threshold for correlation, dataframes must have same columns
#output: a list of all columns which have a correlation coeffecient > threshold

diag2<-function(long, short, threshold){
  y<- average(long, short)
  corr<-cor(long, y)
  diag<-as.data.frame(diag(corr))
  diag<-subset(diag, diag$`diag(corr)`>threshold)
  return(diag)
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

cells_site8<-cells%>%
  filter(Metadata_Site==8)
cells_site9<-cells%>%
  filter(Metadata_Site==9)
cells_site10<-cells%>%
  filter(Metadata_Site==10)
cells_site11<-cells%>%
  filter(Metadata_Site==11)
cells_site12<-cells%>%
  filter(Metadata_Site==12)
cells_site13<-cells%>%
  filter(Metadata_Site==13)
nuclei_site14<-nuclei%>%
  filter(ImageNumber==14)
nuclei_site15<-nuclei%>%
  filter(ImageNumber==15)
nuclei_site10<-nuclei%>%
  filter(Metadata_Site==10)
nuclei_site11<-nuclei%>%
  filter(Metadata_Site==11)

# compute mean, median, variance for each cell feature

cells_mean<- aggregate(.~Metadata_Site, FUN=mean, data=cells)
cells_median<- aggregate(.~Metadata_Site, FUN=median, data=cells)
cells_var<- aggregate(.~Metadata_Site, FUN=var, data=cells)
mean_var<-as.data.frame(sapply(cells_mean, var))

cells_mean2<- aggregate(.~Metadata_PositionY, FUN=mean, data=cells)
cells_median2<- aggregate(.~Metadata_PositionY, FUN=median, data=cells)
cells_var2<- aggregate(.~Metadata_PositionY, FUN=var, data=cells)
mean_var3<-as.data.frame(sapply(cells_mean2, var))

#compute features with averages that don't vary much across sites 
mean_var2<- mean_var %>%
  rownames_to_column('rownames') %>%
  filter(`sapply(cells_mean, var)`>1000) %>%
  column_to_rownames('rownames')

mean_var5<- mean_var3 %>%
  rownames_to_column('rownames') %>%
  filter(`sapply(cells_mean2, var)`>1000) %>%
  column_to_rownames('rownames')

 #scale median dataframe to 0-10, create median profile heatmap
ces <- data.frame(lapply(cells_median[,-7], function(x) (x-min(x))/(max(x) - min(x)) * 10))
ces<-cbind(cells_median[,7], ces)
rownames(ces) <- cells_median[,7]
heatmap.2(as.matrix(ces), scale = "none", dendrogram='none', Rowv=FALSE, Colv=FALSE, 
          col = colorpanel(1000, "white", "blue"), trace = "none", density.info = "none", 
          xlab="median profile", ylab="condition", main= "Median Profiles", margins = c(12,5))

#find correlation across sites
control_cells_site10<-control_cells%>%
  filter(Metadata_Site==10)
control_cells_site11<-control_cells%>%
  filter(Metadata_Site==11)
control_s10_vs_s11<-diag2(control_cells_site11, control_cells_site10, 0.6)

control_nuclei_site10<-control_nuclei%>%
  filter(Metadata_Site==10)
control_nuclei_site11<-control_nuclei%>%
  filter(Metadata_Site==11)
controlnuc_s10_vs_s11<-diag2(control_nuclei_site11, control_nuclei_site10, 0.6)

#find correlation across conditions
DMSO_vs_control<-diag2(control_cells, DMSO_cells, 0.1)
DMSO_cells_s10<-DMSO_cells%>%
  filter(Metadata_Site==10)
DMSOs10vs_control<-diag2( control_cells_site10, DMSO_cells_s10, 0.1)

#create histograms of distributions for a few features
Cells_text<-cells$Texture_Correlation_Actin_10_00
hist(Cells_text, breaks=100)
cells_solidity<-cells$AreaShape_Solidity
hist(cells_solidity, breaks=100)
nuclei_neighbors<-nuclei$Neighbors_NumberOfNeighbors_5
hist(nuclei_neighbors, breaks=10)
nuclei_intensity<-nuclei$Intensity_MinIntensity_Actin
hist(nuclei_intensity, breaks=200)

cells_text<- control_cells$AreaShape_FormFactor
hist(cells_text, breaks=100)

#create distance matrix and heatmap of cell feature means
dists <- dist(cells_mean)
dists<- as.matrix(dists)
dists<-round(dists, 2)
heatmap.2(dists, Rowv=FALSE, Colv=FALSE, col = colorpanel(1000, "red", "white", "blue"),  margins=c(3,3), cellnote=dists,
          notecex=1.0, notecol="black",na.color=par("bg"), main= "Site Distance Heatmap")

#create correlation matrix and heatmap of cell feature means
corr<-cor(t(cells_mean))
corr<-round(as.matrix(corr),3)
heatmap.2(corr, Rowv=FALSE, Colv=FALSE, col = colorpanel(1000, "red", "white", "blue"),  margins=c(3,3), cellnote=corr,         
          notecex=1.0, notecol="black",na.color=par("bg"), main= "Site Correlation Heatmap")


#running Shapiro-Wilks test
Shapiro_control_nuclei2<-(shapiro(control_nuclei_site11, 0.05))
Shapiro_control_nuclei3<-(shapiro2(control_nuclei_site11, 0.05))
shaprio_s11_nuclei<-shapiro2(nuclei_site11,0.005)

shap_control<-shapiro2(control_cyto[0:5000,],0.0005)
shaprio_s10_nuclei<-shapiro2(nuclei_site10,0.05)

shapiro_control<-shapiro2(control_cells[0:5000,], 0.0000000005)

write.csv(shapiro_control_cytoplasm,"C:\\Users\\xinli\\Desktop\\shapiro_control_cytoplasm.csv", row.names = TRUE)

#running KS test
h<- ks.test(control_cells$AreaShape_FormFactor, "dnorm")

b<-as.data.frame(t(colnames(control_cells)))
for(i in names(control_cells)){
  b[2,i]<-ks.test(control_cells[,i])
}

h<- lapply(control_cells, function(i){ # Evaluate E at each level of each column
  x <- factor(control_cells[,i])
  ks.test(x, "dnorm")
}) 

h<-colKS(control_cells, interpret=TRUE)

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
Count[,"Nuc_Cyto"] <- Count[,4]/Count[,3]

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
ggplot(Count, aes(x=Metadata_PositionY, y=AreaOccupied_AreaOccupied_Cytoplasm)) + geom_bar(stat="identity") + 
  labs(x="Experimental Condition", y="Area Occupied by Cytoplasm")
ggplot(Count, aes(x=Metadata_PositionY, y=Count_Cells)) + geom_bar(stat="identity") + 
  labs(x="Experimental Condition", y="Cell count")
ggplot(Count, aes(x=Metadata_PositionY, y=Nuc_Cyto)) + geom_bar(stat="identity") + 
  labs(x="Experimental Condition", y="Nucleus to Cytoplasm ratio")

CHIR_images<-images%>%
  filter(Metadata_PositionY==8)
DMSO_images<-images%>%
  filter(Metadata_PositionY==10)
control_images<-images%>%
  filter(Metadata_PositionY==11)

images_area <-CHIR_images$AreaOccupied_AreaOccupied_Nuclei
hist(images_area, breaks=100)

control_area<-control_images$AreaOccupied_AreaOccupied_Nuclei
hist(control_area, breaks=6)
nuclei_intensity<-nuclei$Intensity_MinIntensity_Actin
hist(nuclei_intensity, breaks=200)

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

#calculate control feature-feature correlation matrix

corr_cells<-cor(control_cells)
corr_cells<-corr_cells[-2,-2]


dists <- pdist(t(DMSO_cells), t(control_cells))
dists<- as.matrix(dists)
heatmap.2(dists, Colv=NA, Rowv=NA, col=cm.colors(256), scale="row", margins=c(3,3), xlab= "DMSO", ylab="Control",
          main= "Feature-Feature Correlation Heatmap")
          
          
