
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
  x<-x[rowSums(is.na(x))<ncol(x)*0.5,colSums(is.na(x))<nrow(x)*0.5]
  return (x)
}


#normalize function
#Parameter: dataframe to normalize
#Returns: normalized dataframe obtained by subtracting the mean of the variable and dividing the result 
#by the variable's standard deviation, first five columns are ignored because they contain information
#about the metadata of the features, and are not features themselves

normalize<-function(x){
  x<-x[,colSums(is.na(x))<nrow(x)]
  x<-x[rowSums(is.na(x))<ncol(x),]
  x1<-x[,c(1:5)]
  x2<- as.data.frame(scale(x[,c(5:ncol(x))]))
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

#function that adds rows with the column's average so the shorter dataframe matches the dimensions of the 
#longer dataframe
#Parameter: two dataframes with the same number of columns
#Returns: the originally shorter dataframe with the additional rows
average<-function(longer,shorter){
  longer<-normalize(longer)
  shorter<-normalize(shorter)
  average<-as.data.frame.list(colMeans(shorter))
  average<-average%>% slice(rep(1, each=nrow(longer)-nrow(shorter)))
  shorter<-rbind(average, shorter)
  return(shorter)
}

#function that produces heatmap of correlation matrix between two dataframes 
#Parameters: two datframes with equal dimensions
#Returns: heatmap of correlation matrix (Method: Pearson)
heatmap_function<-function(x,y){
  corr<-cor(x,y)
  corr<-corr[rowSums(is.na(corr))<ncol(corr)*0.5,colSums(is.na(corr))<nrow(corr)*0.5]
  heatmap1<- heatmap.2(corr, scale = "none", dendrogram='none', Rowv=TRUE, Colv=TRUE, 
              col = colorpanel(1000, "red", "white", "red"), trace = "none", density.info = "none", 
              xlab=x, ylab=y, main= "Correlation matrix heatmap", margins = c(12,12))
  return (heatmap1)
}

#import feature dataframes and remove rows and columns with more than half empty cells

setwd("C:/Users/xinli/Desktop/Work/Spreadsheets/Features")
cells <- list.files(pattern="*Cells.csv")
cells <- import(cells)
cytoplasm <- list.files(pattern="*Cytoplasm.csv")
cytoplasm <-import(cytoplasm)
nuclei <- list.files(pattern="*Nuclei.csv")
nuclei <- import(nuclei)

#import image dataframes

setwd("C:/Users/xinli/Desktop/Work/Spreadsheets/Image Data")

images<-list.files(pattern="*Image.csv")
images<-images%>% map_dfr(read.csv)

#separate dataframes into experimental conditions and normalize data
control_cells<-cells%>% 
  filter(Metadata_PositionY==2|Metadata_PositionY==11)
control_cytoplasm<-cytoplasm%>% 
  filter(Metadata_PositionY==2|Metadata_PositionY==11)
control_nuclei<-nuclei%>% 
  filter(Metadata_PositionY==2|Metadata_PositionY==11)

DMSO_cells <- cells%>%
  filter(Metadata_PositionY==3|Metadata_PositionY==10)
DMSO_nuclei <- nuclei%>%
  filter(Metadata_PositionY==3|Metadata_PositionY==10)
DMSO_cytoplasm <- cytoplasm%>%
  filter(Metadata_PositionY==3|Metadata_PositionY==10)

CHIR_cells<-cells%>%
  filter(Metadata_PositionY==5|Metadata_PositionY==8)
CHIR_nuclei<-nuclei%>%
  filter(Metadata_PositionY==5|Metadata_PositionY==8)
CHIR_cytoplasm<-cytoplasm%>%
  filter(Metadata_PositionY==5|Metadata_PositionY==8)

Rock_cells<-cells%>%
  filter(Metadata_PositionY==4|Metadata_PositionY==9)
Rock_nuclei<-nuclei%>%
  filter(Metadata_PositionY==4|Metadata_PositionY==9)
Rock_cytoplasm<-cytoplasm%>%
  filter(Metadata_PositionY==4|Metadata_PositionY==9)

IWP_cells<-cells%>%
  filter(Metadata_PositionY==6|Metadata_PositionY==7)
IWP_nuclei<-nuclei%>%
  filter(Metadata_PositionY==6|Metadata_PositionY==7)
IWP_cytoplasm<-cytoplasm%>%
  filter(Metadata_PositionY==6|Metadata_PositionY==7)

#normalize data (doesn't work yet)
#cell_list<-list(control_cells, DMSO_cells, CHIR_cells, Rock_cells)
#lapply(cell_list,normalize)
#nuclei_list<-list(control_nuclei, DMSO_nuclei, CHIR_nuclei, Rock_nuclei)
#lapply(nuclei_list, normalize)
#CHIR_nuclei<-normalize(CHIR_nuclei)
#cytoplasm_list<-list(control_cytoplasm, DMSO_cytoplasm, CHIR_cytoplasm, Rock_cytoplasm)
#lapply(cytoplasm_list,normalize)
#CHIR_cytoplasm<-normalize(CHIR_cytoplasm)

# to create a heatmap of a correlation matrix: 
heatmap1<-heatmap_function(average(get_shorter(DMSO_cells,CHIR_cells),get_longer(DMSO_cells,CHIR_cells)),
                           get_longer(DMSO_cells,CHIR_cells))








#OTHER CODE (WORK IN PROGRESS IGNORE FOR NOW)


#c<-c[, colSums(is.na(c)) < nrow(c) * 0.5]

Count<-images[,c(6:8,116)]

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
        main= "Feature-Feature Correlation Heatmap")

