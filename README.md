# Cell-Painting
This image analysis pipeline was developed to differentiate heart cells during their transition from iPSCs to cardiomyocytes. The goal was to extract as many features as possible from a set of 4000 images taken of a 96 well plate.

## More about each file:
### Analysis w mask.cpipe: 
CellProfiler pipeline with module added to remove midlayer of cells
### Analysis.cpipe: 
CellProfiler pipeline for general feature extraction
### AnalysisQC.cpipe: 
CellProfiler pipeline for extracting quality control features at image level
### Illumniation.cpipe: 
Calculates and exports illumination correction function
### Python_Analysis.ipynb: 
PCA and other data analysis 
### R Analysis.R: 
Data analysis
### R Analysis.ipynb: 
same as R but converted to notebook format for better readability 
