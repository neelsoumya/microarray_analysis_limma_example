Microarray data analysis example using limma

adapted from https://www.bioconductor.org/help/course-materials/2005/BioC2005/labs/lab01/estrogen/


Usage:
    nohup R --no-save < microarray_analysis_example.R


Installation:

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("affy")
 
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("hgu95av2cdf")

Download data from https://www.bioconductor.org/help/course-materials/2005/BioC2005/labs/lab01/Data/estrogen.zip
We assume that the data are unzipped and saved in the folder "~/Downloads/estrogen"


