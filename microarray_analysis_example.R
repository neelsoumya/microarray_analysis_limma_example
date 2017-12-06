###########################################################################################################
# Microarray analysis example using limma
# 
# adapted from https://www.bioconductor.org/help/course-materials/2005/BioC2005/labs/lab01/estrogen/
#
# Usage:
#    nohup R --no-save < microarray_analysis_example.R
#
# Installation:
#
# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("affy")
# 
# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("hgu95av2cdf")
#
# Download data from https://www.bioconductor.org/help/course-materials/2005/BioC2005/labs/lab01/Data/estrogen.zip
#   We assume that the data are unzipped and saved in the folder "~/Downloads/estrogen"
#
###########################################################################################################

#####################################
# Load libraries
#####################################
library(limma)
library(genefilter)
library(annotate)
library(affy)
library(hgu95av2cdf)
library(VennDiagram)


##########################################################################
# Load data
# Download data from https://www.bioconductor.org/help/course-materials/2005/BioC2005/labs/lab01/Data/estrogen.zip
##########################################################################
str_data_folder = "~/Downloads/estrogen"
setwd(str_data_folder)

targets <- readTargets(file = "estrogen.txt", sep = "")
targets


#####################################
# Normalize data
#####################################
abatch = ReadAffy(filenames = targets$filename)
# normalize
eset = rma(abatch)


##########################################################################
# Build the design matrix
# 
# We have four pairs of replicate arrays so we should estimate 
# four parameters in the linear model. There are many valid ways to choose
# a design matrix, but perhaps the simplest is to make each column 
# correspond to a particular treatment combination:.
# The four columns of the matrix correspond to absent10, present10, 
# absent48 and present48, respectively. Another way to specify the design
# matrix is described in the Limma User's Guide.
# This design matrix given above can be computed in R as follows:
##########################################################################
f <- paste(targets$estrogen, targets$time.h, sep = "")
f <- factor(f)
f
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
design


##########################################################################
# Fit the linear model
##########################################################################
fit <- lmFit(object = eset, design = design)
names(fit)


##########################################################################
# Build contrast matrix
#
# The idea now is to use contrasts to make any comparisons of interest
# between the four treatment combination. Contrasts are linear 
# combinations of parameters from the linear model fit.
# where beta.png is a vector of contrasts for gene g.png, C.png is the 
# contrasts matrix, and  alpha.png is a vector of coefficients 
# (estimated log fold changes), obtained from a linear model fit.
# We will estimate three contrasts (so our contrasts matrix will have 
# three columns). The first contrast is an estrogen effect at time 10
# hours, the second as an estrogen effect at time 48 hours and the third
# is the time effect in the absence of estrogen. These are not all the
# comparisons which might have been made.
##########################################################################
cont_matrix <- makeContrasts(e10 = "present10 - absent10", 
                             e48 = "present48 - absent48",
                             time = "absent48 - absent10",
                             levels = design)
cont_matrix
#Contrasts
#Levels      e10 e48 time
#absent10   -1   0   -1
#absent48    0  -1    1
#present10   1   0    0
#present48   0   1    0


##########################################################################
# Extract the linear model fit for the contrasts
##########################################################################
fit_contrasts = contrasts.fit(fit = fit,
                              contrasts = cont_matrix)

# Given a microarray linear model fit, compute moderated t-statistics, 
# moderated F-statistic, and log-odds of differential expression by
# empirical Bayes moderation of the standard errors towards a common value.
fit_contrasts = eBayes(fit = fit_contrasts)


##########################################################################
# Assessing differential expression
##########################################################################

# We now use the function topTable to obtain a list genes differentially 
# expressed between Estrogen-Present and Estrogen-Absent at time 10 hours, 
# followed by a list of genes differentially expressed between 
# Estrogen-Present and Estrogen-Absent at time 48 hours.
colnames(fit_contrasts)
topTable(fit = fit_contrasts, coef = 2, adjust.method = "BH")

# The function decideTests() provides a variety of ways to assign statistical
# significance to the contrasts while controlling for multiple testing.
results = decideTests(fit_contrasts)
summary(results)

##########################################################################
# Venn diagram of  differential expression
##########################################################################

vennDiagram(results)

