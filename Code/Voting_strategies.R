################################################################################
#                                  LIBRARIES                                   #
################################################################################
# Load necessary libraries for spatial and statistical operations
library(sf)      # For handling spatial data (v.1.0.16)
library(terra)   # For raster data manipulation (v.1.7.55)

# Code released in Nov 2024. Check library versions.

################################################################################
#                             SYSTEM PARAMETERS                                #
################################################################################
# Disable spherical geometry for sf operations to simplify spatial calculations
sf::sf_use_s2(FALSE) 

# Set the working directory
WD <- "F:/DST"
setwd(WD)

################################################################################
#                                FUNCTIONS                                     #
################################################################################

# Valid_Classif : Validation function
Valid_Classif = function(LifeForm.field.vc, LifeForm.Classif.vc){
  
  ##############################################################################
  # INPUT :                                                                    #
  # - LifeForm.field.vc : A vector corresponding to reference data             #
  # - LifeForm.Classif.vc : A vector correspinding to classif data             #
  #                                                                            #
  # OUTPUT : A list validation indices                                         #
  ##############################################################################
  
  # Compute confusion matrix
  ConfMat      = table(LifeForm.field.vc,LifeForm.Classif.vc)
  GroundNames  = dimnames(ConfMat)[[1]]
  ClassifNames = dimnames(ConfMat)[[2]]
  
  # Adding missing classes (columns)
  MissingColumns = GroundNames[!(GroundNames %in% ClassifNames)]
  zeros          = matrix(0, ncol = length(MissingColumns), nrow = nrow(ConfMat), dimnames=list(GroundNames, MissingColumns))
  ConfMat        = cbind(ConfMat, zeros)
  ClassifNames   = dimnames(ConfMat)[[2]]
  
  # Adding missing classes (rows)
  MissingRows    = ClassifNames[!(ClassifNames %in% GroundNames)]
  zeros          = matrix(0, nrow = length(MissingRows), ncol = ncol(ConfMat), dimnames=list(MissingRows, ClassifNames))
  ConfMat        = rbind(ConfMat, zeros)
  GroundNames    = dimnames(ConfMat)[[1]]
  
  # Reordering confusion matrix
  ConfMat       = ConfMat[sort(GroundNames),]
  ConfMat       = ConfMat[,sort(GroundNames)]
  
  # Basic stats
  N      = sum(ConfMat)
  SumRow = rowSums(ConfMat)
  SumCol = colSums(ConfMat)
  
  # Extraction of diagonal
  Diag   = sapply(1:ncol(ConfMat),function(i){ConfMat[i,i]})
  
  # Overall accuracy
  OverAcc = sum(Diag)/N
  
  # Precision and recall
  Precision = sapply(1:ncol(ConfMat),function(i){Diag[i]/(SumCol[i])})
  Recall = sapply(1:ncol(ConfMat),function(i){Diag[i]/(SumRow[i])})
  
  # Kappa index
  A     = sum(SumRow * SumCol)
  Kappa = (N*sum(Diag) - A)/ (N^2-A)
  
  # User's accuracy and producer's accuracy
  ProdAcc = Diag/SumRow
  UserAcc = Diag/SumCol
  
  #compute F score for each class
  F.score= 2*(Precision*Recall)/(Precision+Recall)
  
  return(list("ConfMat" = ConfMat, "OverAcc" = OverAcc, "Kappa" = Kappa, "ProdAcc" = ProdAcc, "UserAcc" = UserAcc, "F1-score" = F.score))
}

#Voting_fusion : Function to combine map classifications by voting strategies
Voting_fusion <- function(classifications, weight){
  
  ##############################################################################
  # INPUT :                                                                    #
  # - classifications : an sf object corresponding to a sample of points       #
  # and providing for each point a class according to the classifications      #
  # to be merged                                                               #
  #                                                                            #
  # OUTPUT : An sf object providing a classifications according to voting      # 
  # strategies approach                                                        #
  ##############################################################################
  
  # Retrieve and store sf object geometry 
  geometry.classif <- st_geometry(classifications)
  
  # Remove geometry of the sf object to get a dataframe
  classifications = st_drop_geometry(classifications)
  
  # Check the dimensions between input vectors
  if (length(classifications) != length(weight)) {
    stop("The number of classifications must correspond to the number of weights.")
  }
  
  # Normalize weights (sum = 1)
  weights.norm <- weight / sum(weight)
  
  # Weighted fusion
  fusion <- rowSums(mapply(function(classe, p) classe * p, classifications, weights.norm))
  
  # Application of decision rule (majority)
  # If >= 0.5, class OP (1), otherwise NOP (0)
  fusionned.classif <- ifelse(fusion >= 0.5, 1, 0)
  
  # Creation of an sf object with merged classification and corresponding geometries
  fusionned.classif <- st_sf(classif = fusionned.classif, geometry = geometry.classif)
  
  return(fusionned.classif)
}

################################################################################
#                                MAIN                                          #
################################################################################

# Load required spatial data
train.points = st_read("./DST/Sample_train_mdpi.shp")
test.points = st_read("./DST/Sample_test_mdpi.shp")

# Validation process for maps to be fused
Valid.D <- Valid_Classif(test.points$Ref_bin, test.points$D)         # Validation for D
Valid.I <- Valid_Classif(test.points$Ref_bin, test.points$I)         # Validation for I
Valid.X <- Valid_Classif(test.points$Ref_bin, test.points$X)         # Validation for X
Valid.MB <- Valid_Classif(test.points$Ref_bin, test.points$MB)       # Validation for MB

# Setting weights
weight.equal.vc <- c(1, 1, 1, 1)  
weight.Kappa.vc = c(Valid.D$Kappa, Valid.I$Kappa, Valid.X$Kappa, Valid.MB$Kappa)

# Application of the Voting_fusion function
Applied.equal.weighted = Voting_fusion(test.points[,c(4:7)], weight.equal.vc)
Applied.kappa.weighted = Voting_fusion(test.points[,c(4:7)], weight.Kappa.vc)

# Statistical validation of results
Valid.egal.weighted<- Valid_Classif(test.points$Ref_bin, Applied.equal.weighted$classif)      
Valid.Kappa.weighted<- Valid_Classif(test.points$Ref_bin, Applied.kappa.weighted$classif)  
