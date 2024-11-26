################################################################################
#                                  LIBRARIES                                   #
################################################################################
# Load necessary libraries for spatial and statistical operations
library(sf)      # For handling spatial data (v.1.0.16)
library(terra)   # For raster data manipulation (v.1.7.55)
library(pracma)  # For numerical and mathematical functions (v.2.3.8)

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
  
  return(list("ConfMat" = ConfMat, "OverAcc" = OverAcc, "Kappa" = Kappa, "ProdAcc" = ProdAcc, "UserAcc" = UserAcc, "F1-score" = F.score, "Cell"=i))
}

# Mass_assign : # Function to calculate mass assignments based on classification validation results
Mass_assign = function(Valid.Classif){
  
  ##############################################################################
  # INPUT :                                                                    #
  # - Valid.Classif : A list object derived from Valid_Classif function        #
  #                                                                            #
  # OUTPUT : A dataframe of masses needed for DST application                  #
  ##############################################################################
  
  # Extract the confusion matrix from the validation results
  CM.tmp = Valid.Classif$ConfMat
  
  # Extract true positives (TP), true negatives (TN), false negatives (FN), and false positives (FP)
  TP = CM.tmp[2,2] # Bottom-right cell: True Positives
  TN = CM.tmp[1,1] # Top-left cell: True Negatives
  FN = CM.tmp[2,1] # Bottom-left cell: False Negatives
  FP = CM.tmp[1,2] # Top-right cell: False Positives
  
  # Extract the Kappa statistic from the validation results
  Kappa.tmp = Valid.Classif$Kappa
  
  # Create a data frame to store mass assignment values
  # Columns: NOP (No Oil Palm), OP (Oil Palm), Unc (Uncertainty)
  source.df = data.frame(NOP = c(0,	0), OP = c(0,	0), Unc=c(0,0))
  
  # Row 1: Mass assignments for the "Oil Palm" class (OP)
  source.df[1,1] = Kappa.tmp # Probability of correctly identifying OP (Kappa score)
  source.df[1,2] = (1-Kappa.tmp)*(FN/(FN+TP)) # Misclassification of OP as NOP
  source.df[1,3] = 1-(source.df[1,1]+source.df[1,2]) # Remaining probability assigned to uncertainty
  
  # Row 2: Mass assignments for the "No Oil Palm" class (NOP)
  source.df[2,1] = (1-Kappa.tmp)*(FP/(FP+TN)) # Misclassification of NOP as OP
  source.df[2,2] = Kappa.tmp # Probability of correctly identifying NOP (Kappa score)
  source.df[2,3] = 1-(source.df[2,1]+source.df[2,2]) # Remaining probability assigned to uncertainty
  
  # Return the data frame containing mass assignments
  return(source.df)
}

#Combin_process : Function to process combinations using Dempster-Shafer Theory (DST)
Combin_process = function(Nb.class, Nb.sources, List.sources, sampling){
  
  ##############################################################################
  # INPUT :                                                                    #
  # - Nb.class : Number of class available in the question  (OP, NOP)          #
  # - Nb.sources : Number of sources available (DESCALS, IIASA, XU, MB)        #
  # - List.sources : A list object derived from Mass_assign function and       #
  # corresponding to the mass assignment for each sources                      #
  #                                                                            #
  # OUTPUT : A list of sf point object corresponding to sampling data and      #
  # a dataframe of all available combinations                                  #
  ##############################################################################
  
  # Generate all possible combinations of 0 and 1 for the given number of sources
  Combi.df <- expand.grid(rep(list(0:(Nb.class - 1)), Nb.sources))
  names(Combi.df) = c("D", "I", "X", "MB") # Naming the columns for clarity
  
  # Initialize vectors to store computed probabilities for each combination
  DR.OP = numeric(nrow(Combi.df))
  DR.NOP = numeric(nrow(Combi.df))
  DST.UNC = numeric(nrow(Combi.df))
  DST.CONF = numeric(nrow(Combi.df))
  
  ############################## DST rule ######################################
  
  # Iterate through each combination
  for (i in 1:nrow(Combi.df)){
    
    # Initialize source1
    source1 = List.sources[[1]]
    # Define the mass to be chosen depending on the pixel value
    if (Combi.df[i,1] == 0){ row1 = 1} # NOP
    if (Combi.df[i,1] == 1){ row1 = 2} # OP
    source1 = source1[row1,]
    
    # Iterate through remaining sources (starting from the second source)
    for (j in 2:Num.sources){
      
      # Set new source
      source2 = List.sources[[j]]
      
      # Define the mass to be chosen depending on the pixel value
      if (Combi.df[i,j] == 0){ row2 = 1} # NOP
      if (Combi.df[i,j] == 1){ row2 = 2} # OP
      source2 = source2[row2,]
      
      # Calculate the conflict between the source
      conflit = 1-(source1[1]*source2[2]+source2[1]*source1[2])
      
      # Compute new masses based on DST rules
      m.NOP = (source1[1]*source2[1] + source1[1]*source2[3] + source2[1]*source1[3])/conflit
      m.OP  = (source1[2]*source2[2] + source1[2]*source2[3] + source2[2]*source1[3])/conflit  
      m.unc = (source1[3] * source2[3])/conflit
      
      # Update source1 with the new masses for further combination
      source1 = c(m.NOP,m.OP,m.unc)
    }
    
    # Decision making 
    Bel.OP   = m.OP # OP Belief function    
    Plaus.OP = m.OP + m.unc # OP Plausibility function
    Pign.OP = m.OP + 0.5*m.unc # OP Pignistic function
    PCR.OP  = m.OP + (1-conflit) + m.unc # OP PCR modified
    PCR.NOP = m.NOP # NOP probabilities
    
    # Store results in the vectors
    DR.OP[i] = as.numeric(PCR.OP)
    DR.NOP[i] = as.numeric(PCR.NOP)
    DST.UNC[i] = as.numeric(m.unc)
    DST.CONF[i] = as.numeric(1-conflit)
  }
  
  # Add computed probabilities to the combinations data frame
  Combi.df$DR.OP = DR.OP
  Combi.df$DR.NOP = DR.NOP
  Combi.df$Unc = DST.UNC
  Combi.df$Conflict = DST.CONF
  
  ############################## Sampling Processing ###########################
  
  # Initialize vectors to store probabilities for all samples
  DR.OP.all.vc = numeric(nrow(sampling))
  DR.NOP.all.vc = numeric(nrow(sampling))
  DST.Unc.all.vc = numeric(nrow(sampling))
  DST.Conf.all.vc = numeric(nrow(sampling))
  
  # Match the combinations with sample data and assign probabilities
  for (i in 1:nrow(Combi.df)){
    ind.tmp = which(sampling$D == Combi.df$D[i]
                    & sampling$I == Combi.df$I[i]
                    & sampling$X == Combi.df$X[i]
                    & sampling$MB == Combi.df$M[i])
    
    DR.OP.all.vc[ind.tmp] = Combi.df$DR.OP[i]
    DR.NOP.all.vc[ind.tmp] = Combi.df$DR.NOP[i]
    DST.Unc.all.vc[ind.tmp] = Combi.df$Unc[i]
    DST.Conf.all.vc[ind.tmp] = Combi.df$Conflict[i]
  }
  
  # Add the computed values to the sampling data frame
  sampling$DR.OP.DST_4s = DR.OP.all.vc
  sampling$DR.NOP.DST_4s = DR.NOP.all.vc
  sampling$Uncertainty = DST.Unc.all.vc
  sampling$Conflict = DST.Conf.all.vc
  
  # Return the combinations data frame and updated sampling as a list
  return(list(Combi.df = Combi.df, Sampling = sampling))
}

# Perform_decision : Function to perform decision-making based on DST and validate results
Perform_decision <- function(Test.sample, Applied.sample, Combi.df) {
  
  ##############################################################################
  # INPUT :                                                                    #
  # - Test.sample : Testing data from Prepare_sample function                  #
  # - Applied.sample : Sample to apply the DST and validate the results        #
  # - Combi.df : A dataframe with all available combinations                   #
  #                                                                            #
  # OUTPUT : A list with all validation indices for each dataset               #
  ##############################################################################
  
  # Iterate over unique OP values to compute classification and validation metrics
  for (i in seq_along(val.un.OP)) {
    
    # Reclassify samples based on the current threshold value
    tmp.reclass.vc <- numeric(nrow(Applied.sample))
    tmp.reclass.vc[Applied.sample$DR.OP.DST_4s >= val.un.OP[i]] <- 1
    cat("\n i = ", i, " - Seuil = ", val.un.OP[i], " - Sum = ", sum(tmp.reclass.vc))
    
    # Retrieve combination indices for the current threshold
    D.classif <- which(val.un.OP[i] == Combi.df$DR.OP)
    I.classif <- which(val.un.OP[i] == Combi.df$DR.OP)
    X.classif <- which(val.un.OP[i] == Combi.df$DR.OP)
    M.classif <- which(val.un.OP[i] == Combi.df$DR.OP)
    
    # Concatenate the combinations of D, I, X, and MB for the current threshold
    Combi[i] <- paste0(
      Combi.df[D.classif, ]$D, 
      Combi.df[I.classif, ]$I, 
      Combi.df[X.classif, ]$X, 
      Combi.df[M.classif, ]$M, 
      sep = ""
    )
    
    # Validate classification and compute Kappa value
    tmp.valid <- Valid_Classif(Ref.all, tmp.reclass.vc)
    K.all.vc[i] <- tmp.valid$Kappa
  }
  
  # Compile results into a dataframe
  Result.df.DIXM <- data.frame(
    Decision = val.un.OP,        # Decision thresholds (OP)
    Decision_NOP = val.un.NOP,   # Decision thresholds (NOP)
    Kall = K.all.vc,             # Kappa values
    Combinaison = Combi          # Combinations of D, I, X, MB
  )
  
  # Determine final decisions based on the higher of DR.OP or DR.NOP
  Decision.finale.DIXM <- which(Applied.sample$DR.OP.DST_4s > Applied.sample$DR.NOP.DST_4s)
  Applied.sample$DST.DIXM <- numeric(nrow(Applied.sample))  # Initialize decision column
  Applied.sample$DST.DIXM[Decision.finale.DIXM] <- 1        # Assign final decisions
  
  # Validate the results using the reference sample
  Valid.D <- Valid_Classif(Test.sample$Ref_bin, Test.sample$D)         # Validation for D
  Valid.I <- Valid_Classif(Test.sample$Ref_bin, Test.sample$I)         # Validation for I
  Valid.X <- Valid_Classif(Test.sample$Ref_bin, Test.sample$X)         # Validation for X
  Valid.MB <- Valid_Classif(Test.sample$Ref_bin, Test.sample$MB)       # Validation for MB
  Valid.DST <- Valid_Classif(Test.sample$Ref_bin, Applied.sample$DST.DIXM)  # Validation for DST
  
  # Return validation results as a list
  return(list(
    D = Valid.D,
    I = Valid.I,
    X = Valid.X,
    MB = Valid.MB,
    DST = Valid.DST
  ))
}

################################################################################
#                                MAIN                                          #
################################################################################

# Alternatively, load pre-saved training and testing samples
train.sample <- st_read("./DST/Sample_train_mdpi.shp")
test.sample <- st_read("./DST/Sample_test_mdpi.shp")

# Perform validation for each source (D, I, X, MB) on the training sample
i <- 1  # Initialize index for tracking
Valid.D <- Valid_Classif(train.sample$Ref_bin, train.sample$D)  # Validate source D
Valid.I <- Valid_Classif(train.sample$Ref_bin, train.sample$I)  # Validate source I
Valid.X <- Valid_Classif(train.sample$Ref_bin, train.sample$X)  # Validate source X
Valid.MB <- Valid_Classif(train.sample$Ref_bin, train.sample$MB)  # Validate source MB

# Assign mass functions for each source
Mass.D <- Mass_assign(Valid.D)  # Mass assignment for source D
Mass.I <- Mass_assign(Valid.I)  # Mass assignment for source I
Mass.X <- Mass_assign(Valid.X)  # Mass assignment for source X
Mass.MB <- Mass_assign(Valid.MB)  # Mass assignment for source MB

# Create a list of mass functions for combination processing
List.sources <- list(sourceD = Mass.D, 
                     sourceI = Mass.I, 
                     sourceX = Mass.X, 
                     sourceM = Mass.MB
                     )

Num.sources <- length(List.sources)  # Total number of sources

# Perform the combination process using Dempster-Shafer Theory (DST)
Applied.sample <- Combin_process(2, 4, List.sources, test.sample)
Combi.df <- Applied.sample$Combi.df  # Extract combination dataframe
Applied.sample <- Applied.sample$Sampling  # Update Applied.sample with combined sampling data

# Compute unique decision values for OP (Object Presence) and NOP (No Object Presence)
val.un.OP <- sort(unique(st_drop_geometry(Applied.sample[, 9])), decreasing = TRUE)  # OP thresholds
val.un.NOP <- sort(unique(st_drop_geometry(Applied.sample[, 10])), decreasing = TRUE)  # NOP thresholds

# Initialize variables for validation results
K.all.vc <- numeric(length(val.un.OP))       # To store Kappa values for all thresholds
Ref.all <- numeric(length(Applied.sample$Ref))  # Initialize reference data
Ref.all[Applied.sample$Ref == 1] <- 1  # Map reference values for class 1
Ref.all[Applied.sample$Ref == 2] <- 1  # Map reference values for class 2
Combi <- numeric(length(val.un.OP))         # To store combination results

Res.DST = Perform_decision(test.sample, Applied.sample, Combi.df)

################################################################################
#                             ENDING                                           #
################################################################################
