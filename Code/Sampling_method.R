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

# DST_sampling : Function to create a stratified random sample
DST_sampling <- function(Grid.roi, IOP.field, SHOP.field) {
  
  ##############################################################################
  # INPUT :                                                                    #
  # - Grid.roi : A grid intersecting the entire ROI. In our case each cell has #
  # a square size of 0.5 arc_degree                                            #
  # - IOP.field : Reference data for industrial oil palm plantations           #
  # - SHOP.field : Reference data for smallholder oil palm plantations         #
  #                                                                            #
  # OUTPUT : A list of sf points object corresponding to a sample of           #
  # 2500 points/grid cell with a "Reference" column corresponding to the data  # 
  # from Gaveau et al., 2022 These 2500 are distributed according to the       #
  # proportion of areas in each class studied.                                 #
  ##############################################################################
  
  # Display the processing code for each cell
  cat("\n Sampling cell = ", i, " / ", nrow(Grid.roi),"\n")
  
  # Display the processing code for each cell
  grid.tmp = Grid.roi[i,]
  Area.grid= as.numeric(st_area(grid.tmp))
  
  # Intersect IOP and SHOP fields with the grid cell and calculate their areas
  IOP.Gaveau.grid = st_union(st_intersection(IOP.field, grid.tmp))
  Area.IOP = as.numeric(st_area(IOP.Gaveau.grid))
  
  #shop selection and area
  SHOP.Gaveau.grid = st_union(st_intersection(SHOP.field, grid.tmp))
  Area.SHOP = as.numeric(st_area(SHOP.Gaveau.grid))
  
  # Calculate proportions of IOP, SHOP, and Other
  Prop.IND = as.numeric(100*sum(Area.IOP)/Area.grid)
  Prop.SH = as.numeric(100*sum(Area.SHOP)/Area.grid)
  Prop.Other = as.numeric(100 - (Prop.IND + Prop.SH))
  
  # Determine sample size allocation based on proportions
  N.sample = 2500
  N.IOP.Grid = round(N.sample*(Prop.IND/100))
  N.SHOP.Grid = round(N.sample*(Prop.SH/100))
  N.OLU.Grid = round(N.sample*(Prop.Other/100))
  
  # Generate random points within the grid cell
  All.sample.pts.grid = st_sample(grid.tmp, 10000, type="random")
  
  # Identify which points fall into IOP and SHOP areas  
  IOP.sample.pts.grid = st_intersects(All.sample.pts.grid, IOP.Gaveau.grid)
  SHOP.sample.pts.grid = st_intersects(All.sample.pts.grid, SHOP.Gaveau.grid)
  
  # Convert intersect results to binary vectors
  IOP.sample.pts.grid.vc = sapply(IOP.sample.pts.grid, function(x) length(x))
  SHOP.sample.pts.grid.vc = sapply(SHOP.sample.pts.grid, function(x) length(x))
  
  # Classify points based on overlaps
  Grid.Gaveau.OLU = which(IOP.sample.pts.grid.vc == 0 & SHOP.sample.pts.grid.vc == 0)
  Grid.Gaveau.IOP = which(IOP.sample.pts.grid.vc == 1)
  Grid.Gaveau.SHOP = which(SHOP.sample.pts.grid.vc == 1)
  
  # Extract classified points  
  Grid.Gaveau.OLU.pts = All.sample.pts.grid[Grid.Gaveau.OLU]
  Grid.Gaveau.IOP.pts = All.sample.pts.grid[Grid.Gaveau.IOP]
  Grid.Gaveau.SHOP.pts = All.sample.pts.grid[Grid.Gaveau.SHOP]
  
  # Sample points for each class and convert to sf objects with metadata
  
  # OLU
  if(N.OLU.Grid > 0){
    Grid.Gaveau.OLU.pts.grid = sample(Grid.Gaveau.OLU.pts, as.integer(N.OLU.Grid))
    OLU.Gaveau.samples.pts.sf.Grid = st_as_sf(Grid.Gaveau.OLU.pts.grid)
    OLU.Gaveau.samples.pts.sf.Grid$Ref = rep(0,nrow(OLU.Gaveau.samples.pts.sf.Grid)) 
    OLU.Gaveau.samples.pts.sf.Grid$Grid = i
  }else{
    Grid.Gaveau.OLU.pts.grid = sample(Grid.Gaveau.OLU.pts, 0)
    OLU.Gaveau.samples.pts.sf.Grid = st_as_sf(Grid.Gaveau.OLU.pts.grid)
    OLU.Gaveau.samples.pts.sf.Grid$Ref = rep(0,nrow(OLU.Gaveau.samples.pts.sf.Grid)) 
  }
  
  # IOP
  if(N.IOP.Grid > 0){
    Grid.Gaveau.IOP.pts.grid = sample(Grid.Gaveau.IOP.pts, as.integer(N.IOP.Grid))
    IOP.Gaveau.samples.pts.sf.Grid = st_as_sf(Grid.Gaveau.IOP.pts.grid)
    IOP.Gaveau.samples.pts.sf.Grid$Ref = rep(1,nrow(IOP.Gaveau.samples.pts.sf.Grid)) 
    IOP.Gaveau.samples.pts.sf.Grid$Grid = i
  }else{
    Grid.Gaveau.IOP.pts.grid = sample(Grid.Gaveau.IOP.pts, 0)
    IOP.Gaveau.samples.pts.sf.Grid = st_as_sf(Grid.Gaveau.IOP.pts.grid)
    IOP.Gaveau.samples.pts.sf.Grid$Ref = rep(1,nrow(IOP.Gaveau.samples.pts.sf.Grid)) 
  }
  
  # SHOP
  if(N.SHOP.Grid > 0){
    Grid.Gaveau.SHOP.pts.grid = sample(Grid.Gaveau.SHOP.pts, as.integer(N.SHOP.Grid))
    SHOP.Gaveau.samples.pts.sf.Grid = st_as_sf(Grid.Gaveau.SHOP.pts.grid)
    SHOP.Gaveau.samples.pts.sf.Grid$Ref = rep(2,nrow(SHOP.Gaveau.samples.pts.sf.Grid)) 
    SHOP.Gaveau.samples.pts.sf.Grid$Grid = i
  }else{
    Grid.Gaveau.SHOP.pts.grid = sample(Grid.Gaveau.SHOP.pts, 0)
    SHOP.Gaveau.samples.pts.sf.Grid = st_as_sf(Grid.Gaveau.SHOP.pts.grid)
    SHOP.Gaveau.samples.pts.sf.Grid$Ref = rep(2,nrow(SHOP.Gaveau.samples.pts.sf.Grid)) 
  }
  
  # Combine all sampled points into one data frame
  All.Gaveau.samples.pts.sf.Grid = rbind(OLU.Gaveau.samples.pts.sf.Grid, 
                                         IOP.Gaveau.samples.pts.sf.Grid, 
                                         SHOP.Gaveau.samples.pts.sf.Grid
  )
  # Add Binarized reference column (0 = NOP ; 1 = OP)
  All.Gaveau.samples.pts.sf.Grid$Ref_bin = All.Gaveau.samples.pts.sf.Grid$Ref 
  All.Gaveau.samples.pts.sf.Grid$Ref_bin[All.Gaveau.samples.pts.sf.Grid$Ref_bin == 2]=1
  
  # Ouput : Sample with Ref, Ref_bin colmuns
  return(All.Gaveau.samples.pts.sf.Grid)
}

# Prepare_sample : Function to prepare training and testing samples from spatial data
Prepare_sample = function(Ls.sample, Grid.roi){
  
  ##############################################################################
  # INPUT :                                                                    #
  # - Ls.sample : A list of sample for each cell grid derived from             #
  # DST_sampling function                                                      #
  # - Grid.roi : A grid intersecting the entire ROI. In our case each cell has #
  # a square size of 0.5 arc_degree                                            #
  #                                                                            #
  # OUTPUT : A list of training and testing sample for each dataset used.      #
  # The sample was divided in 70/30 train and test data.                       #
  ##############################################################################
  
  # Display the processing code for each cell
  cat("\n Extract cell = ", i, " / ", nrow(Grid.roi),"\n")
  
  # Extract the current sample from the list of samples
  Sample.tmp = Ls.sample[[i]]
  
  # Retrieve the list of raster files to extract tile ID from file names 
  MB.ls.files = list.files("./DST/Grid/MB/Sumatra_renamed/")
  ID.tiles = gsub("[^0-9]", "", MB.ls.files[i])
  
  # Load the corresponding raster layers for each dataset
  MB.tmp = rast(paste0("./DST/Grid/MB/Sumatra_renamed/MB_",ID.tiles,".tif"))
  D.tmp = rast(paste0("./DST/Grid/Descals/Sumatra_renamed/Descals_",ID.tiles,".tif"))
  I.tmp = rast(paste0("./DST/Grid/IIASA/Sumatra_renamed/IIASA_",ID.tiles,".tif"))
  X.tmp = rast(paste0("./DST/Grid/Xu/Sumatra_renamed/XU_", ID.tiles, ".tif"))
  # Here it is possible to add or remove sources (Proba.DE)
  
  # Extract data for each dataset and remove NA 
  
  # MAPBIOMAS
  Sample.MB = extract(MB.tmp, Sample.tmp)
  Sample.MB$classification_2019[is.na(Sample.MB$classification_2019)] = 0
  
  # DESCALS
  Sample.D = extract(D.tmp, Sample.tmp)
  Sample.D$Descals_raw[is.na(Sample.D$Descals_raw)] = 0
  
  #IIASA
  Sample.I = extract(I.tmp, Sample.tmp)
  Sample.I$Sumatra_raw_IIASA[is.na(Sample.I$Sumatra_raw_IIASA)] = 0
  
  # XU
  Sample.X = extract(X.tmp, Sample.tmp)
  Sample.X$`2016_op`[is.na(Sample.X$`2016_op`)] = 0
  
  # Combine the extracted values with the original sample
  Sample.tmp = cbind(Sample.tmp, 
                     Sample.D$Descals_raw, 
                     Sample.I$Sumatra_raw_IIASA, 
                     Sample.X$`2016_op`, 
                     Sample.MB$classification_2019
  )
  
  # Rename the columns to have descriptive names
  names(Sample.tmp)[names(Sample.tmp) != "x"] <- c("Ref", 
                                                   "Grid", 
                                                   "Ref_bin", 
                                                   "D", 
                                                   "I", 
                                                   "X", 
                                                   "MB"
  )
  
  # Set seed for reproducibility when splitting data
  set.seed(123)  
  
  # Calculate the sample size for 70% of the data
  sample_size <- floor(0.7 * nrow(Sample.tmp))
  
  # Randomly select 70% of the indices for training data
  train_indices <- sample(seq_len(nrow(Sample.tmp)), size = sample_size)
  
  # Split the data into training (70%) and testing (30%) sets
  Data.train <- Sample.tmp[train_indices, ]  # Training data
  Data.test <- Sample.tmp[-train_indices, ]  # Testing data
  
  # Return the training and testing datasets as a list
  return(list(Sample.train = Data.train, Sample.test = Data.test))
}

################################################################################
#                                MAIN                                          #
################################################################################

# Initialization of lists to store data at various stages
sample.ls <- list()           # List to store samples for each grid
Prep.sample.ls <- list()      # List to store prepared samples
Data.test.ls <- list()        # List to store test samples
Data.train.ls <- list()       # List to store training samples

# Load required spatial data
Grid.roi <- st_read("./DST/Grid/Shape/Grid_test.shp")               # Grid regions of interest
IOP.field <- st_read("./Gaveau/Data/Validation/Ind_OP_Sumatra.shp")  # Field data for "Ind_OP"
SHOP.field <- st_read("./Gaveau/Data/Validation/SH_OP_Sumatra.shp")  # Field data for "SH_OP"

# Loop through each region of interest (ROI) in the grid
for (i in 1:nrow(Grid.roi)) {
  set.seed(i)  # Ensure reproducibility for sampling
  sample.ls[[i]] <- DST_sampling(Grid.roi, IOP.field, SHOP.field)    # Perform sampling
  
  set.seed(i)  # Reset seed for consistency
  Prep.sample.ls[[i]] <- Prepare_sample(sample.ls, Grid.roi)         # Prepare the sample
  
  set.seed(i)  # Reset seed again
  Data.test.ls[[i]] <- Prep.sample.ls[[i]]$Sample.test               # Extract test data
  Data.train.ls[[i]] <- Prep.sample.ls[[i]]$Sample.train             # Extract training data
}

# Combine individual test and training samples into unified datasets
train.sample <- do.call(rbind, Data.train.ls)  # Merge all training samples
test.sample <- do.call(rbind, Data.test.ls)    # Merge all test samples

# Optionally save the datasets as shapefiles (uncomment if needed)
# st_write(train.sample, "./DST/Sample_train_mdpi.shp")
# st_write(test.sample, "./DST/Sample_test_mdpi.shp")