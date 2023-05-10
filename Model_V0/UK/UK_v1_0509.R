#install.packages("sp")
#install.packages("gstat")

library(sp)
library(gstat)

# Read the input file into a dataframe
df <- read.delim("input_file.txt", header=TRUE, stringsAsFactors=FALSE)

# View the first few rows of the dataframe
head(df)

# Create a SpatialPointsDataFrame with 3D coordinates
data <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6), z = c(7, 8, 9))
coordinates(data) <- c("x", "y", "z")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")

# Define the variogram model
vgm_model <- vgm(psill = 1, model = "Exp", range = 5000, nugget = 0.1)

# Fit the variogram to the data
variogram_obj <- variogram(z ~ 1, data)
fit_variogram <- fit.variogram(variogram_obj, vgm_model)

# Create a prediction grid
grid <- expand.grid(x = 1:10, y = 1:10, z = 1:10)
coordinates(grid) <- c("x", "y", "z")
gridded(grid) <- TRUE

# Perform kriging on the grid
kriging_result <- krige(z ~ 1, data, grid, fit_variogram)

# Visualize the results
library(rgl)
plot3d(kriging_result$var1.pred, type = "s", col = "blue")
