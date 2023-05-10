#Original code from https://github.com/IBIC/virtualPathology/bin/
# reassemble a mask from a slices directory with an edited slice
##Rniftilib is deprecated, RNifti is a similar library however all functions need to be swapped.
library("RNifti")
library("png")
#library("Rniftilib")

args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Usage: reassemble.R brainimage roi")
} else {
    brainfile=args[1]
    roi = args[2]
}

brainnifti <- readNifti(brainfile)
braindata <- brainnifti[,,]
                                        # zero
braindata[,,] <- 0

editedfile <- Sys.glob(paste("slices/*", roi, ".png", sep=""))
origfile <- gsub(paste("_", roi, sep=""), "", editedfile)

editedpng <- readPNG(editedfile)
if (length(dim(editedpng)) > 2) {
    editedpng <- editedpng[,,1]
}
origpng <- readPNG(origfile)
mask <- origpng -editedpng

rotate <- function(x) t(apply(x,2,rev))
 #binarize mask
mask[mask !=0] <- 1
# flip left and right
y <- dim(mask)[2]
mask <- mask[,c(y:1)]
# rotate
mask <- rotate(mask)

sliceno <- gsub("slices/s.", "", origfile)
sliceno <- as.numeric(gsub(".png", "", sliceno))

braindata[,sliceno-1,] <- mask
braindata[,sliceno-2,] <- mask
braindata[,sliceno,] <- mask
braindata[,sliceno+1,] <- mask
braindata[,sliceno+2,] <- mask
braindata[,sliceno+3,] <- mask
braindata[,sliceno+4,] <- mask
braindata[,sliceno+5,] <- mask

nim <- asNifti(brainnifti)
nim$dim <- (niftiHeader(braindata))$dim
nim[,,] <- braindata
filename <- paste(roi, "_mask.nii.gz",sep="")
writeNifti(nim, filename, datatype = "int16")

