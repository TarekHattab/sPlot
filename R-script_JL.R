library(spdep)
library(maptools)
library(rgdal)
library(raster)
library(fBasics)
library(RColorBrewer)
library(ade4)
library(ncdf4)
library(hexbin)

rootpath <- "C:/Users/Admin/Desktop/sPlot2.0 database 08.08.2016"
setwd(rootpath)

# Importing the sPlot header data file containing plot (relevé) information
plot_data <- read.csv("data/sPLOT.csv", sep=";", na.strings=c("", "NA"))
str(plot_data)
plot_data <- plot_data[, c(1:3)]
names(plot_data) <- c("plot_id", "plot_lon", "plot_lat")
any(is.na(plot_data)) # FALSE: all relevés have coordinates
str(plot_data) # 1117382 relevés

# Importing the sPlot species data file that contains the composition and the relative cover of each species for each relevé
sp_data <- read.csv("data/plots_without_unified_cover_scales_20160120a.csv")
str(sp_data) # 24241941 occurrences
length(unique(sp_data$PlotObservationID)) # 1117940 relevés
length(unique(sp_data$PlotObservationID))-length(unique(plot_data$plot_id)) # 558 relevés in sp_data but missing from plot_data
posit <- match(sp_data$PlotObservationID, plot_data$plot_id)
any(is.na(posit)) # TRUE: some relevés (n = 558) in sp_data are not in plot_data
length(sp_data$PlotObservationID[is.na(posit)]) # 11369 rows in sp_data corresponding to the species lists of the 558 relevés missing from plot_data
sp_data <- sp_data[is.finite(posit), ] # 24230572 occurrences
length(unique(sp_data$PlotObservationID)) # 1117382 relevés
posit <- match(sp_data$PlotObservationID, plot_data$plot_id)
any(is.na(posit)) # FALSE: all relevés in sp_data are also in plot_data

# Importing the sPlot taxonomic backbone
load("data/backbone.v.2.splot.try3.Rdata")
ls()
tax_back <- backbone.splot.try3
str(tax_back) # 122901 taxonomic records
length(levels(sp_data$Matched.concept)) # 87424 taxonomic entities in sp_data
posit <- match(sp_data$Matched.concept, tax_back$names.sPlot.TRY)
any(is.na(posit)) # FALSE: all taxonomic entities in sp_data are also in tax_back
length(levels(sp_data$species)) # 60908 species in sp_data (cf. some taxonomic entities like subspecies are aggregated at the species level)
posit <- match(sp_data$species, tax_back$name.short.correct)
any(is.na(posit)) # FALSE: all species in sp_data are also in tax_back

# Importing the TRY species data file that contains trait data for each species
load("data/TRY.all.mean.sd.2.Rdata") 
ls()
trait_data <- TRY.all.mean.sd.2
str(trait_data) # 40791 taxonomic records
length(levels(trait_data$ StandSpeciesName)) # 40790 taxonomic entities in TRY
length(trait_data$StandSpeciesName[is.na(trait_data$StandSpeciesName)]) # 1 taxonomic record is NA
trait_data <- trait_data[is.finite(trait_data$StandSpeciesName), ]
str(trait_data) # 40790 taxonomic records
posit <- match(trait_data$StandSpeciesName, tax_back$names.sPlot.TRY)
any(is.na(posit)) # TRUE: some taxonomic entities (n = 3958) in trait_data are not in the taxonomic backbone of sPlot
length(trait_data$StandSpeciesName[is.na(posit)]) # 3958 taxonomic entities
trait_data$sp_name <- tax_back$name.short.correct[posit] # Add a new column with the taxonomic name aggregated at the species level
trait_data$sp_name <- as.factor(trait_data$sp_name)
length(levels(trait_data$sp_name)) # 36811 species with trait data

# Symplifying sp_data and cleaning typing mistakes in a few species names
sp_data <- sp_data[, c(1, 3, 8, 12, 13)]
names(sp_data) <- c("plot_id", "tax_grp", "veg_lyr", "sp_name", "rel_cov")
sp_data$veg_lyr <- as.factor(sp_data$veg_lyr)
head(levels(sp_data$sp_name), 20)
length(which(sp_data$sp_name==" Astragalus")) # 283 occurrences
length(which(sp_data$sp_name=="Astragalus")) # 845 occurrences
length(which(sp_data$sp_name=="\"Algae\"")) # 905 occurrences
length(which(sp_data$sp_name=="Algae")) # 0 occurrence
sp_data$sp_name <- as.character(sp_data$sp_name)
sp_data[which(sp_data$sp_name==" Astragalus"), "sp_name"] <- "Astragalus"
sp_data[which(sp_data$sp_name=="\"Algae\""), "sp_name"] <- "Algae"
sp_data$sp_name <- as.factor(sp_data$sp_name)
length(which(sp_data$sp_name=="Astragalus")) # 1128 occurrences
length(which(sp_data$sp_name=="Algae")) # 905 occurrences
any(is.na(sp_data$sp_name)) # TRUE: some species (n = 20376) are unidentified 
length(sp_data$sp_name[is.na(sp_data$sp_name)]) # 20376 non identified species
sp_data <- sp_data[is.finite(posit), ] # 21879465 occurrences
length(unique(sp_data$plot_id)) #  1113833 relevés with taxa all identified at the species level
posit <- match(sp_data$sp_name, trait_data$sp_name)
length(posit[is.finite(posit)]) # 18988700 occurrences out of 21879465 (i.e. about 87%) in sp_data with trait data

# Adjusting plot_data to delete all relevés that are only made of unidentified species names 
posit <- match(plot_data$plot_id, sp_data$plot_id)
any(is.na(posit)) # TRUE: some relevés (n = 3549) in plot_data are not in sp_data anymore
length(plot_data$plot_id[is.na(posit)]) # 3549 relevés (1113833+3549=1117382)
plot_data <- plot_data[is.finite(posit), ] # 1113833 relevés
posit <- match(plot_data$plot_id, sp_data$plot_id)
any(is.na(posit)) # FALSE: all relevés in plot_data are also in sp_data

# Computing species richness for each relevé
comm_list <- split(as.character(sp_data$sp_name), sp_data$plot_id)
plot_data$sp_rich <- as.double(lapply(comm_list, function (x) length(unique(x))))
any(is.na(plot_data$sp_rich)) # FALSE: all relevés in plot_data have at least 1 species in the community
min(plot_data$sp_rich) # 1 species in the community
max(plot_data$sp_rich) # 662 species in the community: THIS IS HUGE!!!
hist(plot_data$sp_rich) # There is clearly some outlier relevés in sPLot: how do they look like?
posit <- which(plot_data$sp_rich==max(plot_data$sp_rich)) # row 83295 
plot_data$plot_id[posit] # Relevé number 83657
comm_list[posit] # Looks like a very rich tropical plot but how big is this relevé (cf. its surface area)?
ramp_hexbin <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
hexbinplot(sp_rich~plot_lat, data=plot_data, trans=log, inv=exp, xbins=30, border="white", xlab="Latitude (°)", ylab="Species richness", colramp=ramp_hexbin) # Looks like a mid-domain effect!!!

# Computing community weighted means (CWMs) for some trait data
posit <- match(sp_data$sp_name, trait_data$sp_name)
sp_data$SLA_mean <- trait_data$SLA.mean[posit]
sp_data$PH_mean <- trait_data$PlantHeight.mean[posit]
sp_data$SM_mean <- trait_data$SeedMass.mean[posit]
trait_list <- split(sp_data$SLA_mean, sp_data$plot_id)
weight_list <- split(sp_data$rel_cov, sp_data$plot_id) 
trait_list <- Map("*", trait_list, weight_list)
plot_data$CWM_seff <- as.double(lapply(trait_list, function (x) length(x[is.finite(x)])))
length(which(plot_data$CWM_seff<3)) #  68325 relevés for which the CWM value is based on less that three species
length(which(plot_data$CWM_seff>2)) #  1045508 relevés for which the CWM value is based on more that two species
plot_data$SLA_mean <- as.double(lapply(trait_list, function (x) mean(x, na.rm=TRUE)))



# CAUTION: WORK ON WEIGHTED MEAN FORMULA
trait_list <- split(sp_data$PH_mean, sp_data$plot_id)
trait_list <- Map("*", trait_list, weight_list)
plot_data$PH_mean <- as.double(lapply(trait_list, function (x) mean(x, na.rm=TRUE)))
trait_list <- split(sp_data$SD_mean, sp_data$plot_id)
trait_list <- Map("*", trait_list, weight_list)
plot_data$SD_mean <- as.double(lapply(trait_list, function (x) mean(x, na.rm=TRUE)))

# Checking for relevés with identical coordinates (just for information) 
plot_data$lonlat_id <- paste(plot_data$plot_lon, plot_data$plot_lat, sep=":")
length(unique(plot_data$lonlat_id)) # Only 523435 out of 1113833 (i.e. about 47%) have unique coordinates (cf. pseudoreplication in space is high but could potentially be different communities)

# SHALL WE DO THAT: subsampling (randomly) one relevé out of those sharing exactly the same coordinates?
# uniq <- unique(plot_data$lonlat_id)
# uniq_plot_list <- c()
# for (i in 1:length(uniq)){
#   posit_i <- which(plot_data$lonlat_id==uniq[i])
#   uniq_plot_list <- c(uniq_plot_list, as.character(sample(plot_data$plot_id[posit_i], size=1)))
# }
# posit <- which(as.character(plot_data$plot_id) %in% uniq_plot_list)
# plot_data <- plot_data[posit, ]

# Making plot_data a spatial data frame
CRSlonlat <- CRS("+proj=longlat +datum=WGS84")
coords <- cbind(plot_data$plot_lon, plot_data$plot_lat)
coords <- SpatialPoints(coords, proj4string=CRSlonlat)
plot_data <- SpatialPointsDataFrame(coords, plot_data, proj4string=CRSlonlat)
# Plot the data on the world map
data(wrld_simpl)
plot(wrld_simpl)
points(coordinates(plot_data), cex=0.1, col="red", pch=16)



# cartogramme 
library(Rcartogram)
library(getcartr)
library(fftw)

proj4string(wrld_simpl) <- proj4string(plot_data)
overRes <- over(plot_data,wrld_simpl)
NbplotCountry <- table(overRes$ISO2)
IdCountry <- match(wrld_simpl@data$ISO2,names(NbplotCountry))
wrld_simpl@data$NbPlot <- NbplotCountry[IdCountry]
wrld_simpl@proj4string <- CRS(as.character(NA))

tiff("carto_sPlot1.tiff",width = 1000, height = 1000, pointsize = 40)
sPlot.carto.1<- quick.carto(wrld_simpl, wrld_simpl@data$NbPlot+1, blur = 1)
spplot(sPlot.carto.1["NbPlot"])
dev.off()

tiff("carto_sPlot2.tiff",width = 1000, height = 1000, pointsize = 40)
sPlot.carto.2<- quick.carto(wrld_simpl, 1/(wrld_simpl@data$NbPlot+1), blur = 1)
spplot(sPlot.carto.2["NbPlot"])
dev.off()


# Importing the functional biome classification from Higgins et al. (2016)
setwd(paste(rootpath, "/Higgins_al_GCB_functional_biomes", sep = ""))
biom_nc <- nc_open("veg_mean.nc")
print(biom_nc)
lon <- ncvar_get(biom_nc, "longitude")
lat <- ncvar_get(biom_nc, "latitude")
biom <- ncvar_get(biom_nc, "layer")
str(lon)
str(lat)
str(biom)

# Plotting the world terrestrial biomes
r <- raster(nrows=length(lat), ncols=length(lon), xmn=-180, xmx=180, ymn=-90, ymx=90)
biom_r <- setValues(r, as.vector(biom))
ramp_biom <- c(brewer.pal(3, "Greys"), brewer.pal(3, "Blues"), brewer.pal(6, "YlOrRd"), brewer.pal(6, "RdPu"), rev(brewer.pal(6, "BrBG")[1:3]), brewer.pal(3, "Greens"))
code_biom <- c("SLC", "SMC", "SHC", "TLC", "TMC", "THC", "SLD", "SMD", "SHD", "TLD", "TMD", "THD", "SLB", "SMB", "SHB", "TLB", "TMB", "THB", "SLN", "SMN", "SHN", "TLN", "TMN", "THN")
plot(biom_r, col=ramp_biom, legend=FALSE, main="Functional biomes from Higgins et al. (2016)")
plot(biom_r, col=ramp_biom, legend.only=TRUE, legend.width=1, legend.shrink=0.75, axis.args=list(at=seq(1, 24, length.out=24), labels=code_biom, cex.axis=0.6), legend.args=list(text="Biomes", side=3, font=2, line=0.4, cex=0.8))

# Extracting biome information for each plot in the sPlot database
plot_data$biom_class <- extract(biom_r, coordinates(plot_data))
plot_data@data$biom_class <- as.factor(plot_data@data$biom_class)
str(plot_data@data)
levels(plot_data@data$biom_class) # All terrestrial biomes are represented in sPlot except SHC but SHC is among the least frequent terrestrial biomes on Earth
plot_data@data$biom_name <- NA
for (i in 1:length(code_biom)){
  plot_data@data[which(plot_data@data$biom_class==i), "biom_name"] <- code_biom[i]
}
plot_data@data$biom_name <- as.factor(plot_data@data$biom_name)
str(plot_data@data)

# Plotting the sampling effort per terrestrial biome next to the relative freqeuncy of biome classes
par(mfrow=c(1, 2))
barplot(table(factor(plot_data@data$biom_class, levels=1:24)), col=ramp_biom, names.arg=code_biom, las=1, horiz=TRUE)
title(main="Number (n) of plots \nfrom sPlot per biome class")
barplot(biom_r, col=ramp_biom, names.arg=code_biom, las=1, horiz=TRUE, axes=FALSE)
title(main="Relative frequency \nof biome classes on Earth")
par(mfrow=c(1, 1))

# Importing bioclim grids at 10 arc-minutes resolution and resampling each grid at 30-m resolution to match the grid of biome classes
setwd(paste(rootpath, "/data/bio_10m", sep = ""))
file_list <- dir()[grep(".bil", dir())]
file_list <- file_list[c(1, 12:19, 2:11)]
bioclim_r <- raster(file_list[1])
bioclim_r <- resample(bioclim_r, biom_r, "bilinear")
bioclim_r <- mask(bioclim_r, biom_r)
bioclim_wrld <- xFromCell(bioclim_r, cell=c(1:ncell(bioclim_r)))
bioclim_wrld <- as.data.frame(bioclim_wrld)
names(bioclim_wrld) <- c("lon_30m")
bioclim_wrld$lat_30m <- yFromCell(bioclim_r, cell=c(1:ncell(bioclim_r)))
bioclim_wrld$bio1_30m <- getValues(bioclim_r)
for (i in 2:length(file_list)){
  print(paste(i, "of", length(file_list), sep = " "))
  temp <- raster(file_list[i])
  temp <- resample(temp, biom_r, "bilinear")
  temp <- mask(temp, biom_r)
  bioclim_wrld[, paste(strsplit(file_list, ".bil")[[i]], "_30m", sep="")] <- getValues(temp)
}
bioclim_wrld$biom_class <- getValues(biom_r) 
str(bioclim_wrld) # 259200 cells
bioclim_wrld <- na.omit(bioclim_wrld)
str(bioclim_wrld) # 62608 cells
bioclim_wrld$biom_class <- as.factor(bioclim_wrld$biom_class)
bioclim_wrld$biom_name <- NA
for (i in 1:length(code_biom)){
  bioclim_wrld[which(bioclim_wrld$biom_class==i), "biom_name"] <- code_biom[i]
}
bioclim_wrld$biom_name <- as.factor(bioclim_wrld$biom_name)
str(bioclim_wrld)

# Making it a spatial data frame
coords <- cbind(bioclim_wrld$lon_30m, bioclim_wrld$lat_30m)
coords <- SpatialPoints(coords, proj4string=CRSlonlat)
bioclim_wrld <- SpatialPointsDataFrame(coords, bioclim_wrld, proj4string=CRSlonlat)

# Running a principal component analysis (PCA) based on bioclim variables
pca_wrld <- dudi.pca(bioclim_wrld@data[, c(3:21)], scannf=FALSE, nf=3)
sum(((pca_wrld$eig/sum(pca_wrld$eig))*100)[1:2]) # 77% loaded on the first two PCA axes
sum(((pca_wrld$eig/sum(pca_wrld$eig))*100)[1:3]) # 85% loaded on the first three PCA axes
bioclim_wrld$pc1_val <- pca_wrld$li$Axis1
bioclim_wrld$pc2_val <- pca_wrld$li$Axis2
bioclim_wrld$pc3_val <- pca_wrld$li$Axis3
str(bioclim_wrld@data)

# Mapping the two first PCA axes
temp <- getValues(bioclim_r)
temp[which(is.finite(temp))] <- bioclim_wrld@data$pc1_val
pc1_r <- setValues(bioclim_r, temp)
ramp_pc1 <- divPalette(n=100, name=c("RdYlBu"))
breaks_pc1 <- quantile(pc1_r, probs=seq(0, 1, 0.01))
temp[which(is.finite(temp))] <- bioclim_wrld@data$pc2_val
pc2_r <- setValues(bioclim_r, temp)
ramp_pc2 <- colorRampPalette(brewer.pal(11, "Spectral"))
breaks_pc2 <- quantile(pc2_r, probs=seq(0, 1, 0.01))
par(mfrow=c(2, 1))
plot(pc1_r, breaks=breaks_pc1, col=ramp_pc1, main="First PCA axis on bioclim variables at 30 arc-minute resolution", legend=FALSE)
plot(pc1_r, col=ramp_pc1, legend.only=TRUE, legend.width=1, legend.shrink=0.75, axis.args=list(at=seq(min(breaks_pc1), max(breaks_pc1), length.out=10), cex.axis=0.6), legend.args=list(text="PC1", side=3, font=2, line=0.4, cex=0.8))
plot(pc2_r, breaks=breaks_pc2, col=rev(ramp_pc2(100)), main="Second PCA axis on bioclim variables at 30 arc-minute resolution", legend=FALSE)
plot(pc2_r, col=rev(ramp_pc2(100)), legend.only=TRUE, legend.width=1, legend.shrink=0.75, axis.args=list(at=seq(min(breaks_pc2), max(breaks_pc2), length.out=11), cex.axis=0.6), legend.args=list(text="PC2", side=3, font=2, line=0.4, cex=0.8))
par(mfrow=c(1, 1))

# Plotting the convex hulls for each of the biom classes in the PC1-PC2 space 
convex_hull <- function(xcoord, ycoord, fillcol, borcol=NA){
  hpts <- chull(x=xcoord, y=ycoord)
  hpts <- c(hpts, hpts[1])
  polygon(xcoord[hpts], ycoord[hpts], col=fillcol, border=borcol)
}
plot(bioclim_wrld@data$pc1_val, bioclim_wrld@data$pc2_val, asp=0, type="n", main="Biomes", xlab="PC1", ylab="PC2") 
for (i in rev(1:length(code_biom))){
  sel <- which(bioclim_wrld@data$biom_name==code_biom[i])
  convex_hull(bioclim_wrld@data$pc1_val[sel], bioclim_wrld@data$pc2_val[sel], fill= adjustcolor(ramp_biom[i], alpha.f=0.2), borcol=ramp_biom[i])
}
legend("bottomright", code_biom, fill=adjustcolor(ramp_biom, alpha.f=0.2), border=ramp_biom, bty="n", ncol=4)

# Plotting the density of spatial units (cf. 30 arc-minute resoultion here) available worldwide and the density of sPlot relevés within each bioclimatic unit of the PC1-PC2 space
r <- raster(nrows=100, ncols=100, xmn=min(bioclim_wrld@data$pc1_val), xmx=max(bioclim_wrld@data$pc1_val), ymn=min(bioclim_wrld@data$pc2_val), ymx=max(bioclim_wrld@data$pc2_val))
pca_wrld_r <- rasterize(bioclim_wrld@data[, c("pc1_val", "pc2_val")], r, fun="count")
plot(log(pca_wrld_r), asp=0, col=rev(divPalette(n=100, name="RdBu")), main="Number (n) of 30 arc-minute cells per bin (log(n))", xlab="PC1", ylab="PC2")
convex_hull(bioclim_wrld@data$pc1_val, bioclim_wrld@data$pc2_val, fill= adjustcolor("grey", alpha.f=0.2), borcol="grey") 
plot_data$pc1_val <- extract(pc1_r, coordinates(plot_data))
plot_data$pc2_val <- extract(pc2_r, coordinates(plot_data)) 
pca_sPlot_r <- rasterize(plot_data@data[, c("pc1_val", "pc2_val")], r, fun="count")
windows()
plot(log(pca_sPlot_r), asp=0, col=rev(divPalette(n=100, name="RdBu")), main="Number (n) of sPlots relevés per bin (log(n))", xlab="PC1", ylab="PC2")
convex_hull(bioclim_wrld@data$pc1_val, bioclim_wrld@data$pc2_val, fill= adjustcolor("grey", alpha.f=0.2), borcol="grey") 

# Plotting for each bioclimatic unit of the PC1-PC2 space the ratio between the relative proportion of sPlot relevés and the relative proportion of spatial units available worldwide 
prop_wrld_r <- pca_wrld_r/max(getValues(pca_wrld_r), na.rm=T)
prop_sPlot_r <- pca_sPlot_r/max(getValues(pca_sPlot_r), na.rm=T)
plot(ratio_sPlot/ratio_wrld, asp=0, col=rev(divPalette(n=100, name="RdBu")), main="Oversampled (>1) and undersampled (<1) bioclimatic units", xlab="PC1", ylab="PC2")
convex_hull(bioclim_wrld@data$pc1_val, bioclim_wrld@data$pc2_val, fill= adjustcolor("grey", alpha.f=0.2), borcol="grey") 

# Sensitivity analysis of the effect of the resolution of the PC1-PC2 space on the sampling effort statistics 
res <- seq(10, 500, 10)
ncell_disp <- c()
ncell_samp <- c()
seff_med <- c()
seff_mean <- c()
seff_max <- c()
seff_min <- c() 
for (i in 1:length(res)){
  print(paste(i, "of", length(res), sep = " "))
  r <- raster(nrows=res[i], ncols=res[i], xmn=min(bioclim_wrld@data$pc1_val), xmx=max(bioclim_wrld@data$pc1_val), ymn=min(bioclim_wrld@data$pc2_val), ymx=max(bioclim_wrld@data$pc2_val))
  temp <- rasterize(bioclim_wrld@data[, c("pc1_val", "pc2_val")], r, fun="count")
  ncell_disp <- c(ncell_disp, length(which(getValues(temp)>0)))
  temp <- rasterize(plot_data@data[, c("pc1_val", "pc2_val")], r, fun="count")
  ncell_samp <- c(ncell_samp, length(which(getValues(temp)>0)))
  seff_med <- c(seff_med, median(getValues(temp), na.rm=TRUE))
  seff_mean <- c(seff_mean, mean(getValues(temp), na.rm=TRUE))
  seff_max <- c(seff_max, max(getValues(temp), na.rm=TRUE))
  seff_min <- c(seff_min, min(getValues(temp), na.rm=TRUE))
}
par(mfrow=c(2, 2))
plot(res, seff_med)
plot(res, seff_max)
plot(res, seff_mean)
plot(res, ncell_samp/ncell_disp, ylim=c(0, 1))
par(mfrow=c(1, 1))

# Using the median statistic as the threshold for subsampling each grid cell of the PC1-PC2 space (100*100)

library(Rcpp)
library(bigmemory)
library(RcppArmadillo)
library(RcppParallel)
	
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")

# Build C++ fonctions
setwd("C:/Users/Admin/Desktop/sPlot2.0 database 08.08.2016/Cpp functions")
sourceCpp("bray.part.OpenMP.cpp")
sourceCpp("bray.part.C_RcppParallel.cpp")
sourceCpp("hcr.C.cpp")
sourceCpp("cast.cpp")

# R wrapper for BigBray C++ function
BigBrayPart <- function(bigMat){
    zeros <- big.matrix(nrow = nrow(bigMat),
                        ncol = nrow(bigMat),
                        init = 0,
                        type = typeof(bigMat),
			shared=FALSE,
			backingfile=paste("BrayMatrix_",i,sep=""),
			backingpath=getwd(),
			descriptorfile=paste("BrayMatrix_",i,".desc",sep=""))
    bray_distance_OpenMP(bigMat@address, zeros@address)
    return(zeros)
}

plot_data <- plot_data@data

r <- raster(nrows=100, ncols=100, xmn=min(bioclim_wrld@data$pc1_val), xmx=max(bioclim_wrld@data$pc1_val), ymn=min(bioclim_wrld@data$pc2_val), ymx=max(bioclim_wrld@data$pc2_val))
pca_sPlot_r <- rasterize(plot_data[, c("pc1_val", "pc2_val")], r, fun="count")

cutoff<-median(values(pca_sPlot_r),na.rm=T)

tempZoneOut <- coordinates(pca_sPlot_r) [which(values(pca_sPlot_r)>cutoff), ]
  plotToRemove <- NULL
  for (i in 1:nrow(tempZoneOut)){
  sel.plot <- which(plot_data$pc1_val > tempZoneOut[i,1]-(res(r)[1]/2) & 
                    plot_data$pc1_val < tempZoneOut[i,1]+(res(r)[1]/2) &
                    plot_data$pc2_val > tempZoneOut[i,2]-(res(r)[2]/2) & 
                    plot_data$pc2_val < tempZoneOut[i,2] +(res(r)[2]/2))
    
  idZoneOut <- plot_data[sel.plot, "plot_id"]
  sel.comm <- sp_data[which(sp_data$plot_id %in% idZoneOut),c("plot_id", "sp_name", "rel_cov")]
  sel.comm<-na.omit(sel.comm)
  sel.comm[,2]<-factor(sel.comm[,2], labels=seq(1:length(unique(sel.comm[,2])))) 
  sel.comm[,2]<-as.numeric(sel.comm[,2])
  comm.data <- castC(iD=sel.comm[,1], sp=sel.comm[,2],cov=sel.comm[,3])
  rowNames <- comm.data[,1]
  comm.data <- comm.data[,-1]
  gc()
   if (nrow(comm.data)>9000) {bigComMatrix <- as.big.matrix(comm.data,shared=FALSE,backingfile=paste("Matrix_",i,sep=""),backingpath=getwd(),descriptorfile=paste("Matrix_",i,".desc",sep="")) ; brayBalDist <- BigBrayPart(bigComMatrix)}
  else { brayBalDist <- bray_distance_RcppParallel(comm.data); brayBalDist <- as.big.matrix(brayBalDist)}
  selectedPlot <- HcrCPP(brayBalDist@address, nout=cutoff, nsampl=1000)  
  selectedPlot <- rowNames[selectedPlot] 
  selectedPlotIndex <- which( idZoneOut %in%  selectedPlot)
  plotToRemove <-  c(plotToRemove,idZoneOut[-selectedPlotIndex])
  print(paste(round(i/nrow(tempZoneOut)*100,1),"%","    i =" ,i,sep=""))
  output<-list(i,plotToRemove)
  save(output ,file="plotToRemove.RData")	
  }
  
  







