#primary code for analysis of Bachman's Warbler genomic population structure

library("adegenet")
library("hierfstat")
library("pegas")
library(vcfR)
library(reshape2)
library(ggplot2)
library(poppr)
library(ape)
library(igraph)
library(dplyr)
library(diveRsity)
library(nlme)
library("PMCMRplus")
library(rstatix)
library(Rcpp)
library(vegan)

#palette
myCol=c("slateblue2", "darkgreen", "goldenrod1",  "palegreen2" )
                    

#################################################################################################

#convert vcf to genlight 
#input vcf file

#full
vcf<- read.vcfR("BAWA_fullSNPset.recode.vcf")

#restricted
#vcf<- read.vcfR("BAWA_restrictedSNPset.recode.vcf")

#to genind format
BAWA_gen<-vcfR2genind(vcf)

#read in population file 
snp_pop=read.csv("snp_pop_full.csv")
#snp_pop=read.csv("snp_pop_restricted.csv")

#assign population
pop(BAWA_gen)=snp_pop$Pop

###############################################################

#DACP

#NO a priori grouping 
#identify number of groups with k means 
grp=find.clusters(BAWA_gen, max.n.clust=5)

#look at group assignment 
table(pop(BAWA_gen), grp$grp)

#cross-validation to determine number of principal components to use
set.seed(999)
BAWA1 <- xvalDapc(tab(BAWA_gen, NA.method = "mean"),  grp$grp)

#check results 
BAWA1[-1]

#n.da is number of populations - 1 
dapc1 <- dapc(BAWA_gen, var.contrib = TRUE, n.pca=1, n.da=1, grp$grp)

#print contents of the object
print.dapc(dapc1)
#summary/useful info 
summary.dapc(dapc1)
#predict individual assignment 
predict.dapc(dapc1)

scatter(dapc1,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

compoplot(dapc1, col=myCol,lab="", ncol=2)

loadingplot(dapc1$var.contr, threshold=quantile(dapc1$var.contr,0.75))

###############################################################
#DAPC: 
#a priori grouping 
#cross-validation to determine number of principal components to use
set.seed(999)
BAWA <- xvalDapc(tab(BAWA_gen, NA.method = "mean"), pop(BAWA_gen))
#check results 
BAWA[-1]

#n.da is number of populations - 1 
dapc <- dapc(BAWA_gen, var.contrib = TRUE, n.pca=14, n.da=2, pop(BAWA_gen))

#print contents of the object
print.dapc(dapc)
#summary/useful info 
summary.dapc(dapc)
#predict individual assignment 
predict.dapc(dapc)

scatter(dapc,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

compoplot(dapc, col=myCol,lab="", ncol=2)

loadingplot(dapc$var.contr, threshold=quantile(dapc2$var.contr,0.75))

#plot of population assignment predict.dapc(dapc1)
assignplot(dapc, subset=1:33)

#################################################################
#PCA
#convert to genlight format
BAWA_snp <- vcfR2genlight(vcf)
pop(BAWA_snp)=snp_pop$Pop


#remove migratory birds for breeding only analysis
#BAWA_snp1 <- BAWA_snp[pop(BAWA_snp) != "Migratory"]
#remove migratory pop
#snp_pop1 = filter(snp_pop, Pop != "Migratory")
#migratory removed for PCA
#pop(BAWA_snp1)=snp_pop1$Pop

#identify number of groups with k means 
grp=find.clusters(BAWA_snp, max.n.clust=5)
#look at group assignment 
table(pop(BAWA_snp), grp$grp)

#look at eigenvalues to pick best # PCs
BAWA.pca <- glPca(BAWA_snp, nf = 3)
#summary
BAWA.pca 

#eigenvalues
barplot(100*BAWA.pca$eig/sum(BAWA.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance explained", line = 2)
title(xlab="Eigenvalues", line = 1)

BAWA.pca.scores <- as.data.frame(BAWA.pca$scores)
BAWA.pca.scores$pop <- pop(BAWA_snp)

#plot on first 2 PCs
p <- ggplot(BAWA.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = myCol) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_classic()
p
ggsave("full_PCA.jpg", width= 6, height = 6, dpi = 600 )

#breeding birds only
grp=find.clusters(BAWA_snp1, max.n.clust=4)
table(pop(BAWA_snp1), grp$grp)

#look at eigenvalues to pick best # PCs
BAWA.pca <- glPca(BAWA_snp1, nf = 3)
#summary
BAWA.pca
#eigenvalues
barplot(100*BAWA.pca$eig/sum(BAWA.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance explained", line = 2)
title(xlab="Eigenvalues", line = 1)

BAWA.pca.scores <- as.data.frame(BAWA.pca$scores)
BAWA.pca.scores$pop <- pop(BAWA_snp1)

#plot on first 2 PCs
p <- ggplot(BAWA.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = myCol) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_classic()
p

################################################

#admixture plots

tbl=read.table("plink.3.Q")
indTable = read.csv("snp_pop_full.csv",
                      col.names = c("ID", "Pop"))

#basic tables
mergedAdmixtureTable = cbind(tbl, indTable)
ordered = mergedAdmixtureTable[order(mergedAdmixtureTable$Pop),]

#define population borders
xlabels = aggregate (1:nrow(ordered), by = list (ordered[, "Pop"]), FUN = mean)

sampleEdges <- aggregate(1:nrow(ordered),
                         by = list(ordered[, "Pop"]), FUN = max)

#plot
barplot(t(as.matrix(subset(ordered))), col=myCol, border=NA, xlab="Population", ylab = "Ancestry", 
       axisnames = FALSE, space = 0)
abline(v = sampleEdges$x, lwd = 2)
axis(1, at = xlabels$x - 0.5, labels = xlabels$Group.1)


###########################################################

#snps into genind for pop stats

#create hierfstat object
BAWA <- genind2hierfstat(BAWA_gen)

#genetic diversity 
div <- summary(BAWA_gen)
div

#plot observed het
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
#plot observed vs expected het
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

#W&C FST and FIS 
wc(BAWA_gen)

#pairwise Fst
genet.dist(BAWA_gen, method = "WC84")
#bootstrapped confidence intervals
boot.ppfst(dat=BAWA_gen,nboot=1000,quant=c(0.025,0.975),diploid=TRUE)

#use hierfstat to get basic stats
basicstat <- basic.stats(BAWA_gen, diploid = TRUE, digits = 1) 
names(basicstat)

boot.ppfis(dat=BAWA_gen,nboot=1000,quant=c(0,1.0),diploid=TRUE)

#compiled as basic.stats: 
ho <- basicstat$Ho
write.csv(x=ho, file = "ho")
hs <- basicstat$Hs
write.csv(x=hs, file = "hs")
fis <- basicstat$Fis
write.csv(x=fis, file = "fis")

#allelic richness per locus and population
ar <- allelic.richness(BAWA_gen,min.n=NULL,diploid=TRUE)
write.csv(x=ar, file = "allelic_richness")

##################################################################################

#diveRsity
#basic statistics (allelic richness etc using either rarefaction or bootstrapping)
#fis: inbreeding coefficient
#ar: allelic richness


basicStats(infile = "genepop.txt", outfile = "basicstats_diversity", fis_ci = TRUE,
           ar_ci = TRUE, fis_boots = 1000, ar_boots = 1000,
           mc_reps = 1000, rarefaction = TRUE, ar_alpha = 0.05,fis_alpha = 0.05)

fastDivPart(infile = "genepop.txt", outfile = "ROST_divPart", gp = 3,pairwise = TRUE, fst = TRUE,bs_locus = FALSE,
            bs_pairwise = TRUE,boots = 100, plot = FALSE,para = FALSE)


#################################################################################

#isolation-by-distance

library(gdsfmt)
library(SNPRelate)
library(ade4)

#converting a PLINK text/binary file to a GDS file

bed.fn <- "C:~/full.bed"
fam.fn <- "C:~/full.fam"
bim.fn <- "C:~/full.bim"

bed.fn <- "C:~/breeding_full.bed"
fam.fn <- "C:~/breeding_full.fam"
bim.fn <- "C:~/breeding_full.bim"

# Convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "BAWA.gds")

# Summary
snpgdsSummary("BAWA.gds")
#The total number of samples: 46 
#The total number of SNPs: 6147 
##SNP genotypes are stored in SNP-major mode (Sample X SNP).
#The number of valid samples: 46 
#The number of biallelic unique SNPs: 0 

# Open the GDS file
genofile <- snpgdsOpen("BAWA.gds")

# Get population information
#   if it is stored in a text file "pop.txt"
pop_code <- scan("snp_pop_full.txt", what=character())

table(pop_code)

#PCA
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=FALSE)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# assume the order of sample IDs is as the same as population codes
head(cbind(sample.id, pop_code))

# Make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#n×n matrix of genome-wide average IBS pairwise identities

ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE)


# individulas in the same population are clustered together
pop.idx <- order(pop_code)

image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16))

#multidimensional scaling analysis on the n×n matrix of genome-wide IBS pairwise distances
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code)

plot(x, y, col=race, xlab = "", ylab = "",
     main = "Multidimensional Scaling Analysis (IBS)")
legend("topleft", legend=levels(race), pch="o", text.col=1:nlevels(race))

#matrix of genetic distances as 1-IBS
ibs_m <- 1-ibs$ibs
ibs.dist <- dist(ibs_m)
as.matrix(ibs.dist)[1:29, 1:29]

#read in location data 
location <- read.table("locality_breeding.csv", sep=",", header=T)
location.dist <- dist(cbind(location$lat, location$lon))
as.matrix(location.dist)[1:29, 1:29]

mantel.randtest(ibs.dist, location.dist, nrepet = 10000)
plot(location.dist, ibs.dist, pch = 19,
     xlab = "geographic distance", ylab = "genetic distance")