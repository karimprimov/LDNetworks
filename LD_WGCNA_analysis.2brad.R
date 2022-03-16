require(WGCNA)
require(data.table)
require(flashClust)
require(vegan)
require(qqman)
require(cowplot)
require(tidyverse)

setwd("~/Desktop/Research/AdultJuv_Depth_FL_Set2/mcav.analysis/2brad.cluster/2brad.LD")

### CHOOSE SPECIES & LD PARAMETER ###
#~~~~~~~~~~~~
spp <- 'mcav' #ssid or mcav
ld.stat <- 'rEM' #rEM or DEM (r2 or D calculated with EM algorithm)
#~~~~~~~~~~~~

load('LDsquare_datt_traits.RData')
#if(ld.stat == 'rEM'){rm(ldmat.sq.DEM)} else if(ld.stat == 'DEM'){rm(ldmat.sq.rEM)}

#check that all sites are there
if(sum(!(colnames(datt) %in% names(mafs))) > 0) {
  datt <- datt[,-which(!(colnames(datt) %in% names(mafs)))]
}

# ---------------- pruning physically linked sites from LD table; retaining ones with higher minor allele freq

#remove sites where LD = Inf
prune.inf <- unique(rownames(which(ldmat.sq.rEM==Inf, arr.ind = T)))

mindist=1000

#folding mafs
#here "minor" is defined as the lower frequency allele (as opposed to the non-reference allele)
#if the non-reference allele frequency is greater than 0.5 in the mafs file, we subtract from 1 to get the frequency of the "minor" allele
mafs.fold=mafs
mafs.fold[mafs.fold>0.5]=1-mafs.fold[mafs.fold>0.5]

chrom.tm=factor(sub("\\:.+","",sites))
pos.tm=as.numeric(as.character(sub(".+\\:","",sites)))
sites.order <- sites[order(as.numeric(chrom.tm), pos.tm)]
chrom=factor(sub("\\:.+","",sites.order))
pos=as.numeric(as.character(sub(".+\\:","",sites.order)))

#Create data frame of loci ordered by position within each chromosome
sites.df <- data.frame(sites = sites.order, 
                       chrom = sub("\\:.+","",sites.order), 
                       pos = as.numeric(as.character(sub(".+\\:","",sites.order))),
                       maf = mafs[match(sites.order, names(mafs))]) #which mafs dataset you choose here is consequential (mafs vs. mafs_fold)
sites.df$sites <- as.character(sites.df$sites)

#For de novo 2bRAD genomes, linkage filtering must account for SNPs on the same tag but there is no way of knowing the genomic distance between each tag
if(spp == 'ssid'){
  #Generate successive intervals for each RAD tag in the de novo assembly
  tags.df <- data.frame(beg = seq(1,max(sites.df$pos), by=36), end = seq(36,max(sites.df$pos)+36, by=36))
  tags <- c()
  #Match each SNP to one of the RAD tags based on which interval falls within
  for(i in 1:nrow(sites.df)){tags <- c(tags, which(tags.df$beg <= sites.df$pos[i] & tags.df$end >= sites.df$pos[i]))}
  sites.df$tag <- as.factor(tags)
  #Group SNPs by chromosome and tag number and calculate the max maf for SNPs that fall within the same tag
  sites.df <- sites.df %>% group_by(chrom,tag) %>% mutate(max.maf = max(maf))
  #Identify any SNPs in each tag group that do NOT match the max maf. These will be pruned.
  prune.link <- colnames(datt)[which(sites.df$maf != sites.df$max.maf)]
} else if(spp == 'mcav'){
  #Add column of distance between each site
  sites.df <- sites.df %>% group_by(chrom) %>% mutate(diff = c(Inf, diff(pos)))
  links <- which(sites.df$diff < mindist)
  
  linkGrp <- function(x){
    if(!is.numeric(x)) x <- as.numeric(x)
    n <- length(x)
    y <- x[-1] != x[-n] + 1
    i <- c(which(y|is.na(y)),n)
    list(
      nsites = diff(c(0,i)), #number of sites in each linkage group
      begin = x[head(c(0,i)+1,-1)] - 1 #first index of each linkage group
    )
  }
  link.groups <- linkGrp(links)
  
  #Retain the site within each "linkage" group with the highest allele frequency
  #All others are filtered out of the dataset
  prune.link <- c()
  for(i in 1:length(link.groups$nsites)){
    prune.link <- prune.link
    link.grp <- seq(0,link.groups$nsites[i]) + link.groups$begin[i]
    sites.grp <- sites.df[link.grp,]
    max.maf <- which.max(sites.grp$maf)
    prune.link <- c(prune.link, as.character(sites.grp$sites[-max.maf]))
  }
}

#remove monomorphic sites (this should be NULL; if not, something went wrong)
prune.mono <- names(which(apply(datt, 2, var, na.rm=TRUE) == 0))

prune.all <- c(prune.link, prune.mono)
message(length(prune.all)) #2765 (mcav 2brad transect 2) (2778 in try 1)
if(ld.stat == 'DEM'){
  ldmat.filt <- ldmat.sq.DEM[-which(colnames(ldmat.sq.DEM) %in% prune.all),-which(colnames(ldmat.sq.DEM) %in% prune.all)]
} else if (ld.stat == 'rEM'){
  ldmat.filt <- ldmat.sq.rEM[-which(colnames(ldmat.sq.rEM) %in% prune.all),-which(colnames(ldmat.sq.rEM) %in% prune.all)]
} else if (ld.stat == 'rPrsn'){
  ldmat.filt <- ldmat.sq.rPrsn[-which(colnames(ldmat.sq.rPrsn) %in% prune.all),-which(colnames(ldmat.sq.rPrsn) %in% prune.all)]
}
mafs.fold=mafs.fold[-which(names(mafs.fold) %in% prune.all)]
mafs=mafs[-which(names(mafs) %in% prune.all)]
datt=datt[,-which(colnames(datt) %in% prune.all)]

rm(ldmat.sq.rEM)

ldmat.filt2 <- ldmat.filt
ldmat.filt2[ldmat.filt2<0.1] =0

#--------------- WGCNA

# Lowest power for which the the scale-free topology fit index reaches 0.90
# Co-expression similarity matrix is defined as the absolute value of the correlation coefficient between loci
# i.e., square root of R2 from ngsLD (power = 0.5)
power=1 # this is equivalent to a soft thresholding power of 2 since we are using R2 as input
ldmat.filt2[ldmat.filt2==Inf]=1
TOM = TOMsimilarity(ldmat.filt2^power, TOMType="unsigned");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average")

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
data.frame(table(dynamicColors))
# with min mod size of 30: there are 18 initial clusters (14 with try 1)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0,
                    addGuide = FALSE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes (first principal component of geno matrix subset by each gene module)
eg <- datt[,dynamicColors=='blue']
eg.pca <- rda(eg, scale = T)
head(eg.pca$CA$u[,1])

MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
MEs$MEgrey=NULL
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average")

#####
##### MEList$isHub #####
#####

MEDissThres = 0.3 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
#pdf("clusterModuleEigen.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot
#dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs


# plotting the fabulous ridiculogram (Figure 5A)
#pdf(paste0('set1/', spp, "/wgcna/wgcna_", ld.stat,"_geneModuleTree.pdf"), width = 3, height = 2)
plotDendroAndColors(geneTree,mergedColors,
                    dendroLabels = FALSE, hang = 0.0,
                    addGuide = F, guideHang = 0.05,lwd=0.6)
#dev.off()
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MEs$MEgrey=NULL

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
#	sizeGrWindow(7, 6)
plot(METree, main = "mc_all",
     xlab = "", sub = "")
#dev.off()

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs$MEgrey=NULL

# correlations of genes with eigengenes (same as signedkME output between MEs and datt, aka "module membership")
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

# ------------ Module-trait correlations

if(spp == 'mcav'){
  Traits=traits[,c("is1", "is2", "is4", "isA", "isN", "isO", "isD")]
  traitlabels <- c("C1", "C2", "C4", "Age", "Near", "Off", "Deep")
} else if(spp == 'ssid'){
  Traits=traits[,c("is1", "is2", "is4", "is3", "isA", "isN", "isO", "isD")]
  traitlabels <- c("Shallow1", "Shallow2", "Deep1", "Deep2", "Age", "Near", "Off", "Deep")
}
#Traits$isDA <- as.numeric(traits$site=="D" & traits$age=="A")
#Traits$isDJ <- as.numeric(traits$site=="D" & traits$age=="J")
moduleTraitCor = cor(MEs, Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
head(MEs)

#pdf(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_MEheatmap.pdf"), width = 4.5, height = 3.2)
textMatrix = paste0(signif(moduleTraitCor, 2), "\n",
                    ifelse(moduleTraitPvalue > 0.05, "NS", 
                           ifelse(moduleTraitPvalue > 0.01, "*",
                                  ifelse(moduleTraitPvalue > 0.001, "**", "***"))))
dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3));
# reorder table according to previous order
#	moduleTraitCor= moduleTraitCor[mtrows,]

# Display the correlation values within a heatmap plot (Figure 5B)
# Positive correlation indicates elevated frequency of derived (minor) alleles
# Negative correlation indicates elevated frequency of ancestral (major/reference) alleles
heatlabels <- paste0(Hmisc::capitalize(sub("ME","",row.names(moduleTraitCor))), '\n(n=', as.vector(table(moduleColors)[sub("ME","",row.names(moduleTraitCor))]), ')')
par(mar=c(2,2,2,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = traitlabels,
               xLabelsAngle = 70,
               #xLabelsAdj = c(0.5,0),
               #    yLabels = names(MEs),
               yLabels = heatlabels,
               ySymbols = row.names(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               textAdj = c(0.5, 0.6),
               cex.lab = 0.7,
               zlim = c(-1,1))
#dev.off()

#	mtrows=row.names(moduleTraitCor)
print(data.frame(table(moduleColors))) # number of genes in each module
# MCAV
# moduleColors Freq
#         blue 2059
#        brown  364
#    turquoise 3617
#       yellow  295
# SSID
# moduleColors Freq
#         blue 3891
#        brown  867
#        green  262
#    turquoise 5578

row.names(MEs)=row.names(traits)
colnames(MEs)=sub("ME","",colnames(MEs))
head(MEs)

geneModuleMembership = as.data.frame(signedKME(datt, MEs));
colnames(geneModuleMembership) = colnames(MEs)
save(datt,traits,Traits,MEs,METree,geneTree,moduleColors,moduleLabels,nSamples, nGenes, moduleGeneCor, moduleGenePvalue, moduleTraitPvalue, moduleTraitCor,geneModuleMembership,file=paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

# -------------- Are module SNPs clustered within contigs?

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

modCoord <- data.frame(matrix(unlist(strsplit(rownames(geneModuleMembership), ':')), ncol = 2, byrow = TRUE), moduleColors) %>%
  rename(contig = X1, pos = X2, module = moduleColors) %>%
  mutate(isBlu = ifelse(module == 'blue', 1, 0),
         isBro = ifelse(module == 'brown', 1, 0),
         isTur = ifelse(module == 'turquoise', 1, 0),
         isYel = ifelse(module == 'yellow', 1, 0))
obs <- modCoord %>%
  select(contig, isBlu, isBro, isTur, isYel) %>%
  group_by(contig) %>%
  summarise_all(sum)

# Chi-squared
chisq.byMod <- modCoord %>%
  group_by(contig) %>%
  # Expected proportion of SNPs in each contig (number of SNPs in each contig divided by total)
  summarise(prop = n()/6335) %>%
  # Observed and expected counts
  mutate(isBlu_o = obs$isBlu, isBlu_e = prop*sum(obs$isBlu),
         isBro_o = obs$isBro, isBro_e = prop*sum(obs$isBro),
         isTur_o = obs$isTur, isTur_e = prop*sum(obs$isTur),
         isYel_o = obs$isYel, isYel_e = prop*sum(obs$isYel)) %>%
  # Chi-squared statistic
  mutate(isBlu_chi = (isBlu_o-isBlu_e)^2/isBlu_e,
         isBro_chi = (isBro_o-isBro_e)^2/isBro_e,
         isTur_chi = (isTur_o-isTur_e)^2/isTur_e,
         isYel_chi = (isYel_o-isYel_e)^2/isYel_e) %>%
  select(isBlu_chi, isBro_chi, isTur_chi, isYel_chi) %>%
  colSums()

is.signif.chisq <- sapply(chisq.byMod, function(x) pchisq(x, df = 1376-1, lower.tail = FALSE))


#------------------ density plots of derived/ancestral alleles in modules

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

for(module in names(MEs)){
  pdf(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_alleleDensity_", module, ".pdf"), width = 7.5, height = 5.5)
  par(mfrow=c(3,3))
  kmeCutoff=0.25
  tt=traits[,c("is1","is2","is3","is4","isA","isD","isO","isN")]
  modNames = names(MEs)
  for(i in 1:length(tt)) {
    selTrait = as.data.frame(tt[,i]);
    whichTrait=names(tt)[i]
    geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
    moduleGenes = which(moduleColors==module);
    #	plot(density(geneTraitSignificance[moduleGenes,1]),main=module,xlim=c(-1,1))
    
    #retain top percentage of genes in specified module based on absolute value of module membership
    gm=abs(geneModuleMembership[moduleGenes, module])
    moduleGenesTop=moduleGenes[which(gm>quantile(gm,1-kmeCutoff))]
    
    #subset geno data by NON-module genes and samples associated with the trait
    moduledatt0=datt[as.logical(tt[,whichTrait]),which(moduleColors!=module)]
    
    #and calculate average number of derived alleles for each gene
    meangt0=apply(moduledatt0,2,function(x) mean(x, na.rm = T))
    d0=density(meangt0)
    d0$y=d0$y/max(d0$y)
    #    plot(d0,col="cyan3",main=paste(module,":",whichTrait),xlab="average number of derived alleles",mgp=c(2.3,1,0),bty="n",yaxt="n",ylab="")
    #	mtext(side=2,"Normalized density",cex=0.5)
    
    #subset geno data by module genes and samples NOT associated with the trait
    moduledatt1=datt[!as.logical(tt[,whichTrait]),moduleGenes]
    meangt1=apply(moduledatt1,2,function(x) mean(x, na.rm = T))
    d1=density(meangt1)
    d1$y=d1$y/max(d1$y)
    plot(d1,col="cyan3",main=paste(module,":",whichTrait),xlab="average number of derived alleles",mgp=c(2.3,1,0),bty="n",yaxt="n",ylab="", xlim = c(0,2))
    #    lines(d1,col="grey40", lty=3)
    
    #subset geno data by module genes and samples that ARE associated with the trait
    moduledatt2=datt[as.logical(tt[,whichTrait]),moduleGenes]
    meangt2=apply(moduledatt2,2,function(x) mean(x, na.rm = T))
    d2=density(meangt2)
    d2$y=d2$y/max(d2$y)
    lines(d2,col="goldenrod")
    
    # subset geno data by module genes with top 25% module membership and samples that ARE associated with the trait
    moduledatt3=datt[as.logical(tt[,whichTrait]),moduleGenesTop]
    meangt3=apply(moduledatt3,2,function(x) mean(x, na.rm = T))
    d3=density(meangt3)
    d3$y=d3$y/max(d3$y)
    lines(d3,col="red")
    
    #	abline(v=0,lty=3)
  }
  plot(d0,col="grey50",type="n",xlab="",mgp=c(2.3,1,0),bty="n",yaxt="n",xaxt="n",ylab="",main="")
  #  legend("topright",lty=c(1,1,3,1),col=c("red","goldenrod","grey40","cyan3"),lwd=2,legend=c("top 25% kME genes & trait","module genes & trait","module genes & non-trait","non-module genes & trait"),bty="n")
  legend("topright",lty=c(1,1,1),col=c("red","goldenrod","cyan3"),lwd=2,legend=c("top 25% kME genes & trait","module genes & trait","module genes & non-trait"),bty="n")
  
  dev.off()
}

# ---- Module, kme per gene (-log10 p-value)

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))
genes=read.table(paste0("set1/", spp, "/mcav_gene_regions.tab")) %>%
  rename(contig = V1, start = V2, end = V3, gene = V4) %>%
  filter(end > start) %>%
  mutate(contig = as.character(contig), gene = as.character(gene))
genes_2000=read.table(paste0("set1/", spp, "/mcav_genes.txt")) %>%
  rename(contig = V1, start = V2, end = V3, gene = V4) %>%
  filter(end-2000 > start+2000) %>%
  mutate(contig = as.character(contig), gene = as.character(gene))

modNames = names(MEs)

# Generate list of each module's SNPs with associated contig, position and -log10 p-value of module membership
coords=c()
for(module in modNames){
  moduleGenes = names(moduleGenePvalue[paste("ME",module,sep=""), moduleColors==module])
  coords[[module]]=plyr::ldply(strsplit(moduleGenes,":"))
  names(coords[[module]])=c("contig","pos")
  coords[[module]]$pos=as.numeric(coords[[module]]$pos)
  coords[[module]]$lpv=-log(moduleGenePvalue[paste("ME",module,sep=""),moduleColors==module],10)
  coords[[module]]$lpv=as.numeric(coords[[module]]$lpv)
}

# Generate data frame of module membership (-log10 p-value) for functional genes (max membership if multiple SNPs in gene)
ngenes=nrow(genes)
pb=txtProgressBar(0,ngenes)
Gmax=c()
for (i in 1:ngenes) {
  setTxtProgressBar(pb,i)
  gmax=c()
  for(module in modNames){
    s=subset(coords[[module]],contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
    if (nrow(s)>0) { gmax=c(gmax,max(s$lpv)) } else { gmax=c(gmax,0)} #if there are multiple SNPs within a gene, keep the max log10 p-value of module membership
  }
  Gmax=rbind(Gmax,gmax)
}
Gmax=data.frame(Gmax)
names(Gmax)=modNames
row.names(Gmax)=genes$gene
table(Gmax[,1]>0) # 1382 blue
table(Gmax[,2]>0) # 257 brown
table(Gmax[,3]>0) # 2396 turquoise
table(Gmax[,4]>0) # 189 yellow
save(Gmax,file=paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_modmemByGene.RData"))
head(Gmax)
dim(Gmax)
for(module in modNames){
  togo=data.frame(cbind("gene"=as.character(genes$gene),"moduleLPV"=Gmax[,module]))
  write.csv(togo,row.names=F,quote=F,file=paste0('set1/', spp,"/wgcna/wgcna_",ld.stat,"_",module,".csv"))
}

#----------------- saving data for GO analysis (gene x module kME)

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

# bads=which(is.na(apply(geneModuleMembership,1,mean)))
# geneModuleMembership= geneModuleMembership[-bads,]
# moduleColors= moduleColors[-bads]

#names(geneModuleMembership)=sub("MM","", names(geneModuleMembership))

coords=row.names(geneModuleMembership)
coo=c();i=1
for (i in 1:length(coords)) {
  co=strsplit(coords[i],":")
  coo=rbind(coo,c(co[[1]][1],co[[1]][2]))
}

coo=data.frame(coo)
names(coo)=c("contig","pos")
coo$contig=as.character(coo$contig)
coo$pos=as.numeric(as.character(coo$pos))
str(coo)

# Generate module-specific data frames reporting only kME for each SNP in each corresponding module
# SNPs !in each module are reported as 0
blue=data.frame(cbind(coo,"kme"=abs(geneModuleMembership[,"blue"])))
blue[moduleColors!="blue",3]=0
brown=data.frame(cbind(coo,"kme"=abs(geneModuleMembership[,"brown"])))
brown[moduleColors!="brown",3]=0
turquoise =data.frame(cbind(coo,"kme"=abs(geneModuleMembership[,"turquoise"])))
turquoise[moduleColors!="turquoise",3]=0
yellow =data.frame(cbind(coo,"kme"=abs(geneModuleMembership[,"yellow"])))
yellow[moduleColors!="yellow",3]=0

head(yellow)

plot(density(blue[,3]),ylim=c(0,30))
plot(density(turquoise[,3]),ylim=c(0,30))
plot(density(yellow[,3]),ylim=c(0,30))
plot(density(brown[,3]),ylim=c(0,30))
table(moduleColors)
# blue     brown turquoise    yellow 
# 2059       364      3617       295 

write.table(blue,sep="\t",quote=F,file=paste0('set1/', spp,"/wgcna/wgcna_",ld.stat,"_blueKME.tab"))
write.table(brown,sep="\t",quote=F,file=paste0('set1/', spp,"/wgcna/wgcna_",ld.stat,"_brownKME.tab"))
write.table(turquoise,sep="\t",quote=F,file=paste0('set1/', spp,"/wgcna/wgcna_",ld.stat,"_turquoiseKME.tab"))
write.table(yellow,sep="\t",quote=F,file=paste0('set1/', spp,"/wgcna/wgcna_",ld.stat,"_yellowKME.tab"))

# kme per gene
genes=read.table(paste0("set1/", spp, "/mcav_gene_regions.tab")) %>%
  rename(contig = V1, start = V2, end = V3, gene = V4) %>%
  filter(end > start) %>%
  mutate(contig = as.character(contig), gene = as.character(gene))
genes_2000=read.table(paste0("set1/", spp, "/mcav_genes.txt")) %>%
  rename(contig = V1, start = V2, end = V3, gene = V4) %>%
  filter(end-2000 > start+2000) %>%
  mutate(contig = as.character(contig), gene = as.character(gene))

# How many SNPs are within annotated gene boundaries?
coo.ingene <- data.frame(coo, in_gene = NA, gene = NA)
for(i in 1:nrow(coo)){
  sub <- subset(genes, contig == coo$contig[i] & start <= coo$pos[i] & end >= coo$pos[i])
  coo.ingene$in_gene[i] <- ifelse(nrow(sub)==0, 'n', 'y')
  coo.ingene$gene[i] <- ifelse(nrow(sub)==0, NA, sub$gene)
}
n.snps.ingene <- sum(coo.ingene$in_gene == 'y')
perc.snps.ingene <- n.snps.ingene/nrow(coo)*100
message(paste0(n.snps.ingene, ' SNPs (', round(perc.snps.ingene, 2), '%) are within annotated gene boundaries.'))

# For each module find SNPs within bounds of annotated genes
# If there is more than one SNP within a gene, report the max kME for that module
# In many cases, SNPs associated with different modules occur within the same gene so kME of each module is reported
kme.list <- list()
for (module in gsub('MM','',colnames(geneModuleMembership))){
  dat=get(module)
  i=1;kme=c();ns01=0
  print(module)
  pb=txtProgressBar(0,nrow(genes))
  for (i in 1:nrow(genes)) {
    setTxtProgressBar(pb,i)
    sub=subset(dat,contig==genes$contig[i] & pos>=genes$start[i] & pos<=genes$end[i])
    if (nrow(sub)==0) { 
      kme=append(kme,NA)
    } else {
      kme=append(kme,max(sub$kme))
      ns01=ns01+nrow(sub)
    }
  }
  table(is.na(kme))
  plot(density(na.omit(kme)))
  
  res=data.frame(cbind("gene"=as.character(genes$gene),kme))
  res=na.omit(res)
  res$kme=as.numeric(as.character(res$kme))
  kme.list[[module]]=res
}

bygene=data.frame(cbind("blue"=kme.list$blue$kme,"yellow"=kme.list$yellow$kme,"turquoise"= kme.list$turquoise$kme,"brown"=kme.list$brown$kme))
row.names(bygene)=kme.list$blue$gene

pairs(bygene,pch=16,col=rgb(0,0,0,alpha=0.5),cex=0.7)

save(bygene,file=paste0("set1/",spp,"/wgcna/wgcna_",ld.stat,"_bygene_kme.RData"))

write.csv(kme.list$turquoise,file=paste0("set1/",spp,"/wgcna/wgcna_",ld.stat,"_turquoise.bygene.csv"),row.names=F, quote=F)
write.csv(kme.list$blue,file=paste0("set1/",spp,"/wgcna/wgcna_",ld.stat,"_blue.bygene.csv"),row.names=F, quote=F)
write.csv(kme.list$yellow,file=paste0("set1/",spp,"/wgcna/wgcna_",ld.stat,"_yellow.bygene.csv"),row.names=F, quote=F)
write.csv(kme.list$brown,file=paste0("set1/",spp,"/wgcna/wgcna_",ld.stat,"_brown.bygene.csv"),row.names=F, quote=F)

#------------ Module membership vs. Fst (Figure 5C)

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

baye.fst = read.table(paste0("set1/", spp, "/mc_0618.baye_fst_pos.txt"), header = T)
baye.fst$pos = paste(baye.fst$CHROM, baye.fst$POS, sep = ":")
baye.fst$out = baye.fst$qval<0.05
geneMod = geneModuleMembership %>% mutate(pos = rownames(geneModuleMembership))
geneModBaye = left_join(geneMod, baye.fst, by = 'pos') %>% select(-c("CHROM","POS"))

# How many Bayescan outliers fall within annotated gene boundaries?
baye.outs <- data.frame(subset(baye.fst, out == T), in_gene = NA, gene = NA)
for(i in 1:nrow(baye.outs)){
  sub <- subset(genes, contig == baye.outs$CONTIG[i] & start <= baye.outs$POS[i] & end >= baye.outs$POS[i])
  baye.outs$in_gene[i] <- ifelse(nrow(sub)==0, 'n', 'y')
  baye.outs$gene[i] <- ifelse(nrow(sub)==0, NA, sub$gene)
}
n.outs.ingene <- sum(baye.outs$in_gene == 'y')
perc.outs.ingene <- n.outs.ingene/nrow(baye.outs)*100
message(paste0(n.outs.ingene, ' outlier SNPs (', round(perc.outs.ingene, 2), '%) are within annotated gene boundaries.'))

# Identify the maximum gene-module correlation (absolute value) across the four modules
geneModBaye$color = apply(geneModBaye[,colnames(geneModuleMembership)], 1, function(x) gsub("MM", "", colnames(geneModuleMembership)[which.max(abs(x))]))
#geneModBaye$color[geneModBaye$color == 'brown'] <- 'tan3'
geneModBaye$alpha = ifelse(geneModBaye$out, 1, 0.2)
geneModBaye$size = ifelse(geneModBaye$out, 1.1, 1)
geneModBaye$stroke = ifelse(geneModBaye$out, 0.5, 0)

geneModBaye$abs.max = apply(geneModBaye[,colnames(geneModuleMembership)], 1, function(x) max(abs(x)))
geneModBaye$true.max = apply(geneModBaye[,colnames(geneModuleMembership)], 1, function(x) x[which.max(abs(x))])

#pdf(paste0('set1/', spp, "/wgcna/wgcna_", ld.stat,"_geneModuleFst.pdf"), width = 5.5, height = 3.5)
ggplot(geneModBaye, aes(x=true.max, y=fst, fill = color, size = size, alpha = alpha, stroke = stroke))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle("Gene-Module Membership vs. Fst")+
  xlab("Maximum Gene-Module Correlation")+
  ylab("Alpha (Locus-specific component of Fst)")+
  geom_vline(xintercept = 0, linetype = 2)+
  geom_point(shape = 21)+
  scale_y_continuous(limits = c(0,0.6), expand = c(0,0))+
  scale_alpha_identity(guide = "none")+
  scale_size_identity(guide = "none")+
  scale_fill_identity(guide = "legend")
dev.off()

summary(lm(fst~abs.max, geneModBaye))

mod2plot <- 'yellow'
geneModBaye.byMod <- geneModBaye[moduleColors == mod2plot,] %>%
  select(pos, fst, blue, brown, turquoise, yellow) %>%
  gather('module', 'kME', 3:6)
ggplot(geneModBaye.byMod, aes(x = kME, y = fst, group = module))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())+
  geom_vline(xintercept = 0, linetype = 'dashed')+
  geom_point(aes(fill = ifelse(module == mod2plot, mod2plot, 'grey40'),
                 size = ifelse(module == mod2plot, 1.1, 1),
                 alpha = ifelse(module == mod2plot, 0, 0.2)))+
  geom_point(aes(fill = ifelse(module == mod2plot, mod2plot, 'grey40'),
                 size = ifelse(module == mod2plot, 1.1, 1),
                 stroke = ifelse(module == mod2plot, 0.5, 0),
                 alpha = ifelse(module == mod2plot, 1, 0)), shape = 21)+
  scale_size_identity(guide = "none")+
  scale_alpha_identity(guide = "none")+
  scale_fill_identity()+
  scale_x_continuous(limits = c(-1,1))+
  scale_y_continuous(limits = c(0,0.6), expand = c(0,0))

# Pairwise per-locus Fst
fst.files <- grep(list.files(paste0('set1/', spp, "/wgcna/")), pattern = glob2rx('p*.fst'), value = T)
pop.fsts <- list()
for (i in fst.files) {
  pop.fsts[[gsub('\\.fst', '', i)]] <- read.table(paste0('set1/', spp, "/wgcna/", i))
}
pop.fsts.df <- do.call(rbind, pop.fsts) %>%
  mutate(pos = paste(V1, V2, sep = ':'), V3 = ifelse(V3 < 0, 0, V3), pair = gsub('\\..+', '', rownames(.))) %>%
  mutate(fst = V3/V4) %>%
  select(pos, fst, pair) %>%
  right_join(data.frame(pos = geneMod$pos, module = moduleColors), by = 'pos') %>%
  group_by(pos, module) %>%
  summarise(pop1.fst = mean(fst[grepl('1', pair)], na.rm = T),
            pop2.fst = mean(fst[grepl('2', pair)], na.rm = T),
            pop3.fst = mean(fst[grepl('3', pair)], na.rm = T),
            pop4.fst = mean(fst[grepl('4', pair)], na.rm = T)) %>%
  left_join(geneMod, by = c('pos')) %>%
  rename(topMod = 'module') %>%
  gather('module', 'kME', 7:10)

pop.fsts.plot <- function(df, mod2plot, pop2plot, kmeCutoff){
  if(missing(pop2plot)){
    pop2plot <- unname(which.max(moduleTraitCor[paste0('ME', mod2plot), 1:4]))
  }
  df2plot <- pop.fsts.df[pop.fsts.df$module == mod2plot,]
  
  if(missing(kmeCutoff)){
    kmeCutoff <- 1
  }
  moduleGenes <- which(moduleColors==mod2plot)
  memStrength <- abs(geneModuleMembership[moduleGenes, mod2plot])
  moduleGenesTop <- moduleGenes[which(memStrength > quantile(memStrength, 1-kmeCutoff))]
  
  df2plot.quant <- df2plot[match(rownames(geneModuleMembership), df2plot$pos),] %>% mutate(quant = FALSE)
  df2plot.quant$quant[moduleGenesTop] <- TRUE
  
  scatter <- ggplot(df2plot.quant, aes_string(x = 'kME', y = paste0('pop', as.character(pop2plot), '.fst')))+
    theme_bw()+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'cm'))+
    geom_vline(xintercept = 0, linetype = 'dashed')+
    geom_point(aes(fill = ifelse(quant, mod2plot, 'grey40'),
                   size = ifelse(quant, 1.1, 1),
                   alpha = ifelse(quant, 0, 0.2)))+
    geom_point(aes(fill = ifelse(quant, mod2plot, 'grey40'),
                   size = ifelse(quant, 1.1, 1),
                   stroke = ifelse(quant, 0.5, 0),
                   alpha = ifelse(quant, 1, 0)), shape = 21)+
    scale_size_identity(guide = "none")+
    scale_alpha_identity(guide = "none")+
    scale_fill_identity()+
    scale_x_continuous(limits = c(-1.1,1.1), expand = c(0,0))+
    scale_y_continuous(limits = c(0,1.1), expand = c(0,0))
  
  x.plot <- ggplot(df2plot.quant, aes_string('kME')) + 
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          plot.margin = unit(c(0,0,0,32), 'pt')) +
    geom_density(aes(fill = ifelse(quant, mod2plot, 'grey40')), alpha=.5) + 
    scale_fill_identity() +
    scale_x_continuous(limits = c(-1.1,1.1), expand = c(0,0))
  
  y.plot <- ggplot(df2plot.quant, aes_string(paste0('pop', as.character(pop2plot), '.fst'))) + 
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          plot.margin = unit(c(0,0,24,0), 'pt')) +
    geom_density(aes(fill = ifelse(quant, mod2plot, 'grey40')), alpha=.5, color = 'black') + 
    scale_fill_identity() +
    scale_x_continuous(limits = c(0,1.1), expand = c(0,0)) +
    coord_flip()
  
  # Combine all plots together and flatten density plots with rel_heights
  combine.plot <- plot_grid(x.plot, NULL, scatter, y.plot, ncol = 2, rel_widths = c(3, 1), rel_heights = c(1, 3))
  #  print(combine.plot)
  save_plot(paste0('set1/', spp, "/wgcna/popFsts/wgcna_",ld.stat,"_popFsts_",mod2plot,"_q",kmeCutoff,"_pop",pop2plot,".pdf"), combine.plot, base_width = 3, base_height = 2)
}

pop.fsts.plot(pop.fsts.df, 'turquoise')

for(mod in c('blue', 'brown', 'yellow', 'turquoise')){
  for(pop in 1:4){
    pop.fsts.plot(pop.fsts.df, mod2plot = mod, pop2plot = pop)
  }
}


# Heterozygosity ----------------------------------------------------------

# Save lists of module genes
load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))
for(mod in unique(moduleColors)){
  mod.sites <- rownames(geneModuleMembership)[moduleColors == mod]
  mod.sites.df <- matrix(unlist(strsplit(mod.sites, ':')), ncol = 2, byrow = T)
  write.table(mod.sites.df, file = paste0('set1/', spp, '/wgcna/', mod, '.sites'), sep = '\t', row.names = F, col.names = F, quote = F)
}

# Once sample-specific heterozygosity is calculated with ANGSD/realSFS... (see README)
# Join heterozygosity values with traits df
het <- read.table(paste0("set1/", spp, "/wgcna/het_byIndModule_", spp), header = T) %>%
  full_join(data.frame(ind = rownames(traits), k = traits$k), by = 'ind') %>%
  mutate(pop = as.character(pop), k = as.character(k))
het.byPop <- filter(het, gsub("c", "", pop) == k)

ggplot()+
  geom_boxplot(data = het.byPop, aes(x = mod, y = het, fill = mod))+
  scale_fill_manual(values = c('blue', 'brown', 'turquoise', 'yellow'))

het.byMod <- list()
for(modcol in c('turq', 'blue', 'brown', 'yellow')){
  het.byMod[[modcol]] <- filter(het, mod == modcol) %>%
    mutate(inPop = ifelse(gsub("c", "", pop) == k, 'in', 'out')) %>%
    filter(!is.na(k))
}
het.byMod <- do.call(rbind, het.byMod)
het.byMod$plotCol <-  ifelse(het.byMod$inPop == 'in', as.character(het.byMod$mod), 'grey70')
het.byMod$plotCol[het.byMod$plotCol == 'turq'] <- 'turquoise'

levels(het.byMod$mod) <- c('Blue','Brown','Turquoise','Yellow')
heteroz <- ggplot(data = het.byMod, aes(x = inPop, y = het, fill = plotCol))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank())+
  ylab(expression("Heterozygosity" ~ (H[O])))+
  geom_violin(draw_quantiles = c(0.5), size = 0.2)+
  scale_fill_identity()+
  scale_y_continuous(limits = c(0,0.75))+
  facet_wrap(~mod, nrow = 1)
heteroz
save_plot(paste0('set1/', spp, '/wgcna/heteroz_', spp, '.pdf'), heteroz, base_width = 4, base_height = 3)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# EXTRA ANALYSES ----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Discriminant analysis of individuals based on module eigengenes
rr=rda(MEs,scale=F)
#biplot(rr, col=traits$k)
plot(rr$CA$u,pch=as.numeric(traits$k)-1, col="grey50", xlim=c(-0.25,0.2),ylim=c(-0.25,0.4))
arrowScale=6
mod.cols <- gsub("ME","",colnames(MEs))
nmods=length(mod.cols)
arrows(rep(0,nmods),rep(0,nmods),rr$CA$v[c(1:nmods),1]/arrowScale,rr$CA$v[c(1:nmods),2]/arrowScale,length=0.1,lwd=2,col=mod.cols)
legend("bottomleft",pch=unique(as.numeric(traits$k)-1),legend=unique(traits$k),bty="n",title="cluster",cex=0.9)
dev.off()

conds2=cbind(traits,MEs)
conds2$k=factor(conds2$k)

pdf(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_eigPopCluster.pdf"), width = 8, height = 6)
par(mfrow=c(2,3))
for(i in 1:nmods){plot(get(mod.cols[i])~k, conds2, main=mod.cols[i], mgp=c(2.1,1,0), ylab = "Eigengene")}
dev.off()

# Discriminant analysis of SNPs based on module membership/correlation
par(mfrow=c(1,1))
rr=rda(t(moduleGeneCor))
plot(rr$CA$eig)
arrowScale=20
axes2plot=c(1,2)
plot(rr$CA$u[,axes2plot], pch=16, cex=0.5, col=rgb(0,0,0,alpha=0.1), mgp=c(2.3,1,0), ylim = c(-0.05, 0.05))
arrows(rep(0,nmods), rep(0,nmods), rr$CA$v[c(1:nmods), axes2plot[1]]/arrowScale, rr$CA$v[c(1:nmods), axes2plot[2]]/arrowScale, length=0.1, lwd=2, col=mod.cols)
dev.off()


#------------ Manhattan plot

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

scaff.num <- as.factor(sub("\\:.+", "", rownames(geneModuleMembership)))
levels(scaff.num) <- seq(1:length(levels(scaff.num)))
scaff.num <- as.numeric(scaff.num)
kmeCutoff <- 0.25

trace(manhattan, edit = T)
# Edit line 112 to from col = "green3" to col = module2plot

topModMem.all <- list()
for(module2plot in names(MEs)){
  # Identify SNPs in top 25% percentile wrt SNP-module correlation ('geneModuleMembership')
  gm=abs(geneModuleMembership[which(moduleColors==module2plot), module2plot])
  moduleGenesTop=which(moduleColors==module2plot)[which(gm>quantile(gm,1-kmeCutoff))]
  # Plot -log10(p-value) for all SNPs, highlighting those identified above
  # In some cases, the highest SNP-module correlation does not reflect the highest -log10(p-value)
  mandf <- data.frame(SNP = seq(1, nGenes),
                      CHR = scaff.num,
                      BP = as.numeric(sub(".+\\:", "", rownames(geneModuleMembership))),
                      P = moduleGenePvalue[paste0("ME",module2plot),])
  pdf(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_manhattan_", module2plot, ".pdf"), width = 12, height = 4)
  manhattan(mandf, main = toupper(module2plot), highlight = moduleGenesTop, suggestiveline = FALSE, genomewideline = FALSE)
  dev.off()
}


#------------ By-trait scatters

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))
library(WGCNA)
library(plotrix) # for color.scale
length(moduleColors)
data.frame(table(moduleColors))

nGenes = ncol(datt);
nSamples = nrow(datt);
# names (colors) of the modules
modNames = names(MEs)
geneModuleMembership = as.data.frame(signedKME(datt, MEs)); # same as moduleGeneCor
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)); # same as moduleGenePvalue
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

for(whichTrait in c("is1","is2","is3","is4","isA","isD","isO","isN")){
  selTrait = as.data.frame(traits[,whichTrait]);
  names(selTrait) = whichTrait
  
  # correlation of genes to traits
  # (positive value suggests membership to trait category is correlated with more derived alleles at that site)
  geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
  names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
  
  pdf(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_geneSigModuleCorr_", whichTrait, ".pdf"), width = 7.5, height = 5.5)
  par(mfrow=c(2,3))
  counter=0
  for(module in modNames){
    counter=counter+1
    if (counter>9) {
      quartz()
      par(mfrow=c(3,3))
      counter=1
    }
    column = match(module, modNames);
    moduleGenes = which(moduleColors==module );
    moduleMafs = mafs[moduleColors==module]; #use mafs here instead of mafs.fold to represent frequency of derived allele
    moduleGenes=moduleGenes[order(moduleMafs)]
    length(moduleGenes)
    moduleMafs=moduleMafs[order(moduleMafs)]
    
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       ylim=c(-1,1),xlim=c(0,1),
                       geneTraitSignificance[moduleGenes, 1],
                       #	abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste(module,"kME"),
                       ylab = paste("GS for", whichTrait),
                       #	col = rgb(0,0,0,alpha=0.1),
                       col=color.scale(moduleMafs,c(0,1),c(0.8,0),c(1,0),alpha=0.5),
                       pch=16,mgp=c(2.3,1,0)
    )
    abline(0,1,lty=3,col="grey50")
    abline(0,-1,lty=3,col="grey50")
  }
  
  # plotting color scale bar
  colscale=data.frame(cbind("color"=1,"maf"=seq(0,1,0.01)))
  plot(color~maf,colscale,pch=16,cex=3,col=color.scale(seq(0,1,0.01),c(0,1),c(0.8,0),c(1,0)))
  
  dev.off()
}

# Plot interpretation:
# isD trait:
# For loci in the yellow module (n = 295) with a HIGH minor allele frequency across all samples (more red)...
# more "yellow-ness" = stronger POSITIVE correlation between # of derived alleles and being a deep individual
# For loci in the yellow module (n = 295) with a LOW minor allele frequency across all samples (more blue)...
# more "yellow-ness" = stronger NEGATIVE correlation between # of derived alleles and being a deep individual

# For ancestral loci (LOW minor allele frequency) in the yellow module, deep samples are more ancestral
# For derived loci (HIGH minor allele frequency) in the yellow module, deep samples are more derived


#------------ Genotype frequency (if datt includes hard calls)

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

plot_list <- list()
for(module in names(MEs)){
  kmeCutoff=0.25
  tt=traits[,c("is1","is2","is3","is4")]
  modNames = names(MEs)
  
  selTrait = as.data.frame(tt[,1])
  firstTrait=names(tt)[1]
  geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
  moduleGenes = which(moduleColors==module)
  
  #retain top percentage of genes in specified module based on absolute value of module membership
  gm=abs(geneModuleMembership[moduleGenes, module])
  moduleGenesTop=moduleGenes[which(gm>quantile(gm,1-kmeCutoff))]
  
  #subset geno data by module genes and samples associated with each trait
  modulegeno1=data.frame(table(datt_hard[as.logical(tt[,firstTrait]),moduleGenes]), pop=rep(firstTrait,3))
  for(each in names(tt)[names(tt) != firstTrait]){
    modulegeno1=rbind(modulegeno1,data.frame(table(datt_hard[as.logical(tt[,each]),moduleGenes]), pop=rep(each,3)))
  }
  modulegeno1 <- modulegeno1 %>% 
    group_by(pop) %>%
    mutate(Freq = Freq/sum(Freq),)
  
  # subset geno data by module genes with top 25% module membership and samples associated each the trait
  modulegeno2=data.frame(table(datt_hard[as.logical(tt[,firstTrait]),moduleGenesTop]), pop=rep(firstTrait,3))
  for(each in names(tt)[names(tt) != firstTrait]){
    modulegeno2=rbind(modulegeno2,data.frame(table(datt_hard[as.logical(tt[,each]),moduleGenesTop]), pop=rep(each,3)))
  }
  modulegeno2 <- modulegeno2 %>% 
    group_by(pop) %>%
    mutate(Freq = Freq/sum(Freq))
  
  topCorr <- names(which.max(moduleTraitCor[module,names(tt)]))
  
  plot_list[[module]] <- ggplot(data = modulegeno1)+
    theme_bw()+
    theme(panel.grid = element_blank(), legend.position = "none")+
    ggtitle(paste0(toupper(module), " Module"), subtitle = paste0(as.character(length(moduleGenes))," SNPs"))+
    labs(x="Genotype", y="Frequency") +
    geom_col(aes(x=factor(modulegeno1$Var1), y=Freq, fill=pop), position = "dodge", col = 'black', size = 0.2)+
    scale_fill_manual(values = ifelse(names(tt)==topCorr,module,"grey90"))+
    scale_x_discrete(labels=c("0" = "Homozygote\n(Ancestral)", 
                              "1" = "Heterozygote", 
                              "2" = "Homozygote\n(Derived)"))+
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  plot_list[[paste0(module,"top25")]] <- ggplot(data = modulegeno2)+
    theme_bw()+
    theme(panel.grid = element_blank(), legend.position = c(0.87,0.77), legend.key.size = unit(0.2, "in"))+
    ggtitle(paste0(toupper(module), " Module - top 25% kME"), subtitle = paste0(as.character(length(moduleGenesTop))," SNPs"))+
    labs(x="Genotype", y="Frequency") +
    geom_col(aes(x=factor(modulegeno1$Var1), y=Freq, fill=pop), position = "dodge", col = 'black', size = 0.2)+
    scale_fill_manual(values = ifelse(names(tt)==topCorr,module,"grey90"))+
    scale_x_discrete(labels=c("0" = "Homozygote\n(Ancestral)", 
                              "1" = "Heterozygote", 
                              "2" = "Homozygote\n(Derived)"))+
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
}

cp <- plot_grid(plotlist=plot_list, ncol = 2)
save_plot(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_genoFrequency.pdf"), cp, base_width = 8.5, base_height = 16)
dev.off()


#------------ Heatmap of LD values for sample of loci (n=1000)

library(viridis)

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

LD.long <- read.table(paste0("set1/", spp, "/wgcna/test_sample_Prsn.txt"), header = T)[,c(1,2,4)]
LD.long.inv <- data.frame(V1 = LD.long$V2, V2 = LD.long$V1, V4 = LD.long$V4)
LD.long.sq <- rbind(LD.long, LD.long.inv)
LD.long.sq <- LD.long.sq[LD.long.sq$V4 <= 0.2,]
LD.long.sq$V4 <- log10(LD.long.sq$V4 + 1)

ggplot(LD.long.sq, aes(V1, V2)) + 
  geom_tile(aes(fill = V4)) + 
  scale_fill_viridis(option = 'inferno', name = 'R2') +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

loc.mean <- read.table(paste0("set1/", spp, "/wgcna/meanLD_byLocus_Prsn.txt"), header = T)
colnames(loc.mean) <- c('pos','LD')
geneModBayeLD <- merge(geneModBaye, loc.mean)
geneModBayeLD <- geneModBayeLD[is.finite(geneModBayeLD$LD),]
ggplot(data = geneModBayeLD, aes(x=fst, y=LD))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "lm")

ggplot(geneModBayeLD, aes(x=fst, y=LD))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle("LD vs. Fst")+
  xlab("Fst")+
  ylab("LD (R2-EM)")+
  geom_point(shape = 21, aes(fill = color), stroke = 0.1, size = 1, alpha = 0.5)+
  geom_smooth(method = 'lm', color = 'black', lwd = 0.5)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_identity(guide = "legend")+
  facet_wrap(.~color)


#------------ Why does TURQUOISE module show opposite pattern from all others in SSIDs?

load(paste0('set1/', spp, "/wgcna/wgcna_",ld.stat,"_p1.RData"))

geneMod$module <- moduleColors
ggplot()+
  theme_bw()+
  ggtitle('Membership to TURQUOISE module')+
  xlab('Module SNPs')+
  ylab('kME')+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_boxplot(data = geneMod, aes(x = module, y = turquoise, fill = module))+
  scale_fill_identity()
a <- geneMod %>% 
  filter(turquoise > 0.5) %>%
  mutate(topMod = ifelse(module == 'turquoise', turquoise, 
                         ifelse(module == 'blue', blue, 
                                ifelse(module == 'yellow', yellow, brown))))
b <- geneMod %>% 
  filter(turquoise < -0.5) %>%
  mutate(topMod = ifelse(module == 'turquoise', turquoise, 
                         ifelse(module == 'blue', blue, 
                                ifelse(module == 'yellow', yellow, brown))))
ggplot(a)+
  theme_bw()+
  ggtitle('TURQUOISE kME > 0.5')+
  xlab('kME to Assigned Module')+
  ylab('# Sites')+
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = 0.5, linetype = 'dotted')+
  geom_histogram(aes(x=topMod, fill = module), binwidth = 0.04)+
  scale_fill_identity(guide = 'legend')+
  scale_x_continuous(limits = c(-1.1,1.1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,35), expand = c(0,0))+
  guides(fill=guide_legend(title="Assigned module"))