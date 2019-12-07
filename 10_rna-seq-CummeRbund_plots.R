if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cummeRbund")

library(cummeRbund)

setwd("output/10_cuffdiff/")

wd <- getwd()

list.files(wd)

cuff<-readCufflinks()
cuff


disp<-dispersionPlot(genes(cuff))

png("./dispersion.png",   # cria arquivo do tipo PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels de largura
    height = 5*300,       # 5 x 300 pixels de altura
    res = 300,            # 300 pixels por polegada
    pointsize = 8)        # tamanho da fonte

disp

dev.off()


png("./genes.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv
dev.off()

png("./dens.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
dens<-csDensity(genes(cuff))
dens
dev.off()

png("./densRep.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
densRep<-csDensity(genes(cuff),replicates=T)
densRep
dev.off()

png("./b.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
b<-csBoxplot(genes(cuff))
b
dev.off()

png("./s.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
s<-csScatterMatrix(genes(cuff))
s
dev.off()

png("./v.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
v<-csVolcano(genes(cuff),"SAMPLEA","SAMPLEB")
v
dev.off()

data(sampleData)

myGeneIds <- read.table("../../ref/sampleIDs")
myGeneIds <- as.character(myGeneIds$V1)

myGenes<-getGenes(cuff,myGeneIds)
myGenes

png("./h.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
h<-csHeatmap(myGenes,cluster='both')
h
dev.off()

png("./h.rep.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
h.rep<-csHeatmap(myGenes,cluster='both',replicates=T)
h.rep
dev.off()

png("./b2.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
b<-expressionBarplot(myGenes)
b
dev.off()

png("./gl.rep.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
gl.rep<-expressionPlot(myGenes,replicates=TRUE)
gl.rep
dev.off()

png("./pca.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.PCA.rep
dev.off()

png("./genes.MDS.rep.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)
genes.MDS.rep
dev.off()


ic<-csCluster(myGenes,k=4)
head(ic$cluster)

png("./icp.png",width = 5*300,height = 5*300,res = 300,pointsize = 8)
icp<-csClusterPlot(ic)
icp
dev.off()
