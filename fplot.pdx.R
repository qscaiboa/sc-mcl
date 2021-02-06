---
title: "Untitled"
output: html_document
---
#conda activate sc-r351-s234



```{r}
install.packages('gridBase')
install.packages('gridExtra')
install.packages('ggplot2')
install.packages('igraph')
install.packages('packcircles')
install_github('hdng/trees')
```

```{r}
rm(list=ls())
setwd("C:/Users/qcai1/Downloads/sciclone")
library(sciClone)
library(clonevol)
library(fishplot)
clone.colors <- c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e', "#f0aa08", "#faebd7", "#b8ead2")
```

# VAF reading
```{r}
library(qmrparser)
readvaf <- function(fileloc){
    sample <- read.table(fileloc, fill = T) #, stringsAsFactors=FALSE)

#sample <- sample[sample$V7 %in% "",]
sample <- sample[substring(sample$V1,1,3) %in% 'chr',]
sample <- sample[, c(1,2,5,6)]
colnames(sample) <- c("V1", "V2", "V3", "V4")
sample[,3:4] <- data.frame(lapply(sample[,3:4], as.character), stringsAsFactors=FALSE)

sample$V6 <- NA
i <- 1
while (i < nrow(sample)){
 sample[i,6]  <- isLetter(sample[i,3])
  i <- i+1
}
sample <- sample[!sample$V6.1,]

sample<- sample[substring(sample$V1,1,3) %in% 'chr',]
sample$V2 <- as.numeric(sample$V2)
sample$V3 <- as.numeric(sample$V3)
sample$V4 <- as.numeric(sample$V4)
sample$V5 =sample$V4 / (sample$V3+sample$V4)  *100
sample[,c(1:4,7)]
}
```



```{r}
PDX3 <- readvaf("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/Z2.PDX3.fishplot")
PDX56 <- readvaf("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/Z2.PDX56.fishplot")
#MCL726s<- readvaf("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/Z2.MCL726s.fishplot")
MCL726s1<- readvaf("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/Z2.MCL726s1.fishplot")
MCL726s2<- readvaf("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/Z2.MCL726s2.fishplot")
MCL726s3<- readvaf("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/Z2.MCL726s3.fishplot")
MCL726s4<- readvaf("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/Z2.MCL726s4.fishplot")


length(intersect(PDX3$V2 ,MCL726s1$V2))
length(intersect(MCL726s1$V2 ,PDX56$V2))

length(intersect(PDX3$V2 ,MCL726s2$V2))
length(intersect(MCL726s2$V2 ,PDX56$V2))

length(intersect(PDX3$V2 ,MCL726s3$V2))
length(intersect(MCL726s3$V2 ,PDX56$V2))

length(intersect(PDX3$V2 ,MCL726s4$V2))
length(intersect(MCL726s4$V2 ,PDX56$V2))

a = merge(MCL726s1[,1:4], MCL726s2[,1:4], by = c("V1", "V2"), all =T)
a = merge(a, MCL726s3[,1:4], by = c("V1", "V2"), all =T)
a = merge(a, MCL726s4[,1:4], by = c("V1", "V2"), all =T)


a[is.na(a)] <- 0
a$V3 <- a[,3]+ a[,5]+ a[,7]+ a[,9]
a$V4 <- a[,4]+ a[,6]+ a[,8]+ a[,10]
a <- a[,c("V1", "V2", "V3", "V4")]
a$V5 =a$V4 / (a$V3+a$V4)  *100

length(intersect(PDX3$V2 ,a$V2))
length(intersect(a$V2 ,PDX56$V2))

a <- a[a$V2 %in% c(PDX3$V3, PDX56$V2),]

MCL726s <- a

PDX3_cnv<- read.table("Y:/dna-seq-gatk-SNV-71-paired-tacc/call_CNV/sandbox/PDX3-tumor.cr.igv.seg/PDX3-tumor.cr.igv.seg", header = T)
PDX3_cnv$Segment_Mean <- round(2^(PDX3_cnv$Segment_Mean*2),2)*2

PDX56_cnv<- read.table("Y:/dna-seq-gatk-SNV-71-paired-tacc/call_CNV/sandbox/PDX56-tumor.cr.igv.seg/PDX56-tumor.cr.igv.seg", header = T)
PDX56_cnv$Segment_Mean <- round(2^(PDX56_cnv$Segment_Mean*2),2)*2

MCL726s_cnv<- read.table("Y:/dna-seq-gatk-SNV-71-paired-tacc/call_CNV/sandbox/MCL726s-tumor.cr.igv.seg/MCL726s-tumor.cr.igv.seg", header = T)
MCL726s_cnv$Segment_Mean <- round(2^(MCL726s_cnv$Segment_Mean*2),2)*2
```

```{r}
v1 = PDX3
v2 = PDX56
v3 = MCL726s

v1 <- v1[substring(v2$V1,1,3) %in% 'chr',]
v2 <- v2[substring(v2$V1,1,3) %in% 'chr',]
v3 <- v3[substring(v3$V1,1,3) %in% 'chr',]

cn1 = PDX3_cnv[,c(2:4,6)]
cn2 = PDX56_cnv[,c(2:4,6)]
cn3 = MCL726s_cnv[,c(2:4,6)]

samples = c("PDX3","PDX56","MCL726s")
```

```{r}
#cluster

sc_b = sciClone(vafs=list(v1,v2,v3), sampleNames=samples, copyNumberCalls=list(cn1,cn2,cn3), doClusteringAlongMargins=FALSE, maximumClusters=10, clusterMethod= 'binomial.bmm' )
#saveRDS(sc_b, file="rdata.0204PDX3_PDX56_MCL726_binomial.bmm.RDS")
sc_b  <- readRDS(file="rdata.0204PDX3_PDX56_MCL726_binomial.bmm.RDS")
sc <- sc_b
sc@clust$cluster.means

writeClusterTable(sc,"cluster")
writeClusterSummaryTable(sc,"cluster.summary")

sc.plot2d(sc,"figure5-sc_b-PDX3-PDX56-MCL726sp.pdf", singlePage=TRUE, scale=1.8)
```


```{r}

vafs = data.frame(cluster=sc@vafs.merged$cluster,
                  PDX3=sc@vafs.merged$PDX3.vaf,
                  PDX56=sc@vafs.merged$PDX56.vaf,
                  MCL726s=sc@vafs.merged$MCL726s.vaf,
                  stringsAsFactors=F)
vafs = vafs[!is.na(vafs$cluster) & vafs$cluster > 0,]
vafs_bmm <- vafs

plot.variant.clusters(vafs,
                      cluster.col.name = 'cluster',
                      show.cluster.size = FALSE,
                      cluster.size.text.color = 'blue',
                      vaf.col.names = samples,
                      vaf.limits = 30,
                      sample.title.size = 20,
                      violin = FALSE,
                      box = FALSE,
                      jitter = TRUE,
                      jitter.shape = 1,
                      jitter.color = clone.colors,
                      jitter.size = 3,
                      jitter.alpha = 1,
                      jitter.center.method = 'median',
                      jitter.center.size = 1,
                      jitter.center.color = 'darkgray',
                      jitter.center.display.value = 'none',
                      highlight = 'is.driver',
                      highlight.shape = 21,
                      highlight.color = 'blue',
                      highlight.fill.color = 'green',
                      highlight.note.col.name = 'gene',
                      highlight.note.size = 2,
                      order.by.total.vaf = FALSE)
# plot clusters pairwise-ly
plot.pairwise(vafs, col.names =samples,
              out.prefix = 'variants.pairwise.plot',
              colors = clone.colors)

# plot mean/median of clusters across samples (cluster flow)
plot.cluster.flow(vafs, vaf.col.names =samples,
                  sample.names = samples,
                  colors = clone.colors)
```

```{r}
res = infer.clonal.models(variants = vafs,
                          cluster.col.name = 'cluster',
                          vaf.col.names = samples,
                          subclonal.test = 'bootstrap',
                          subclonal.test.model = 'non-parametric',
                          num.boots = 1000,
                          founding.cluster = 1,
                          cluster.center = 'mean',
                          ignore.clusters = NULL,
                          min.cluster.vaf = 0.01,
                          sum.p = 0.05,
                          alpha = 0.05)


# new clonevol
res = convert.consensus.tree.clone.to.branch(res, branch.scale='sqrt')
```

```{r}
plot.clonal.models(res,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   #fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = 'output-sc_b-PDX3-PDX56-MCL726sp',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 8,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,2))

```
```{r}
## create a list of fish objects - one for each model (in this case, there's only one)
f = generateFishplotInputs(results=res)
fishes = createFishPlotObjects(f)

## plot each of these with fishplot


for (i in 1:length(fishes)){
    fish = layoutClones(fishes[[i]])
    fish = setCol(fish,f$clonevol.clone.colors)
    fishPlot(fish,shape="spline", title.btm="PatientID", cex.title=0.7,
             vlines=seq(1, length(samples)), vlab=samples[1:3], pad.left=0.5)
}

save(PDX3, PDX56, MCL726s, sc_b, file="rdata.0205PDX3_PDX59_MCL726_bmm.Rdata")
```
```{r}

vafs.merged <- sc@vafs.merged

#vafs.merged<- vafs.merged[!is.na(vafs.merged$cluster),]
```


annotion 

```{bash}
#cat  final.frame1-24.txt | awk -F" " '{print $1 " " $2 " " $3 " " $4  " " $5 " " $7}' > final.frame1-24forfishplot.txt

cat 01.PDX3_*_9_somatic_oncefiltered.ann.vcf | grep "missense_variant" | awk -F" " '{print $1 " "  $2 " " $4 " " $5 " " $8}'  |   awk -F"ANN=" '{print $1 " " $2 " " }'  | awk -F" " '{print $1 " "  $2 " " $3 " " $4 " " $6}' |  awk -F"|"   '{print $1 " "  $4 " " $10 " " $11 " "}' > PDX3.missense.ann

cat 01.PDX56_*_9_somatic_oncefiltered.ann.vcf | grep "missense_variant" | awk -F" " '{print $1 " "  $2 " " $4 " " $5 " " $8}'  |   awk -F"ANN=" '{print $1 " " $2 " " }'  | awk -F" " '{print $1 " "  $2 " " $3 " " $4 " " $6}' |  awk -F"|"   '{print $1 " "  $4 " " $10 " " $11 " "}' > PDX56.missense.ann

cat 01.MCL726s1_*_9_somatic_oncefiltered.ann.vcf | grep "missense_variant" | awk -F" " '{print $1 " "  $2 " " $4 " " $5 " " $8}'  |   awk -F"ANN=" '{print $1 " " $2 " " }'  | awk -F" " '{print $1 " "  $2 " " $3 " " $4 " " $6}' |  awk -F"|"   '{print $1 " "  $4 " " $10 " " $11 " "}' > MCL726s1.missense.ann

cat 01.MCL726s2_*_9_somatic_oncefiltered.ann.vcf | grep "missense_variant" | awk -F" " '{print $1 " "  $2 " " $4 " " $5 " " $8}'  |   awk -F"ANN=" '{print $1 " " $2 " " }'  | awk -F" " '{print $1 " "  $2 " " $3 " " $4 " " $6}' |  awk -F"|"   '{print $1 " "  $4 " " $10 " " $11 " "}' > MCL726s2.missense.ann

cat 01.MCL726s3_*_9_somatic_oncefiltered.ann.vcf | grep "missense_variant" | awk -F" " '{print $1 " "  $2 " " $4 " " $5 " " $8}'  |   awk -F"ANN=" '{print $1 " " $2 " " }'  | awk -F" " '{print $1 " "  $2 " " $3 " " $4 " " $6}' |  awk -F"|"   '{print $1 " "  $4 " " $10 " " $11 " "}' > MCL726s3.missense.ann

cat 01.MCL726s4_*_9_somatic_oncefiltered.ann.vcf | grep "missense_variant" | awk -F" " '{print $1 " "  $2 " " $4 " " $5 " " $8}'  |   awk -F"ANN=" '{print $1 " " $2 " " }'  | awk -F" " '{print $1 " "  $2 " " $3 " " $4 " " $6}' |  awk -F"|"   '{print $1 " "  $4 " " $10 " " $11 " "}' > MCL726s4.missense.ann
```


```{r}

MCL726s1.ann <- read.table("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/MCL726s1.missense.ann", fill = T)
MCL726s2.ann <- read.table("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/MCL726s2.missense.ann", fill = T)
MCL726s3.ann <- read.table("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/MCL726s3.missense.ann", fill = T)
MCL726s4.ann <- read.table("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/MCL726s4.missense.ann", fill = T)

MCL726s.ann <- rbind(MCL726s1.ann ,MCL726s2.ann )
MCL726s.ann <- rbind(MCL726s.ann ,MCL726s3.ann )
MCL726s.ann <- rbind(MCL726s.ann ,MCL726s4.ann )

PDX3.ann <- read.table("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/PDX3.missense.ann", fill = T)

PDX56.ann <- read.table("C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/PDX56.missense.ann", fill = T)

PDX.ann <- rbind(PDX3.ann ,PDX56.ann )
colnames(PDX.ann)[1:2] <- c("chr", "st")
PDX.ann  <- PDX.ann[!duplicated(PDX.ann$st), ]




merged<- merge(vafs.merged, PDX.ann, by = c("chr", "st"))
merged  <- merged[!duplicated(merged$st), ]
merged  <- merged[!is.na(merged$cluster), ]

merged[,c("chr","st", "PDX3.vaf", "PDX56.vaf","MCL726s.vaf","V6", "V8")]
b <- merged[!(merged$PDX3.vaf == 0 & merged$MCL726s.vaf ==0 & !merged$PDX56.vaf ==0),][,c("chr","st", "PDX3.vaf", "PDX56.vaf","MCL726s.vaf","V6", "V8", "cluster")]
b[order(b$cluster),]
```
```{r}
vafs.merged <- sc@vafs.merged
a = merge(MCL726s1[,1:4], MCL726s2[,1:4], by = c("V1", "V2"), all =T)
a = merge(a, MCL726s3[,1:4], by = c("V1", "V2"), all =T)
a = merge(a, MCL726s4[,1:4], by = c("V1", "V2"), all =T)


a[is.na(a)] <- 0
a$V3 <- a[,3]+ a[,5]+ a[,7]+ a[,9]
a$V4 <- a[,4]+ a[,6]+ a[,8]+ a[,10]
a <- a[,c("V1", "V2", "V3", "V4")]
a$V5 =a$V4 / (a$V3+a$V4)  *100

a <- a[!a$V2 %in% vafs.merged$st,]
#colnames(a)[1:2] <- c("chr", "st")
merged.a <- merge(a, MCL726s.ann, by = c("V1", "V2"))
merged.a$depth <- merged.a$V3.x + merged.a$V4.x
merged.a  <- merged.a[!duplicated(merged.a$V8), ]

merged.a[merged.a$V5.x > 10 & merged.a$V5.x < 50 & merged.a$depth > 100,]   # for MCL726s speically
```

vafs.merged <- sc@vafs.merged
vafs.merged <- vafs.merged[vafs.merged$st %in%  PDX3.ann$V2,]

vafs.merged <- vafs.merged[vafs.merged$PDX3.depth > 100 ,]
t1<- vafs.merged[is.na(vafs.merged$cluster),][,c("cluster","PDX3.vaf", "PDX56.vaf","MCL726s.vaf")]
colnames(t1) <- c("cluster","PDX3", "PDX56","MCL726s")
t1 <- t1[t1$PDX3 > 15 & t1$PDX56 ==0 & t1$MCL726s ==0,]
t1$cluster <- 8
vafs3 <- rbind(vafs2,t1)

```{r}
vafs = data.frame(cluster=sc@vafs.merged$cluster,
                  PDX3=sc@vafs.merged$PDX3.vaf,
                  PDX56=sc@vafs.merged$PDX56.vaf,
                  MCL726s=sc@vafs.merged$MCL726s.vaf,
                  stringsAsFactors=F)
vafs = vafs[!is.na(vafs$cluster) & vafs$cluster > 0,]
```


add ann
```{r}
vafs.merged <- sc@vafs.merged
vafs.merged$PDX3 <- vafs.merged$PDX3.vaf
vafs.merged$PDX56 <- vafs.merged$PDX56.vaf
vafs.merged$MCL726s <- vafs.merged$MCL726s.vaf

vafs  <- vafs.merged[,c("chr", "st", "cluster","PDX3.vaf", "PDX56.vaf","MCL726s.vaf")]
colnames(vafs) <- c("chr", "st","cluster","PDX3", "PDX56","MCL726s")
vafs = vafs[!is.na(vafs$cluster) & vafs$cluster > 0,]
vafs
vafs <- merge(vafs,PDX.ann[,c("chr","st", "V6")], by = c("chr","st") ,all.x=T)
colnames(vafs) <- c("chr", "st","cluster","PDX3", "PDX56", "MCL726s", "gene")
vafs



#vafs <- vafs.merged[,c("cluster","PDX3.vaf", "PDX56.vaf","MCL726s.vaf")]
#colnames(vafs) <- c("cluster","PDX3", "PDX56","MCL726s")

add<- merged.a[merged.a$V5.x > 10 & merged.a$V5.x < 50 & merged.a$depth > 100,]
add$cluster <- 7
add$PDX3 <- 0
add$PDX56 <- 0
add <- add[,c("V1", "V2", "cluster","PDX3", "PDX56", "V5.x", "V6")]
colnames(add) <- c("chr", "st","cluster","PDX3", "PDX56", "MCL726s", "gene")

vafs2 <- rbind(vafs,add)

vafs2  <- vafs2[order(vafs2$cluster),]
vafs2
vafs2$is.driver <- FALSE
vafs2[vafs2$cluster == 1,]$is.driver <- TRUE
vafs2
str(vafs)
str(vafs2)
```
```{r}
plot.variant.clusters(vafs2,
                      cluster.col.name = 'cluster',
                      show.cluster.size = FALSE,
                      cluster.size.text.color = 'blue',
                      vaf.col.names = samples,
                      vaf.limits = 30,
                      sample.title.size = 20,
                      violin = FALSE,
                      box = FALSE,
                      jitter = TRUE,
                      jitter.shape = 1,
                      jitter.color = clone.colors,
                      jitter.size = 3,
                      jitter.alpha = 1,
                      jitter.center.method = 'median',
                      jitter.center.size = 1,
                      jitter.center.color = 'darkgray',
                      jitter.center.display.value = 'none',
                      highlight = 'is.driver',
                      highlight.shape = 21,
                      highlight.color = 'blue',
                      highlight.fill.color = 'green',
                      highlight.note.col.name = 'gene',
                      highlight.note.size = 2,
                      order.by.total.vaf = FALSE)
# plot clusters pairwise-ly
plot.pairwise(vafs, col.names =samples,
              out.prefix = 'variants.pairwise.plot',
              colors = clone.colors)

# plot mean/median of clusters across samples (cluster flow)
plot.cluster.flow(vafs2, vaf.col.names =samples,
                  sample.names = samples,
                  colors = clone.colors)
```

```{r}
res = infer.clonal.models(variants = vafs2,
                          cluster.col.name = 'cluster',
                          vaf.col.names = samples,
                          subclonal.test = 'bootstrap',
                          subclonal.test.model = 'non-parametric',
                          num.boots = 1000,
                          founding.cluster = 1,
                          cluster.center = 'mean',
                          ignore.clusters = NULL,
                          min.cluster.vaf = 0.01,
                          sum.p = 0.05,
                          alpha = 0.05)


# new clonevol
res = convert.consensus.tree.clone.to.branch(res, branch.scale='sqrt')
```

```{r}

plot.clonal.models(res,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   
                   
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.2,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   
                   
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 0.5,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 0.5,
                   mtcab.node.text.size = 0.5,
                   
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   
                   
                   # output figure parameters
                   out.dir="C:/Users/qcai1/Downloads/sciclone/sciclone-meta-master/sciclone-meta-master/manuscript/data/figure5/output2",
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 6,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,3))

```


```{r}
plot.clonal.models(res,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,  #1
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = 'output-sc_b-PDX3-PDX56-MCL726sp',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 6,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,2))

```
```{r}
## create a list of fish objects - one for each model (in this case, there's only one)
f = generateFishplotInputs(results=res)
fishes = createFishPlotObjects(f)

## plot each of these with fishplot


for (i in 1:length(fishes)){
    fish = layoutClones(fishes[[i]])
    fish = setCol(fish,f$clonevol.clone.colors)
    fishPlot(fish,shape="spline", title.btm="PatientID", cex.title=0.7,
             vlines=seq(1, length(samples)), vlab=samples[1:3], pad.left=0.5)
}

save(PDX3, PDX56, MCL726s, sc_b,add, vafs2, file="rdata.0205PDX3_PDX59_MCL726_bmm.Rdata")
```


```{r}
vafs2
```

