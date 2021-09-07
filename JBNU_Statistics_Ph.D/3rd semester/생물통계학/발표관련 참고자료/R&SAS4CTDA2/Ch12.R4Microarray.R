
###################################################
### chunk number 2:  eval=FALSE
###################################################
## #line 103 "R4Microarray.rnw"
## source("http://bioconductor.org/biocLite.R")
## biocLite()


###################################################
### chunk number 3: dat
###################################################
#line 185 "R4Microarray.rnw"
dat = read.table("MA.txt", header=T)
dat[,1:2]


###################################################
### chunk number 4: 
###################################################
#line 214 "R4Microarray.rnw"
library(affy)


###################################################
### chunk number 5: read.dat
###################################################
#line 220 "R4Microarray.rnw"
Dat = ReadAffy(filenames=paste(dat$Samples,".cel",sep=""))


###################################################
### chunk number 6: dat.print
###################################################
#line 226 "R4Microarray.rnw"
Dat


###################################################
### chunk number 7: 
###################################################
#line 233 "R4Microarray.rnw"
annotation(Dat)


###################################################
### chunk number 8: 
###################################################
#line 238 "R4Microarray.rnw"
phenoData(Dat)


###################################################
### chunk number 9: 
###################################################
#line 243 "R4Microarray.rnw"
# the original sample name
sampleNames(Dat)
# change it to tumor type 
sample.names = dat$Tumor
colnames(exprs(Dat)) = sample.names


###################################################
### chunk number 10: 
###################################################
#line 254 "R4Microarray.rnw"
e = exprs(Dat)
dim(e)
nrow(Dat)*ncol(Dat)


###################################################
### chunk number 11: 
###################################################
#line 261 "R4Microarray.rnw"
# gene numbers
gnames = geneNames(Dat)
# the total number of genes
length(gnames)
# print the first 20 genes
gnames[1:20]


###################################################
### chunk number 12: MA.fig4box1
###################################################
#line 282 "R4Microarray.rnw"
boxplot(Dat,col=c(rep("Green",6),rep("Blue",16),
			 rep("red",27)))


###################################################
### chunk number 13: 
###################################################
#line 288 "R4Microarray.rnw"
#line 282 "R4Microarray.rnw#from line#288#"
boxplot(Dat,col=c(rep("Green",6),rep("Blue",16),
			 rep("red",27)))
#line 289 "R4Microarray.rnw"


###################################################
### chunk number 14: MA.fig4MAplot1
###################################################
#line 310 "R4Microarray.rnw"
MAplot(Dat[,c(1,8)], pair=T) 


###################################################
### chunk number 15: 
###################################################
#line 316 "R4Microarray.rnw"
#line 310 "R4Microarray.rnw#from line#316#"
MAplot(Dat[,c(1,8)], pair=T) 
#line 317 "R4Microarray.rnw"


###################################################
### chunk number 16: rma
###################################################
#line 338 "R4Microarray.rnw"
# call RMA function
eset = rma(Dat)


###################################################
### chunk number 17: 
###################################################
#line 386 "R4Microarray.rnw"
# dimension 
dim(exprs(eset))
# print the first 20 genes for three patients
exprs(eset)[1:20,1:3]


###################################################
### chunk number 18: MA.fig4box2
###################################################
#line 397 "R4Microarray.rnw"
boxplot(data.frame(exprs(eset)),
col=c(rep("Green",6),rep("Blue",16), rep("red",27)))


###################################################
### chunk number 19: 
###################################################
#line 404 "R4Microarray.rnw"
#line 397 "R4Microarray.rnw#from line#404#"
boxplot(data.frame(exprs(eset)),
col=c(rep("Green",6),rep("Blue",16), rep("red",27)))
#line 405 "R4Microarray.rnw"


###################################################
### chunk number 20: 
###################################################
#line 478 "R4Microarray.rnw"
library(limma)


###################################################
### chunk number 21: 
###################################################
#line 487 "R4Microarray.rnw"
# make the design matrix
design <- model.matrix(~ 0+factor(c(rep("apocrine",6), 
	rep("basal",16), rep("luminal",27))))
# label the design
colnames(design) <- c("apocrine", "basal", "luminal")
# print the design
design


###################################################
### chunk number 22: 
###################################################
#line 498 "R4Microarray.rnw"
# fit the linear model for all genes
fit = lmFit(eset, design)
# make a contrast for pair-wise comparisons
cont.matrix = makeContrasts(Comp2to1=basal-apocrine,
Comp3to1=luminal-apocrine,Comp3to2=luminal-basal,levels=design)
# then contrast fit
fit2 = contrasts.fit(fit, cont.matrix)
# then call ebayes
fit2 = eBayes(fit2)


###################################################
### chunk number 23: 
###################################################
#line 512 "R4Microarray.rnw"
options(digits=3)
topTable(fit2, coef=1, adjust="BH")


###################################################
### chunk number 24: MA.fig4volcano
###################################################
#line 543 "R4Microarray.rnw"
volcanoplot(fit2,coef=2, highlight=10)
abline(v=c(-1,1),col="red")
ilogit = function(p) exp(p)/(1+exp(p))
abline(h=ilogit(0.05), col="blue")


###################################################
### chunk number 25: 
###################################################
#line 552 "R4Microarray.rnw"
#line 543 "R4Microarray.rnw#from line#552#"
volcanoplot(fit2,coef=2, highlight=10)
abline(v=c(-1,1),col="red")
ilogit = function(p) exp(p)/(1+exp(p))
abline(h=ilogit(0.05), col="blue")
#line 553 "R4Microarray.rnw"


###################################################
### chunk number 26: 
###################################################
#line 566 "R4Microarray.rnw"
# call decideTests function
results = decideTests(fit2)
# classification counts in a Venn diagram
venn    = vennCounts(results)
# print the venn table
print(venn)


###################################################
### chunk number 27: MA.fig4venn
###################################################
#line 577 "R4Microarray.rnw"
vennDiagram(results,include=c("up","down"),
		counts.col=c("red","green"))


###################################################
### chunk number 28: 
###################################################
#line 584 "R4Microarray.rnw"
#line 577 "R4Microarray.rnw#from line#584#"
vennDiagram(results,include=c("up","down"),
		counts.col=c("red","green"))
#line 585 "R4Microarray.rnw"


