############################################################################################################
####################  "Gene Affinity Network Fusion Model"  ###############################################

########  Step1:  data loading -> data collapse -> data filtering                                       #
######  Step2:  gene fusing ->  gene clustering -> eigen-gene extraction -> integrative correlation   ##
####  Step3:  module robustness analysis  -> module wellness validation (silhouvette measurement)   ###
##  Step4:  Visulization of modules correlation plot, cytoscape output                            ####
  
####################################################################################################



# Step1: data loading
#-----------------------------------------------------
#Importing expr/meth/pheno data for data integration trial
library(GEOquery)
library(oligo)
library(SNFtool)
library(lumi)
library(ggplot2)
library(data.table)
library(magrittr)
library(reshape2)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(caret)
library(Rtsne)
library(WGCNA)
library(genefilter)
library(dplyr)
library(splitstackshape)




getwd()
setwd("/global/home/hpc4009/MESAdata/dataIntProject/")


folderPath = "~/MESAdata/dataIntProject"


#Load or retrieve downloaded GEO meth data
loadMethData <- function() {
  setwd(folderPath)
  if (file.exists("methEset.txt")) {
    methEset <- as.data.frame(read.table("methEset.txt"))
  } else {
    #Must have downloaded GSE40576_series_matrix.txt.gz from Series Matrix File(s) link on https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40576
    x <- getGEO(filename = "GSE40576_series_matrix.txt.gz", getGPL = FALSE)
    write.exprs(x, file = "methEset.txt")
    methEset <- as.data.frame(exprs(x))
  }
  
  #Convert methylation IDs to match up with expression IDs
  conversion <- read.table("idConversionIndex.txt")
  colnames(methEset) <- conversion$V2
  methEset = methEset[,order(colnames(methEset))]
  methEset = cbind(methEset[,62:194], methEset[1:61])
  
  return(methEset)
}

#Load or retrieve downloaded GEO expr data
loadExprData <- function() {
  setwd(folderPath)
  if (file.exists("exprEset.txt")) {
    exprEset <- as.data.frame(read.table("exprEset.txt"))
  } else {
    #Must have downloaded GSE40732_series_matrix.txt.gz from Series Matrix File(s) link on https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40732
    x <- getGEO(filename = "GSE40732_series_matrix.txt.gz", getGPL = FALSE)
    write.exprs(x, file = "exprEset.txt")
    exprEset <- as.data.frame(exprs(x))
  }
  return(exprEset)
}

#Setting up phenotype data
phenotypeData <- function() {
  setwd(folderPath)
  pheno <- read.table("phenotypeSUM.txt", header = TRUE)
  pheno <- pheno[,2:9]
  asthma <- factor(pheno$asthma)
  subtype <- factor(pheno$subtype)
  sex <- factor(pheno$sex)
  race_aa <- factor(pheno$race_aa)
  race_his <- factor(pheno$race_hispanic)
  race_dom_or_hait <- factor(pheno$dominican_or_haitian)
  race_other <- factor(pheno$race_other)
  design <- model.matrix(~asthma+subtype+sex+race_aa+race_his+race_dom_or_hait+race_other)
  phenotype <- design[,2:10]
  phenotype <- cbind(phenotype, pheno$age)
  colnames(phenotype) <- c("asthma","subtype_lungfunction","subtype_PC20","subtype_reversible","sex_male","race_aa","race_hisp","race_dom/hait", "race_other","age")
  end <- phenotype
  return(end)
}


#Remove sex chromosome probes for Methylation
makeAutosomalMeth <- function(matrix){
  require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  good.probes <- rownames(ann450k)[-which(ann450k$chr == "chrX" | ann450k$chr == "chrY")]
  matrix[rownames(matrix) %in% good.probes,]
}

#Remove sex chromosome probes for Expression
makeAutosomalExpr <- function(matrix){
  setwd(folderPath)
  annExpr <- read.table("expr_annotation.txt", header = TRUE, sep = "|", row.names = 1, fill = TRUE)
  good.probes <- rownames(annExpr)[-which(annExpr$CHROMOSOME == "chrX" | annExpr$CHROMOSOME == "chrY")]
  matrix[rownames(matrix) %in% good.probes,]
}



############################################################## All the data loadings
phenotype <- phenotypeData()
methTable = loadMethData() #size: 485460    194
exprTable = loadExprData() #size:  45033    194
methTableAutosomal <- makeAutosomalMeth(methTable) #size: 473814    194
exprTableAutosomal <- makeAutosomalExpr(exprTable) #size:  43110    194
exprTable <- as.data.frame(exprTable)
methTable <- as.data.frame(methTable)
exprTableAutosomal <- as.data.frame(exprTableAutosomal)
methTableAutosomal <- as.data.frame(methTableAutosomal)
phenotype <- as.data.frame(phenotype)






# Step1: data collapse
#Methylation Gene Level collapsing
#-----------------------------------------------------
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(methTableAutosomal),ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
methyRowName <- rownames(methTableAutosomal)
annotfile <- as.data.frame(cbind(methyRowName,ann450kSub$Relation_to_Island,ann450kSub$UCSC_RefGene_Name))
varianceMeth <- apply(methTableAutosomal, 1, var)
varianceMeth <- as.data.frame(varianceMeth)
annotfile <- cbind(annotfile,varianceMeth)
annotfile[annotfile==""] <- NA
annotfile = na.omit(annotfile)
dim(annotfile)
colnames(annotfile) = c("ProbeID","RelationToIsland","GeneName","ProbeVariance")
annotfile1<-(cSplit(annotfile, "GeneName", ";", "long", makeEqual = FALSE))
probeToGeneMethy <- annotfile[!duplicated(annotfile$ProbeID),]
dim(probeToGeneMethy) # original methylation size: 473,814; after removed the duplicated probes: 356,080


## Collapse probes to gene level based on Island -> Shore -> and then -> highest variance (ordered with rank)
testData <- cbind(annotfile1,matrix(1,nrow(annotfile1),1))
testData[which(annotfile1$RelationToIsland=="Island"),5] <- 4
testData[which(annotfile1$RelationToIsland=="S_Shore" | annotfile1$RelationToIsland=="N_Shore"),5] <- 3
testData[which(annotfile1$RelationToIsland=="S_Shelf" | annotfile1$RelationToIsland=="N_Shelf" | annotfile1$RelationToIsland=="OpenSea"),5] <- 1
colnames(testData)[5] <- "Fun"
df_arranged <- arrange(testData,by=GeneName,desc(Fun),desc(ProbeVariance))
probeTogeneM <- distinct(df_arranged, by=GeneName,.keep_all = T)

indexProbeToGene <- match(probeTogeneM$ProbeID,rownames(methTableAutosomal))
Final_autosomalMeth <- methTableAutosomal[indexProbeToGene,]
rownames(Final_autosomalMeth) <- as.matrix(probeTogeneM[,3])
dim(Final_autosomalMeth) #size: 20,271   194



# Step1: data collapse
#Expression Gene Level collapsing
#-----------------------------------------------------
annot <- read.delim("GPL16025_100718_HG18_opt_expr.ngd", header = FALSE, sep = "", colClasses = c("character", "NULL", "NULL","NULL","NULL","NULL","NULL","NULL","character","NULL","NULL"))
annot <- annot[annot[,3] != "",] #removes blank rows that were created
annot <- annot[,-3]
colnames(annot) <- c("Probe", "Gene")

# get all probes and associating gene names
for (i in 1:nrow(annot)) {
  annot[i,2] <- substr(annot[i,2],2,nchar(annot[i,2])-2)
  if (i %% 1000 ==0) {
    print(i)
  }
}

IndExpr <- match(rownames(exprTableAutosomal),annot[,1])
FinAnnotExpr <- annot[IndExpr,]
#Median absolute deviation: 
#mad_among_Expr_id <- apply(exprTableAutosomal, 1, mad)
var_among_Expr_id <- apply(exprTableAutosomal, 1, var)
Expr_var_id <- data.frame(Probe=rownames(exprTableAutosomal),probeValue=var_among_Expr_id)
Merged_Expr_id_geneVAR <- merge(FinAnnotExpr,Expr_var_id,by = "Probe")

indeNAExp <- which(Merged_Expr_id_geneVAR$Gene=="N/A",arr.ind = T)
Unique_Expr_genes <- Merged_Expr_id_geneVAR[-indeNAExp,]
Unique_Expr_genes<-na.omit(Unique_Expr_genes)

df_arranged1 <- arrange(Unique_Expr_genes,by=Gene,desc(probeValue))
probeTogeneE <- distinct(df_arranged1, by=Gene,.keep_all = T)
dim(probeTogeneE)
indexProbeToGene <- match(probeTogeneE$Probe,rownames(exprTableAutosomal))
Final_autosomalExpr <- exprTableAutosomal[indexProbeToGene,]
rownames(Final_autosomalExpr) <- as.matrix(probeTogeneE[,2])
dim(Final_autosomalExpr) #current: 22,653   194;



# Step1: data filtering
# Pick up the overlapping genes
#-----------------------------------------------------
ComGeneInd <- match(rownames(Final_autosomalMeth),rownames(Final_autosomalExpr))
ComGeneInd <- na.omit(ComGeneInd)
length(ComGeneInd)
Final_Expr_AutosomalData <- Final_autosomalExpr[ComGeneInd,]
ComGeneIndM <- match(rownames(Final_Expr_AutosomalData),rownames(Final_autosomalMeth))
length(ComGeneIndM)  ## 15338 genes in common
Final_Meth_AutosomalData <- Final_autosomalMeth[ComGeneIndM,]
dim(Final_Expr_AutosomalData)
dim(Final_Meth_AutosomalData)


# Step1: data filtering !!!!!! (overlaping + filtering and then overlaping)
# Quantile based var-filtering
#-----------------------------------------------------
objectExpr <- new("ExpressionSet", exprs=as.matrix(Final_Expr_AutosomalData))
objectMeth <- new("ExpressionSet", exprs=as.matrix(Final_Meth_AutosomalData))
dim(objectMeth)
filtered_expr<-varFilter(eset = objectExpr,var.func = var,filterByQuantile = TRUE,var.cutoff = 0.5) # taking top 50% variance of expr
filtered_meth<-varFilter(eset = objectMeth,var.func = var,filterByQuantile = TRUE,var.cutoff = 0.5) # taking top 50% variance of methy
Expr_quantile_filtered <- filtered_expr@assayData$exprs
Meth_quantile_filtered <- filtered_meth@assayData$exprs
dim(Expr_quantile_filtered)
dim(Meth_quantile_filtered)

########## Using filtered Data to build target data_list
EqfInd <- rownames(Expr_quantile_filtered)
MqfInd <- rownames(Meth_quantile_filtered)

FinOverIndE <- match(MqfInd,EqfInd)
FinOverIndE <- na.omit(FinOverIndE)
length(FinOverIndE)
ExprData_filtered <- Expr_quantile_filtered[FinOverIndE,]

FinOverIndM <- match(rownames(ExprData_filtered),rownames(Meth_quantile_filtered))
FinOverIndM <- na.omit(FinOverIndM)
length(FinOverIndM) ## proportion: 3802/15338=0.2478811;

#Final_Expr_AutosomalData <- as.data.frame(Final_Expr_AutosomalData)
MethData_filtered <- Meth_quantile_filtered[FinOverIndM,]
dim(ExprData_filtered)
dim(MethData_filtered)

CellType1Seq <- seq(1, 72, by=3)
CellType2Seq <- seq(2, 72, by=3)
CellType3Seq <- seq(3, 72, by=3)

ExprData_filtered1 <- ExprData_filtered[,CellType1Seq]
ExprData_filtered2 <- ExprData_filtered[,CellType2Seq]
ExprData_filtered3 <- ExprData_filtered[,CellType3Seq]
MethData_filtered1 <- MethData_filtered[,CellType1Seq]
MethData_filtered2 <- MethData_filtered[,CellType2Seq]
MethData_filtered3 <- MethData_filtered[,CellType3Seq]
dim(ExprData_filtered1)
dim(MethData_filtered1)
phenotype1 <- phenotype[CellType1Seq,]
phenotype2 <- phenotype[CellType2Seq,]
phenotype3 <- phenotype[CellType3Seq,]


###------------------------------------------------------------------------------------- Step 1 finished 


# Step2.1: gene fusing ----------------------------------------------------
# SNF-----------------------------------------------------
library(SNFtool)

K <- 35 ##number of neighbors, must be greater than 1. usually (20~50)
alpha <- 0.5 ##hyperparameter, usually (0.3~0.8)
T <-  20 ###Number of Iterations, usually (10~50)

data_filtered_list <- list(ExprData_filtered1, MethData_filtered1)
names(data_filtered_list) <- c("Expression", "Methylation")
dist2 <- SNFtool::dist2
data_dist_list <- lapply(X = data_filtered_list,
                         function(x){(dist2(as.matrix(x),as.matrix(x)))^(1/2)})
diag(data_dist_list$Expression) <- 0
diag(data_dist_list$Methylation) <- 0

data_aff_list <-  lapply(X = data_dist_list,
                         function(x){affinityMatrix(x,K,alpha)})
W1_1 <-  SNF(data_aff_list, K, T)
rm(data_dist_list, data_aff_list)
dim(W1_1)



# Step2.2: gene clustering ------------------------------------------------

library(WGCNA)
geneTreec <- hclust(as.dist(1-W1_1), method = "average")              ## WGCNA
dynamicModsHybrid = cutreeDynamic(dendro = geneTreec, distM = 1-W1_1,method="hybrid",deepSplit = 4, pamRespectsDendro = FALSE,minClusterSize = 30)
dynamicColorsHybrid = labels2colors(dynamicModsHybrid)
plotDendroAndColors(geneTreec, dynamicColorsHybrid, "Hybrid Tree Cut",dendroLabels = FALSE, hang = 0.05,addGuide = TRUE,
                    guideHang = 0.05, main = "Gene modules dendrogram of fused affinity matrix by GANF")
length(table(dynamicModsHybrid))



# Step2.3: eigen-gene extraction ------------------------------------------

# write the Rtsne loop operation into a function
# x represents genes' colors (e.g. dynamicColorsHybrid for W, W is affinity matrix)
# y represents unique(colors) (unique colors' names): ModuleList for W1 
# z represents original data (for example: expression data or methylation data)
DimensionRRtsne <- function(x, y, z){
  GeneData <- matrix(0,nrow = ,ncol = ncol(z))
  RTsneResult <- matrix(0,nrow = ncol(z),ncol = length(y))
  GeneNamesList <- rownames(z)
  print("Extracting eigen-gene from module:")
  for (i in 1:length(y)) {
    GeneinModule <-GeneNamesList[which(x == y[i])] # tells which genes belong to the specific color
    GeneData <- z[match(GeneinModule,rownames(z)),] # extracting those genes' data with specific color
    set.seed(42)
    pplex <- floor((nrow(t(GeneData))-1)/3)
    Result <- Rtsne(t(GeneData), pca = F, perplexity = pplex, theta = 0.0, dims = 1)
    RTsneResult[,i] <- Result$Y
    GeneinModule <- character()
    GeneData <- matrix(0,nrow = ,ncol = ncol(z))
    cat(i,"-")
  }
  RTsneResult <- as.data.frame(RTsneResult)
  colnames(RTsneResult) <- y
  return(RTsneResult)
}

ModuleList <- unique(dynamicColorsHybrid)
length(ModuleList)  

RTsneDataExp <- DimensionRRtsne(dynamicColorsHybrid, ModuleList, data_filtered_list$Expression)
RTsneDataMet <- DimensionRRtsne(dynamicColorsHybrid, ModuleList, data_filtered_list$Methylation)


# Step2.4: integrative correlation ------------------------------------------

# use eigenvectors to calculate the multiple correlation coeffient (correlation between expr, methy and pheno, this can be expand to more than 3 inputs) 
# http://www.real-statistics.com/correlation/multiple-correlation/
## Write into a function that doing multiple correlation coeffient calculation
# x is eigenvectors from expression (RTsneDataExp_wgcnaT1); y is eigenvector from methylation (RTsneDataMet_wgcnaT1); z is target phenotypes (disease on first column)
MultiCorCal <- function(x,y,z){
  r1 <- vector() # cor of xz
  r2 <- vector() # cor of yz
  r3 <- vector() # cor of xy
  r <- vector()  # target cor
  for (i in 1:ncol(x)) {
    r1 <- cor(x[,i],z[,1])
    r2 <- cor(y[,i],z[,1])
    r3 <- cor(x[,i],y[,i])
    r[i] <- sqrt(((r1)^2+(r2)^2-2*r1*r2*r3)/(1-(r3)^2))
  }
  return(r)
}

Rvalue <- MultiCorCal(RTsneDataExp,RTsneDataMet,phenotype1)
names(Rvalue) <- unique(dynamicColorsHybrid)
Pvalue <- corPvalueStudent(Rvalue,ncol(data_filtered_list$Expression))
names(Pvalue) <- unique(dynamicColorsHybrid)
GeneNamesList1 <- colnames(W1_1)
BestAsthmaModule <- GeneNamesList1[which(dynamicColorsHybrid=="lightcyan1")]
length(BestAsthmaModule)
View(BestAsthmaModule)
length(which(Pvalue<=0.05))

### Showing the best module of gene, which has the strongest correlation with asthma
# x is module correlation values from previous step, 
# y is original data with rownames of genes,   z is module colors
# returning best module of gene names
TopSignificantModule <- function(x,y,z){
  BestModule <- list()
  names(x) <- unique(z)
  indexTargetM <- names(x[which(x==max(x),arr.ind = T)])
  GeneNames <- rownames(y)
  BestModule$genes <- GeneNames[which(z==indexTargetM)]
  BestModule$Correlation <- max(x)
  cat("Gene numbers: ",length(BestModule$genes))
  cat("\n")
  print(BestModule)
  return(BestModule)
}

BestModule <- TopSignificantModule(Rvalue,data_filtered_list$Expression,dynamicColorsHybrid)





# Step3.1:  module robustness analysis ------------------------------------
# Step3.1:  module robustness analysis
# x is the best module you've found in prevous step;
# y is the testing arounds;
# z is the similarity matrix (W here); 
# m is the testing range; n is the datasets (as a list which contain expression and methylation data in this case)
# o is phenotype
RunningT <- 10000
testRange <- seq(from=20,to=50)
CheckRobustness <- function(x,y,z,m,n,o){
  set.seed(NULL)
  cat("Input running time:",y, " .\n Input testing parameter range:",m[1], "to", m[length(m)],"\n")
  targetList <- list()
  testList <- matrix(data = NA, nrow = y, ncol = 3)
  parameterList <- matrix(data = NA,nrow = length(m), ncol = 2)
  geneTreec <- hclust(as.dist(1-z), method = "average")
  randomPickup <- sample(m,y,replace = TRUE)
  k=1
  for (i in m[1]:m[length(m)]) {
    dynamicModsHybrid = cutreeDynamic(dendro = geneTreec, distM = 1-z,method="hybrid",deepSplit = 4, pamRespectsDendro = FALSE,minClusterSize = i)
    dynamicColorsHybrid = labels2colors(dynamicModsHybrid)
    cat("Building up all tunning table, parameter choosing from:",m[1], " to", m[length(m)],".\n Test for",i,"\n")
    parameterList[k,1] <- i
    ModuleList <- unique(dynamicColorsHybrid)
    cat("\n Starting calculating Eigen-genes from expression:\n")
    RTsneDataExp <- DimensionRRtsne(dynamicColorsHybrid, ModuleList, n$Expression)
    cat("\n Starting calculating Eigen-genes from methylation:\n")
    RTsneDataMet <- DimensionRRtsne(dynamicColorsHybrid, ModuleList, n$Methylation)
    Rvalue <- MultiCorCal(RTsneDataExp,RTsneDataMet,o)
    names(Rvalue) <- unique(dynamicColorsHybrid)
    parameterList[k,2] <- max(Rvalue)
    colnames(testList) <- c("RunningTime", "minClusterSize", "MaxCorrelation")
    cat("\n --------------------------------------------------Done for testing parameter:", i, ".\n")
    k <- k+1
    set.seed(NULL)
  }
  for (j in 1:y) {
    testList[j,1] <- j
    testList[j,2] <- randomPickup[j]
    testList[j,3] <- parameterList[which(parameterList[,1]==randomPickup[j]),2]
  }
  cat("\n --------------------------------------------------Done for running times:", j, ".\n")
  RobustModuleList <- as.data.frame(table(testList[,3]))
  RobustModuleList <- RobustModuleList[order(RobustModuleList$Freq, decreasing = TRUE),]
  if (x$Correlation==RobustModuleList[1,1]) {
    ratio <- length(which(testList[,3]==x$Correlation,arr.ind = T))/y
    cat("\n Yes, your module is robust and it appears ", ratio*100,"% times. \n")
  } else {
    cat("\n No, your module is not the best. \n The best module is: \n")
    print(testList[which(testList[,3]==RobustModuleList[1,1],arr.ind = T)[1],])
  }
  targetList$ParameterList <- parameterList
  targetList$BestModule <- RobustModuleList
  targetList$TunningProcess <- testList
  return(targetList)
}

# test our best module found in the previous step and setting up 
RobustnessList <- CheckRobustness(BestModule,RunningT,W,testRange,data_filtered_list,phenotype)

RobustTableGanf <- as.data.frame(RobustnessList$TunningProcess)
RobustTableGanf$Output <- rep(1,times=nrow(RobustTableGanf))
notRobIDs <- which(RobustTableGanf$MaxCorrelation!=RobustTableGanf$MaxCorrelation[1],arr.ind = T)
RobustTableGanf$Output[notRobIDs] <- rep(0,times=length(notRobIDs))
View(RobustTableGanf)

# visualization of robustness test
RobustTableGanf %>% 
  mutate(Output1 = case_when(Output==1 ~ "Yes", Output==0 ~ "No")) %>% 
  ggplot(aes(x=RunningTime, y=Output1))+
  geom_jitter(height=.1, size=.5)+theme_bw()+
  theme(
    panel.border = element_blank(),axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size=16),
    plot.subtitle = element_text(hjust = 0.5,size=12)
  )+
  labs(title = "Robustness test on the most relevant module"
       , subtitle = "the most relevant module stays stable with 76.72% chance (Yes), and 23.28% break apart (No)",
       y="Robustness test outcome"
       )
  +
  scale_y_continuous(breaks = seq(0, 1, .2))

# Step3.2: module wellness validation (silhouvette measurement) -------------

# calculate the silhouvette values for each node
# input1: total module colours; input2: similarity matrix
SilhovtForModule <- function(ClusterColors, W1){
  GeneNames <- rownames(W1)
  SumarySilhoutM <- vector()
  # SumarySilhoutM <- matrix(0,1000, length(unique(ClusterColors)))
  for (i in 1:length(unique(ClusterColors))) {
    moduleList <- unique(ClusterColors)
    Module_i <- GeneNames[which(ClusterColors==unique(ClusterColors)[i])]
    IndxModule_i <- match(Module_i,rownames(W1))
    SilhovtMatrix_i <- matrix(0,length(IndxModule_i),length(moduleList))
    SilhovtFinnMatrix_i <- matrix(0,length(IndxModule_i),2)
    cat("Cluster: ", i,". ")
    cat("Cluster size: ", length(IndxModule_i))
    #Calculate the average inner cluster dissimilarity 
    innerClusteri_simila <- rowSums(1-W1[IndxModule_i,IndxModule_i])/length(IndxModule_i)
    cat("\n----------Finish calculating INNER cluster dissimilarity------------\n")
    SilhovtMatrix_i[,1] <- innerClusteri_simila
    #Calculate the average external cluster dissimilarity from other cluster 
    moduleList <- moduleList[-i]
    bValue <- matrix(0,length(IndxModule_i),1)
    SilhovtFinnMatrix_i[,1] <- SilhovtMatrix_i[,1]
    for (j in 1:length(moduleList)) {
      otherModule <- GeneNames[which(ClusterColors==moduleList[j])]
      Indx_otherModule <- match(otherModule,rownames(W1))
      externalClusteri_simila<- rowSums(1-W1[IndxModule_i,Indx_otherModule])/length(IndxModule_i)
      SilhovtMatrix_i[,j+1] <- externalClusteri_simila
    }
    cat("----------Finish calculating EXTERNAL cluster dissimilarity------------\n")
    bValue <- apply(SilhovtMatrix_i[,-1], 1, FUN=min)
    SilhovtFinnMatrix_i[,2] <- bValue
    maxBandA <- apply(SilhovtFinnMatrix_i, 1, FUN=max)
    # Now calculate the Silhovtte value
    SilhovtCluster_i <- (SilhovtMatrix_i[,1]-bValue)/(maxBandA)
    # cat(length(SilhovtCluster_i))
    cat("----------Finish calculating Silhouvtte value------------\n\n")
    SumarySilhoutM[i] <- sum(SilhovtCluster_i)/length(SilhovtCluster_i)
    #SumarySilhoutM[1:length(SilhovtCluster_i),i] <- SilhovtCluster_i
  }
  names(SumarySilhoutM) <- unique(ClusterColors)
  return(SumarySilhoutM)
}


SilhouettofTuning1 <- SilhovtForModule(dynamicColorsHybrid, W1)
View(SilhouettofTuning1)


# Step4: Visulization (cytoscape) Take the top significant asthmatic module as example
# Step4: Visulization (cytoscape) Take the top significant asthmat --------


BestAsthmaModule <- GeneNamesList1[which(dynamicColorsHybrid=="turquoise")]
length(BestAsthmaModule)
IndexGANFM1 <- match(BestAsthmaModule,rownames(W1))
mod1W1 <- W1[IndexGANFM1,IndexGANFM1]
diag(mod1W1) <- 1
dim(mod1W1)

getwd()
cyt = exportNetworkToCytoscape(mod1W1,
                               edgeFile = paste("/global/project/hpcg1553/Yang/CytoscapeInput-edges-", ".txt", sep=""),
                               nodeFile = paste("/global/project/hpcg1553/Yang/CytoscapeInput-nodes-", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.006,
                               nodeNames = BestAsthmaModule
                               #altNodeNames = modGenes,
                               #nodeAttr = mergedColors[inModule]
)


modGenes1 <- c("AADACL3", "ACPP", "ADCY7","AMIGO3","APOA2","AQP9","C11orf34","C18orf1","C1orf161","C1orf186",
              "C1RL","C21orf7","C6orf105","CCR1","CCRL2","CD5L","CDA","CHI3L1","CLEC4A","CLU","CSTA","CTSG",
              "DEFA4","DHRS9","DPPA4","DYNLL1","EHF","EPX","ESM1","F13A1","FCAR","FCER1A","FCER2","FUT4","FYB",
              "GNLY","GPR21","GPR77","GPR97","HAL","HRH4","IGSF6","IL1R1","IL1RL1","IL31RA","LAIR2","LIMS1",
              "LRG1","MPO","NFE2","P2RX7","P2RY12","PF4","PGLYRP4","PILRA","RBBP5","RNASE3","S100A12","SERPINB2",
              "SLC43A2","SLPI","SPARC","SYNJ2","TM4SF4","TMEM100","TREM2","TRIM14","TSC22D1","TUBA8","TUBB1",
              "VNN3","ZNF124")

IndexTargetGeneSimiMatrix <- match(modGenes1,rownames(W))
modTom1 <- W[IndexTargetGeneSimiMatrix,IndexTargetGeneSimiMatrix]
dim(modTom1)
View(modTom1)

getwd()
cyt = exportNetworkToCytoscape(modTom1,
                               edgeFile = paste("/global/project/hpcg1553/Yang/CytoscapeInput-edges1-", ".txt", sep=""),
                               nodeFile = paste("/global/project/hpcg1553/Yang/CytoscapeInput-nodes1-", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.005,
                               nodeNames = modGenes1
                               #altNodeNames = modGenes,
                               #nodeAttr = mergedColors[inModule]
)



# Stpe5: Comparing with WGCNA expression data ------------------------------------



############ Comparing with WGCNA expression data
objectExpr <- new("ExpressionSet", exprs=as.matrix(Final_autosomalExpr))
objectMeth <- new("ExpressionSet", exprs=as.matrix(Final_autosomalMeth))
filtered_expr<-varFilter(eset = objectExpr,var.func = var,filterByQuantile = TRUE,var.cutoff = 0.75) # taking top 25% variance of expr by variable
filtered_meth<-varFilter(eset = objectMeth,var.func = var,filterByQuantile = TRUE,var.cutoff = 0.75) # taking top 25% variance of methy by variable
Expr_quantile_filtered <- filtered_expr@assayData$exprs
Meth_quantile_filtered <- filtered_meth@assayData$exprs
dim(Expr_quantile_filtered)
dim(Meth_quantile_filtered)

Expr <- t(Expr_quantile_filtered)
Methy <- t(Meth_quantile_filtered)
dim(Expr)
dim(Methy)



###### WGCNA for expression dataset
# * WGCNA for expression dataset --------------------------------------------


# Choose a set of soft-thresholding powers 
powers = c(c(1:20), seq(from = 22, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(Expr, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# We now calculate the adjacencies, using the soft thresholding from above
softPower = 32;
adjacency = adjacency(Expr, power = softPower,type = "signed");

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency,TOMType = "signed");
dissTOM = 1-TOM
rownames(TOM) <- colnames(Expr)
colnames(TOM) <- colnames(Expr)

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = 20);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
length(table(dynamicColors))
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



# Calculate eigengenes
MEList = moduleEigengenes(Expr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
#merge = mergeCloseModules(t(data_filtered_list$Expression), dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergeE <- moduleMergeUsingKME(
  Expr, dynamicColors, ME = MEs, 
  threshPercent = 50, mergePercent = 25, 
  reassignModules = TRUE, 
  convertGrey = TRUE, 
  omitColors = "grey", 
  reassignScale = 1,  
  threshNumber = NULL)

# The merged module colors
MEListKME <- moduleEigengenes(Expr, mergeE$moduleColors)
#MEsKME = MEListKME$eigengenes
mergedColorsE = mergeE$moduleColors;
# Eigengenes of the new merged modules:
mergedMEsE = MEListKME$eigengenes;
length(unique(mergedColorsE))
### relating modules to trait ###

# Define numbers of genes and samples; rows are sample, columns are gene.
nGenes = ncol(Expr);
nSamples = nrow(Expr);


moduleTraitCor = cor(mergedMEsE, phenotype, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


# 
# moduleTraitCor <- as.data.frame(moduleTraitCor)
# moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
# indexnegE <- which(moduleTraitCor[,1]<= -0.1, arr.ind = T)
# indexposE <- which(moduleTraitCor[,1]>= 0.1, arr.ind = T)
# indexcomb <- c(indexnegE,indexposE)
# indexcomb <- sort(indexcomb)
# length(indexcomb)
# moduleCortestE <- moduleTraitCor[indexcomb,]
# modulePvaltestE <- moduleTraitPvalue[indexcomb,]
# 
# interested phenotypes
# a <- c(1,5,6,10)
# 
# phenotype <- as.data.frame(phenotype)
# pheno1 <- phenotype[,a]
# 
# moduleCortestE <- as.matrix(moduleCortestE)
# modulePvaltestE <- as.matrix(modulePvaltestE)
# moduleCortestE <- moduleCortestE[,a]
# modulePvaltestE <- modulePvaltestE[,a]
# 
# sizeGrWindow(10,6)
# # Will display correlations and their p-values
# textMatrix1 = paste(signif(moduleCortestE, 2), "\n(",
#                    signif(modulePvaltestE, 1), ")", sep = "");
# dim(textMatrix1) = dim(moduleCortestE)
# par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleCortestE,
               xLabels = colnames(pheno1),
               yLabels = rownames(moduleCortestE),
               ySymbols = rownames(moduleCortestE),
               cex.lab.y = 0.4,
               colorLabels = FALSE,
               cex.lab.x = 1,
               xLabelsAdj = 0.3,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix1,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Significant module-trait relationships of expression data"))


dim(mergedMEs)
length(unique(mergedColors))
GeneNamesList <- colnames(Expr)
BestAsthmaModule <- GeneNamesList[which(mergedColors=="darkturquoise")]
length(BestAsthmaModule)


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 1), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(phenotype),
               yLabels = rownames(moduleTraitCor),
               ySymbols = names(moduleTraitCor),
               colorLabels = FALSE,
               cex.lab.x = 0.7,
               xLabelsAdj = 0.3,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships of expression data"))
#oldMTPV <- moduleTraitPvalue
colnames(moduleTraitPvalue) <- c("Chemoresistance")
Moduletrait_adjusted_pval_chemosensitivity<-as.data.frame(p.adjust(p = moduleTraitPvalue[,"Chemoresistance"], method="fdr"))
#Moduletrait_adjusted_pval_daystodeath<-as.data.frame(p.adjust(p = moduleTraitPvalue[,"DaystoDeath"], method="bonferroni"))
#Moduletrait_adjusted_pval_stage<-as.data.frame(p.adjust(p = moduleTraitPvalue[,"Stage"], method="bonferroni"))


#hubgene
hub_probes<-as.data.frame(chooseTopHubInEachModule(datExpr = Expr_quantile_filtered,colorh = mergedColors,type="signed"))



###### WGCNA for methylation dataset
# * WGCNA for methylation dataset -------------------------------------------


# Choose a set of soft-thresholding powers 
powers = c(c(1:20), seq(from = 22, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(Methy, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# We now calculate the adjacencies, using the soft thresholding from above (9 or 10 at 50%)
softPower = 10;
adjacency1 = adjacency(Methy, power = softPower,type = "signed");

# Turn adjacency into topological overlap
TOM1 = TOMsimilarity(adjacency1,TOMType = "signed");
dissTOM1 = 1-TOM1
rownames(TOM1) <- colnames(Methy)
colnames(TOM1) <- colnames(Methy)

# Call the hierarchical clustering function
geneTree1 = hclust(as.dist(dissTOM1), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree1, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree1, distM = dissTOM1,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = 20);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors1 = labels2colors(dynamicMods)
length(table(dynamicColors1))
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree1, dynamicColors1, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



# Calculate eigengenes
MEList = moduleEigengenes(Methy, colors = dynamicColors1)
MEs1 = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss1 = 1-cor(MEs1);
# Cluster module eigengenes
METree1 = hclust(as.dist(MEDiss1), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree1, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
#merge = mergeCloseModules(t(data_filtered_list$Expression), dynamicColors, cutHeight = MEDissThres, verbose = 3)
merge <- moduleMergeUsingKME(
  Methy, dynamicColors1, ME = MEs1, 
  threshPercent = 50, mergePercent = 25, 
  reassignModules = TRUE, 
  convertGrey = TRUE, 
  omitColors = "grey", 
  reassignScale = 1,  
  threshNumber = NULL)

# The merged module colors
MEListKME <- moduleEigengenes(Methy, merge$moduleColors)
#MEsKME = MEListKME$eigengenes
mergedColors = merge$moduleColors;
# Eigengenes of the new merged modules:
mergedMEs = MEListKME$eigengenes;

### relating modules to trait ###

# Define numbers of genes and samples; rows are sample, columns are gene.
nGenes = ncol(Methy);
nSamples = nrow(Methy);


moduleTraitCorM = cor(mergedMEs, phenotype, use = "p");
moduleTraitPvalueM = corPvalueStudent(moduleTraitCorM, nSamples);



#### Cytoscape
# Cytoscape ---------------------------------------------------------------


dim(mergedMEsE)
length(unique(mergedColorsE))
GeneNamesList <- colnames(Expr)
BestAsthmaModuleE <- GeneNamesList[which(mergedColorsE=="darkred")]
length(BestAsthmaModuleE)
modGenesEm2 <- BestAsthmaModuleE
View(BestAsthmaModuleE)

IndexTargetGeneSimiMatrixEm2 <- match(modGenesEm2,rownames(TOM))
mod2TomE <- TOM[IndexTargetGeneSimiMatrixEm2,IndexTargetGeneSimiMatrixEm2]
dim(mod2TomE)
View(mod2TomE)


GeneNamesListM <- colnames(Methy)
BestAsthmaModuleM <- GeneNamesListM[which(dynamicColors1=="brown")]
length(BestAsthmaModuleM)
View(BestAsthmaModuleM)
modGenesMm1 <- BestAsthmaModuleM

IndexTargetGeneSimiMatrixMm1 <- match(modGenesMm1,rownames(TOM1))
mod1TomM <- TOM1[IndexTargetGeneSimiMatrixMm1,IndexTargetGeneSimiMatrixMm1]
dim(mod1TomM)
View(mod1TomM)


getwd()
cyt = exportNetworkToCytoscape(mod1TomM,
                               edgeFile = paste("/global/project/hpcg1553/Yang/CytoscapeInput-edges1-Mm1", ".txt", sep=""),
                               nodeFile = paste("/global/project/hpcg1553/Yang/CytoscapeInput-nodes1-Mm1", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.01,
                               nodeNames = rownames(mod1TomM)
                               #altNodeNames = modGenes,
                               #nodeAttr = mergedColors[inModule]
)








moduleTraitCorM <- moduleTraitCorM[,a]
moduleTraitPvalueM <- moduleTraitPvalueM[,a]

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrixM = paste(signif(moduleTraitCorM, 2), "\n(",
                   signif(moduleTraitPvalueM, 1), ")", sep = "");
dim(textMatrixM) = dim(moduleTraitCorM)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCorM,
               xLabels = colnames(pheno1),
               yLabels = rownames(moduleTraitCorM),
               ySymbols = rownames(moduleTraitCorM),
               colorLabels = FALSE,
               cex.lab.y = 0.5,
               colors = blueWhiteRed(50),
               textMatrix = textMatrixM,
               setStdMargins = FALSE,
               xLabelsAdj = 0.3,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Significant module-trait relationships of methylation data"))
#oldMTPV <- moduleTraitPvalue
colnames(moduleTraitPvalue) <- c("Chemoresistance")
Moduletrait_adjusted_pval_chemosensitivity<-as.data.frame(p.adjust(p = moduleTraitPvalue[,"Chemoresistance"], method="fdr"))
#Moduletrait_adjusted_pval_daystodeath<-as.data.frame(p.adjust(p = moduleTraitPvalue[,"DaystoDeath"], method="bonferroni"))
#Moduletrait_adjusted_pval_stage<-as.data.frame(p.adjust(p = moduleTraitPvalue[,"Stage"], method="bonferroni"))








# Step6: Comparison results between GANF and  WGCNA --------


# * Step6.1: Barplot of significant modules from GANF and WGCNA -----------
## Barplot of correlation and pvalues for significant modules from GANFM and WGCNA respectively.
BarPlotTable <- cbind(rownames(as.data.frame(Pvalue))[which(Pvalue<=0.05, arr.ind = T)],Rvalue[which(Pvalue<=0.05, arr.ind = T)],Pvalue[which(Pvalue<=0.05, arr.ind = T)])
rownames(BarPlotTable) <- c()
colnames(BarPlotTable) <- c("ModuleName","Correlation","Pvalue")
BarPlotTable <- as.data.frame(BarPlotTable)
BarPlotTable$Correlation <- round(as.numeric(BarPlotTable$Correlation), digits = 2)
BarPlotTable$Pvalue <- round(as.numeric(BarPlotTable$Pvalue), digits = 3)

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
rownames(moduleTraitPvalue)[which(moduleTraitPvalue$asthma<=0.05, arr.ind = T)]
BarTE <- cbind(c("bisque4","brown4","coral1"),round(moduleTraitCor[which(moduleTraitPvalue$asthma<=0.05, arr.ind = T),1],digits = 2),round(moduleTraitPvalue[which(moduleTraitPvalue$asthma<=0.05, arr.ind = T),1],digits = 3))
rownames(BarTE) <- c()
colnames(BarTE) <- c("ModuleName","Correlation","Pvalue")

moduleTraitPvalueM <- as.data.frame(moduleTraitPvalueM)
rownames(moduleTraitPvalueM)[which(moduleTraitPvalueM$asthma<=0.05, arr.ind = T)]
BarTM <- cbind(c("brown"),abs(round(moduleTraitCorM[which(moduleTraitPvalueM$asthma<=0.05, arr.ind = T),1],digits = 2)),round(moduleTraitPvalueM[which(moduleTraitPvalueM$asthma<=0.05, arr.ind = T),1],digits = 3))
colnames(BarTM) <- c("ModuleName","Correlation","Pvalue")


BarPlotTable <- rbind(BarPlotTable,BarTE,BarTM)
BarPlotTable$Model <- c(rep(c("GANFM"),times=length(which(Pvalue<=0.05,arr.ind = T))),"WGCNA_expr","WGCNA_expr","WGCNA_expr","WGCNA_methy") #times value is decided by sig modules from GANF

c1 <- paste(BarPlotTable$Correlation,"(",BarPlotTable$Pvalue,")",sep = "")
XaxisPosition = paste0(seq_along(BarPlotTable$ModuleName), BarPlotTable$ModuleName)
rowNum <- 1:nrow(BarPlotTable)
BarPlotTable <- cbind(rowNum,XaxisPosition,BarPlotTable)

dev.off()
ggplot(data=BarPlotTable, aes(x=reorder(XaxisPosition, rowNum), y=Correlation, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_discrete(labels=BarPlotTable$ModuleName)+
  geom_text(aes(label=c1), vjust=-0.3, size=5, position = position_dodge(1))+
  theme_minimal()   +
  ggtitle("Significant modules captured by GANFM vs WGCNA") +
  labs(x="Module names", y="Correlation (Pvalue)")+
  theme(axis.text.y =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.6, size = 12),
        plot.title = element_text(hjust = 0.5, size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size=14),
        legend.position="top"
        )


# * Step6.2: Silhouette Score comparison ----------------------------------
## Silhouetto Score for WGCNA_expression (dynamicColors, TOM) and WGCNA_methylation (dynamicColors1, TOM1)
# Expression good modules c("bisque4","brown4","coral1")
SilhouettofTuning1E <- SilhovtForModule(dynamicColors, TOM)
View(SilhouettofTuning1E)
GoodModulesE <- c("bisque4","brown4","coral1")
match(GoodModulesE,names(SilhouettofTuning1E))
SilhExp <- as.data.frame(SilhouettofTuning1E[match(GoodModulesE,names(SilhouettofTuning1E))])
SilhExp$Model <- c(rep(c("WGCNA_expr"),times=nrow(SilhExp)))
colnames(SilhExp)[1] <- c("SilhouettoScore")

# Methylation good modules c("brown")
SilhouettofTuning1M <- SilhovtForModule(dynamicColors1, TOM1)
View(SilhouettofTuning1M)
match("brown",names(SilhouettofTuning1M))
SilhMeth <- as.data.frame(SilhouettofTuning1M[match("brown",names(SilhouettofTuning1M))])
SilhMeth$Model <- c(rep(c("WGCNA_methy"),times=nrow(SilhMeth)))
colnames(SilhMeth)[1] <- c("SilhouettoScore")

# GANF good modules 
SilhouettofTuning1 <- SilhovtForModule(dynamicColorsHybrid, W1)
View(SilhouettofTuning1)
GoodModulesGanf <- BarPlotTable$ModuleName[1:13]
HotsID <- match(GoodModulesGanf,names(SilhouettofTuning1))
HotsGanf <- c(rep(c("cool"),times=length(SilhouettofTuning1)))
SilTableGanf <- cbind(as.data.frame(SilhouettofTuning1),HotsGanf)
SilTableGanf$HotsGanf[HotsID] <- c(rep(c("hot"),times=length(HotsID)))
colnames(SilTableGanf) <- c("SilhScore","ModuleType")

SilhGanf <- as_data_frame(SilhouettofTuning1[HotsID])
SilhGanf$Model <- c(rep(c("GANFM"),times=nrow(SilhGanf)))
rownames(SilhGanf) <- GoodModulesGanf
colnames(SilhGanf)[1] <- c("SilhouettoScore")


# density plot for silhouetto score of each module from ganf model
dev.off()
ggplot(SilTableGanf, aes(x=SilhScore, color=ModuleType)) +
  geom_density()+
  scale_color_manual(values=c("blue", "red"))+
  theme_test()+
  ggtitle("Density Plot of Silhouetto Score for each modules by GANF") +
  labs(x="Silhouetto Score")+
  theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16)
        )


# merge all three tables  and compared with their Silhouetto Scores
SilhScoreTable <- rbind(SilhGanf,SilhExp,SilhMeth)
SilhScoreTable_final <- cbind(SilhScoreTable$SilhouettoScore,BarPlotTable)
colnames(SilhScoreTable_final)[1] <- c("SilhouettoScore")


c2=round(SilhScoreTable_final$SilhouettoScore,digits = 3)

ganfAvg <- sum(SilhScoreTable_final$SilhouettoScore[1:13])/13
wgcnaAvg <- sum(SilhScoreTable_final$SilhouettoScore[14:17]/4)
dev.off()
ggplot(SilhScoreTable_final,aes(x=reorder(XaxisPosition, rowNum), y=SilhouettoScore, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_discrete(labels=SilhScoreTable_final$ModuleName)+
  geom_text(aes(label=c2), vjust=-0.3, size=5, position = position_dodge(1))+
  theme_minimal()   +
  ggtitle("Silhouette Score calculated for significant modules from GANFM vs WGCNA") +
  labs(x="Module names", y="Silhouetto Score")+
  theme(axis.text.y =element_blank(),
        axis.ticks.y =element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.6, size = 12),
        plot.title = element_text(hjust = 0.5, size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size=14),
        legend.position="top"
  )+
  geom_hline(yintercept=ganfAvg, color="red")+
  geom_hline(yintercept=wgcnaAvg, color="black")



#+ coord_cartesian(ylim=c(0, 17))
#require(gridExtra)
#grid.arrange(plot1, plot2, ncol=2)












