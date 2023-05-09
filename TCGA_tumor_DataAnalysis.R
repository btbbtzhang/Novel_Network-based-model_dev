#source('https://bioconductor.org/biocLite.R')
#biocLite('Bioconductor/GenomicDataCommons')
#devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")

library(data.table)
library(readr)
library(TCGAbiolinks)
library(DO.db)
library(affy)
library(qqman)

data1 = read.table("/global/project/hpcg1553/Yang/ESNA/Shrey/FinalCode/GWAS/outputs/gwas_recurrent_wheeze.assoc.logistic", header = TRUE)

getwd()
setwd("/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/")

#import data
#GBM <- fread("/global/home/hpc4009/TCGA_SNF/BREAST/GLIO_Survival.txt")
#GBMexpresion <- read.table("/global/home/hpc4009/TCGA_SNF/GBM/GLIO_Gene_Expression.txt")
#GBMmethylation <- read.table("/global/home/hpc4009/TCGA_SNF/GBM/GLIO_Methy_Expression.txt")
#GBMresistant_barcodes <- fread("/global/home/hpc4009/TCGA_SNF/GBM/GLIO_Survival.txt")

#Breast Cancer data
BreastExpression <- read.table("/global/home/hpc4009/TCGA_SNF/breast_SNP/Breast/BREAST_Gene_Expression.txt")
BreastMethylation <- read.table("/global/home/hpc4009/TCGA_SNF/breast_SNP/Breast/BREAST_Methy_Expression.txt")
BreastBarcodes <- fread("/global/home/hpc4009/TCGA_SNF/breast_SNP/Breast/BREAST_Survival.txt")
Token <- fread("/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/gdc-user-token.2018-08-22T14_39_44-04_00.txt")
Annotationsnp <- fread("/global/home/hpc4009/TCGA_SNF/breast_SNP/GenomeWideSNP_6-na35-annot-csv/GenomeWideSNP_6.na35.annot.csv")


BreBarcodes <- BreastBarcodes$PatientID
BreBarcodes <- substr(BreBarcodes,1,12)

#resistant_string <- lapply(resistant_barcodes,as.character)
#sensitive_barcodes <- read.csv("/global/home/hpc3901/cancerDownload/OV/sensitive_barcodes.txt", col_names = FALSE)
#sensitive_string <- lapply(sensitive_barcodes,as.character)
#TCGAquery_SampleTypes(c('TCGA-A1-A0SD', 'TCGA-A2-A04N'),"TP")
QuerySnp1 <- GDCquery(project = "TCGA-BRCA",
                     data.category = "Simple nucleotide variation",
                     data.type = "Genotypes",
                     barcode = BreBarcodes,
                     experimental.strategy = "Genotyping array",
                     legacy = T,
                     platform = "Affymetrix SNP Array 6.0",
                     sample.type = c("Primary solid Tumor")
)


QuerySnp <- GDCquery(project = "TCGA-BRCA",
                     data.category = "Simple nucleotide variation",
                     data.type = "Simple nucleotide variation",
                     barcode = BreBarcodes,
                     experimental.strategy = "DNA-Seq",
                     legacy = T,
                     platform = c("Illumina HiSeq", "Illumina GA"),
                     sample.type = c("Primary solid Tumor")
                     )

#Test if it has different cell types
#barcode <- Breastresistant_barcodes$PatientID
#Sample <- TCGAquery_SampleTypes(barcode,"TP")



# Extracting manifest files
preManifest <- QuerySnp1$results[[1]]
manifest  <- cbind(preManifest$id, preManifest$file_name, preManifest$md5sum, preManifest$file_size, preManifest$state)
colnames(manifest) <- c("id","filename","md5","size","state")
#### paste("a","b","c", sep="-") -> "a-b-c"
getwd()
write.table(manifest, file = "manifest.txt", quote = F, sep = "\t", row.names = F)

#GDCdownload(downloadRnaEx, method = "client", token.file = Token)


######## Load all genotype data into R
tem <- list.files(path = "/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/preData/", recursive = TRUE,pattern = "*.birdseed.data.txt$",full.names = TRUE)
myfiles <- lapply(tem, read.delim)

for (j in 1:length(myfiles)) {
  colnames(myfiles[[j]]) <- c("ProbeName", "Genotype","Confidence")
}

# test if we set up the low confidence row as 0: testlist <- myfiles[[1]]$Confidence[-3]

for ( i in 1:length(myfiles)){
  for(k in 2:nrow(myfiles[[i]])){
    if (myfiles[[i]]$Confidence[k] > 0.5){
      myfiles[[i]]$Genotype[k] <- -1
      }
  }
  print(i) ## Track the progress
}
getwd()
save(myfiles, file = "/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/snpPreData.RData")
load("/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/snpPreData.RData")


#### Merge snpData with annotation file
Snp_Probes <- as.matrix(myfiles[[1]]$ProbeName)
Snp_Probes <- Snp_Probes[-1,]
Snp_Probes <- as.data.frame(Snp_Probes)
colnames(Snp_Probes) <- "Probe Set ID"

#merge list
CombinedList <- merge(Snp_Probes, Annotationsnp, by="Probe Set ID",sort = FALSE)
rm(Snp_Probes)
#take the first four columns of interest
CombinedListDone <- CombinedList[1:4]
Mapfile <- cbind(CombinedListDone[,3],CombinedListDone[,2])
a <- matrix(0,nrow(Mapfile),1)
Mapfile <- cbind(Mapfile,a)
Mapfile <- cbind(Mapfile, CombinedListDone[,4])
colnames(Mapfile) <- c("chromosome", "rs#", "Genetic distance", "Base-pair position")
dim(Mapfile)
View(Mapfile[1:10,])
rm(a)
Mapfile <- as.data.frame(Mapfile)
misinProbesInd <- which(Mapfile$chromosome=="---")
Mapfile <- Mapfile[-misinProbesInd,]
Mapfile <- as.data.frame(Mapfile)

getwd()
######### Save Map file
fwrite(x = Mapfile, file = "plink.map", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE, showProgress = TRUE, nThread = 8, verbose = TRUE)




################################################# Ped file
CombinedCalls <- matrix(0,906601,105)
CombinedCalls[,1] <- myfiles[[1]]$Genotype

for (i in 1:length(myfiles)) {
  CombinedCalls[,i] <- as.numeric(myfiles[[i]]$Genotype) 
}

CombinedCalls <- CombinedCalls[-1,]
colnames(CombinedCalls) <- BreBarcodes
rownames(CombinedCalls) <- CombinedListDone$`dbSNP RS ID`
CombinedCalls <- t(CombinedCalls)
dim(CombinedCalls)
Pedfile <- matrix(0,105,6)
colnames(Pedfile) <- c("Family ID", "Individual ID", "Paternal ID", "Maternal ID", "Sex", "Phenotype")
Pedfile[,1:2] <- BreBarcodes
Pedfile[,3:4] <- matrix(0,105,2)
Pedfile[,5] <- matrix(2,105,1)
Pedfile[,6] <- matrix(-9,105,1)
PedCombined <- cbind(Pedfile,CombinedCalls)
View(PedCombined[,1:10])
rm(CombinedCalls)
PedCombined1 <- PedCombined[,7:ncol(PedCombined)]
PedCombined1 <- PedCombined1[,-misinProbesInd]
dim(PedCombined1)

PedCombined1[PedCombined1==2]<-"BB"
PedCombined1[PedCombined1==1]<-"AB"
PedCombined1[PedCombined1==0]<-"AA"
PedCombined1[PedCombined1==-1]<-"00"
dim(PedCombined1)
PedCombined1 <- cbind(Pedfile,PedCombined1)
View(PedCombined1[,1:10])

PedCombined1<-as.data.frame(PedCombined1)
getwd()
fwrite(x = PedCombined1, file = "plink.ped", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE, showProgress = TRUE, nThread = 8, verbose = TRUE)

######## Load the results from Plink~~~
afterSnp <- fread("/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/plinkData/OutCombined_plink1.raw", header = F)
afterSnp <- afterSnp[-1,]
afterAnno <- fread("/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/plinkData/Out_plink.map")
colnames(afterAnno) <- colnames(Mapfile)

rsNumbers <- as.data.frame(afterAnno$`rs#`)
colnames(rsNumbers) <- "dbSNP RS ID"

##### Extracting rsNumbers from Annotationsnp file
#names(Annotationsnp)


#ind <- matrix(,length(rsNumbers),1)
#for (i in 1:length(rsNumbers)) {
#  ind[i] <- which(Annotationsnp$`dbSNP RS ID`==rsNumbers[i])
#  print(i)
#}

Indexi <- match(afterAnno$`rs#`,Annotationsnp$`dbSNP RS ID`)
length(Indexi)
EnsemblGee <- Annotationsnp[Indexi,]
TranscriptID <- substr(EnsemblGee$`Associated Gene`,1,15)
EnsemblG <- substr(EnsemblGee$`dbSNP RS ID`,1,15)
#EnsemblGee <- as.data.frame(EnsemblGee)
#colnames(EnsemblGee) <- "EnsemblID"
#EnsemblGee <- as.character(EnsemblGee)
#Combine rs# and transcpritID
AsciatedRsTras <- cbind(rsNumbers,as.data.frame(TranscriptID))

#### Making Bed file. 4Cols: Chr, start postion, end position, rsID.
Bedfile <- matrix(0,nrow(EnsemblGee),4)
Bedfile[,1] <- paste("chr",EnsemblGee$Chromosome, sep = "")  
Bedfile[,2] <- EnsemblGee$`Physical Position`
Bedfile[,3] <- Bedfile[,2]
Bedfile[,4] <- EnsemblGee$`dbSNP RS ID`
Bedfile <- as.data.frame(Bedfile)
getwd()
fwrite(Bedfile,file = "YangSnpData.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE, showProgress = TRUE, nThread = 8, verbose = TRUE)
############### getting associating gene info if you have rs# or transcrpit id or whatever
library(biomaRt)
ensembl1 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
TargetGenes1 <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), 
                     filters ='ensembl_transcript_id', values =TranscriptID, mart = ensembl1)


ensembl = useEnsembl(biomart="snp", dataset="hsapiens_snp")
TargetGenes <- getBM(attributes=c('refsnp_id','chr_name','chrom_start','chrom_end','variation_names','associated_gene','ensembl_gene_stable_id','ensembl_transcript_stable_id' ), 
                     filters ='snp_filter', values =EnsemblG, mart = ensembl)
IndexiGe <- match(TargetGenes1$ensembl_transcript_id, AsciatedRsTras$TranscriptID)
AsciatedRsTrasGe <- AsciatedRsTras[IndexiGe,]
colnames(AsciatedRsTrasGe)[2] <- c("EnsemblTranscript")
rm(AsciatedRsTras)
IndexiGeSnp <- match(AsciatedRsTrasGe$`dbSNP RS ID`, afterAnno$`rs#`)
IndexiGeSnp <- as.double(IndexiGeSnp+6)



afterSnp<-as.data.frame(afterSnp)
colnames(afterSnp)[1:6] <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
#colnames(afterSnp)[7] <- paste("rs10458597")
colnames(afterSnp)[7:ncol(afterSnp)] <- afterAnno$`rs#`

View(afterSnp[,1:10])
SnpGeneData <- afterSnp[,IndexiGeSnp]

rownames(SnpGeneData) <- afterSnp$IID
SnpGeneData <- t(SnpGeneData)
View(SnpGeneData[,1:100])
SNPGene <- cbind(TargetGenes1,AsciatedRsTrasGe$`dbSNP RS ID`,SnpGeneData)
View(SNPGene[1:100,])
assGene <- TargetGenes1$hgnc_symbol
inde <- which(assGene=="")
FinSNPGeneT <- SNPGene[-inde,]
View(FinSNPGeneT[1:100,])

colnames(FinSNPGeneT)[7] <- "rsID"
sort_SnpGeneT <- FinSNPGeneT[order(FinSNPGeneT$hgnc_symbol),]
View(sort_SnpGeneT[1:100,])
#test first
sort_SnpGeneT <- as.matrix(sort_SnpGeneT)
sort_SnpGeneTOri <- sort_SnpGeneT[,8:ncol(sort_SnpGeneT)]
sort_SnpGeneT <- rbind(sort_SnpGeneT,matrix(0,1,ncol(sort_SnpGeneT))) 
dupGenInd <- vector()
dim(sort_SnpGeneT)
b <- which(is.na(sort_SnpGeneT)==TRUE,arr.ind = TRUE)
naGeneNameInd <- b[1:12,1]
sort_SnpGeneT <- sort_SnpGeneT[-naGeneNameInd,]

#### impute missing genomic value with majority genotype value of this gene
###      
for (i in 1:nrow(sort_SnpGeneTOri)) {
  b <- which(is.na(sort_SnpGeneTOri[i,])==TRUE)
  c <- as.data.frame(table(sort_SnpGeneTOri[i,]))
  d <- which.max(c$Freq)
  sort_SnpGeneTOri[i,b] <- c[d,1]
  if (i %% 1000==0) {
    print(i)
  }
}

sort_SnpGeneT[,8:ncol(sort_SnpGeneT)] <- sort_SnpGeneTOri
## Testing if dataset still has NA values
x <- which(is.na(sort_SnpGeneT)==TRUE,arr.ind = T)


##### Calculation for Mutation weight per gene
for (i in 1:nrow(sort_SnpGeneT)) {
  if (sort_SnpGeneT[i,3]==sort_SnpGeneT[i+1,3]) {
    sort_SnpGeneT[i+1,8:112] <- as.numeric(sort_SnpGeneT[i,8:112])+as.numeric(sort_SnpGeneT[i+1,8:112])
    cat("Duplicated gene:",i,"|")
    dupGenInd[i] <- i
  } else{
    sort_SnpGeneT[i,] <- sort_SnpGeneT[i,]
  }
}
dupGenInd <- dupGenInd[!is.na(dupGenInd)]
FinSort_SnpGeneT <- sort_SnpGeneT[-dupGenInd,]
FinSort_SnpGeneT <- FinSort_SnpGeneT[-nrow(FinSort_SnpGeneT),]

dim(FinSort_SnpGeneT)
dim(sort_SnpGeneTOri)
View(FinSort_SnpGeneT[1:100,])
View(sort_SnpGeneT[1:100,])

#### ExpresAft <- lumiN(RoutSamExpres, method='rsn') # Using Robust Spline Normalization (RSN)
rm(SnpGeneData)
FinSort_SnpGeneT <- as.data.frame(FinSort_SnpGeneT)
GenoData <- cbind(FinSort_SnpGeneT$hgnc_symbol,FinSort_SnpGeneT$start_position,FinSort_SnpGeneT$end_position,FinSort_SnpGeneT[8:112])
GenoData <- cbind(matrix(0,nrow(GenoData),1),GenoData)
colnames(GenoData)[1:4] <- c("Gene Length", "Gene Name", "Start P", "End P")
GenoData <- as.matrix(GenoData)
GenoData[,1] <- as.numeric(GenoData[,4])-as.numeric(GenoData[,3])
GenoData[,5:ncol(GenoData)] <- as.numeric(GenoData[,5:ncol(GenoData)])/as.numeric(GenoData[,1])
dim(GenoData)
View(GenoData[1:100,])
max(as.numeric(GenoData[,5:ncol(GenoData)]))

StatAll <- function(x){
  MeanDa <- colMeans(x)
  SDDa <- apply(x, 2, sd)
  MmeanDa <- matrix(rep(apply(x, FUN = mean, 2),nrow(x)),nrow=nrow(x),byrow=T)
  PredistDa <- abs(x-MmeanDa)
  DistDa <- colSums(PredistDa)/nrow(PredistDa)
  StatAll <- cbind(MeanDa, SDDa, DistDa)
  return(StatAll)
}
GeeNumeric <- matrix(as.numeric(GenoData[,5:ncol(GenoData)]),nrow(GenoData),105)
rownames(GeeNumeric) <- GenoData[,2]
colnames(GeeNumeric) <- colnames(GenoData[,5:ncol(GenoData)])
GeeStats <- StatAll(t(GeeNumeric))
View(GeeNumeric[1:100,])
#GeeNumeric <- GeeNumeric+1000
dim(GeeNumeric)

fwrite(as.data.frame(GenoData), file = "/global/home/hpc4009/TCGA_SNF/breast_SNP/Breast/GenoData_continuous.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, verbose = TRUE)
Genotype <- fread("/global/home/hpc4009/TCGA_SNF/breast_SNP/Breast/GenoData_continuous.txt")


library(IlluminaHumanMethylation27k.db)
x <- IlluminaHumanMethylation27kSYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
MethyAnno <- as.list(x[mapped_probes])
MethyAnno <- as.data.frame(MethyAnno)
BreastMethyAnno <- as.matrix(cbind(rownames(MethyAnno),as.character(MethyAnno$V1)))
View(BreastMethyAnno)
BreastMethylation <- as.matrix(BreastMethylation)
IndeMethy <- match(rownames(BreastMethylation),BreastMethyAnno[,1])
BreastMethyAnno <- BreastMethyAnno[IndeMethy,]
View(BreastMethyAnno)
rownames(BreastMethylation) <- BreastMethyAnno[,2]
View(BreastMethylation[1:100,1:100])
colnames(BreastMethylation) <- colnames(GenoData[,5:ncol(GenoData)])
colnames(BreastExpression) <- colnames(GenoData[,5:ncol(GenoData)])
common <- intersect(rownames(BreastMethylation),rownames(GeeNumeric))
common <- intersect(common,rownames(BreastExpression))
commonIndEx <- match(common,rownames(BreastExpression))
commonIndMe <- match(common,rownames(BreastMethylation))
commonIndGe <- match(common,rownames(GeeNumeric))
############ Final Datasets of expression, methylation, genotype
FinalBreastExpre <- BreastExpression[commonIndEx,]
FinalBreastMethy <- BreastMethylation[commonIndMe,]
FinalBreastGenot <- GeeNumeric[commonIndGe,]
FinalBreastExpre <- as.data.frame(FinalBreastExpre)
FinalBreastMethy <- as.data.frame(FinalBreastMethy)
FinalBreastGenot <- as.data.frame(FinalBreastGenot)
getwd()
fwrite(FinalBreastExpre, file = "/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/FinalGenesData/FinalBreastExpre_Majo.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, verbose = TRUE)
fwrite(FinalBreastMethy, file = "/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/FinalGenesData/FinalBreastMethy_Majo.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, verbose = TRUE)
fwrite(FinalBreastGenot, file = "/global/home/hpc4009/TCGA_SNF/breast_SNP/snpdata/FinalGenesData/FinalBreastGenot_Majo.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote=FALSE, verbose = TRUE)


x <- methTable[2,]
x <- t(x)
hist(x, breaks=20)

methStats <- StatAll(t(methTable))
variancMeth <- apply(methTable, 1, var)
variancMeth <- as.data.frame(variancMeth)
methStats$Variance <- variancMeth
names(variancMeth) <- c("Variance")
max(variancMeth)
min(variancMeth)
distaceM <- max(variancMeth)-min(variancMeth)
# range 1: (0.005989271,2.757109) 0% ~20%
# range 2: (2.757109,5.508229)  20% ~40%
# range 3: (5.508229,8.259349)  40% ~60%
# range 4: (8.259349,11.01047)  60% ~80%
# range 5: (11.01047,13.76159)  80% ~100%

range1_min <- 0.005989271
range1_max <- 2.757109

range5_min <- 11.01047
range5_max <- 13.76159

length(which(variancMeth>=range1_min & variancMeth<=range1_max, arr.ind=T))
length(which(variancMeth>=range5_min & variancMeth<=range5_max, arr.ind=T))
length(which(variancMeth>=range1_min & variancMeth<=range1_max, arr.ind=T))
length(which(variancMeth>=range1_min & variancMeth<=range1_max, arr.ind=T))
length(which(variancMeth>=range1_min & variancMeth<=range1_max, arr.ind=T))

length(variancMeth)

