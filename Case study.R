library(readr)
library(BSgenome.Hsapiens.UCSC.hg19)

#m5C peak region from meRIP-seq
T1 <- read_tsv("GSE53370_consistent_m5C_sites_in_MCF10A.bed")

T1gr <- GRanges(seqnames = T1$`# chr`,
                 ranges = IRanges(start = T1$chromStart,end=T1$chromEnd)
                 ,strand = T1$strand)


# devtools::install_github("ZhenWei10/m6ALogisticModel")
# obtain cytosins from peak region
m5C <- m6ALogisticModel::sample_sequence("C",T1gr,Hsapiens,Fixed=T)

library(BSgenome)
library(dplyr)
library(e1071)
library(dplyr)
library(caret)
library(ROCR)
library(xgboost)
library(plyr)
source("PSNP.R")
#predict substrate of NSUN2 based on the mature mRNA model
Pdata <- readRNAStringSet("/FA/mature_mRNA_positive_trainNSUN2_30")
Ndata <- readRNAStringSet("/FA/mature_mRNA_negative_trainNSUN2_30")

#feature encoding for training data
FTT <- PSNP(Pdata1,Ndata,NI=2)
Pdata_F <- FTT[1:length(Pdata),]
Ndata_F <- FTT[-c(1:length(Pdata)),]

Train_nunber <- c(1:length(Pdata))[c(1:round(length(Pdata)/10*10))]
label_train <- c(rep("motif",length(Train_nunber)),rep("non",length(Train_nunber)))
label_train <- factor(label_train,labels=c("motif","non"))
training <- as.data.frame(rbind(Pdata_F[Train_nunber,],Ndata_F[Train_nunber,]))

BIOmotifvsnontrain <- rbind(Pdata_F,Ndata_F)

#model training for NSUN2 substrate prediction
BIOmotifvsnon_SVM <- svm(BIOmotifvsnontrain,label_train,cross=5, probability = TRUE)

####RNA sequences with 61bp flanking windows 
data_test <- as.character(DNAStringSet(Views(Hsapiens,m5C+30)))
data_test <- gsub("T","U",data_test)

##feature encoding for test data 
FTT <- PSNP(data_test,Ndata,NI=2)
data_test_F <- FTT[1:length(data_test),]


pred <- predict(BIOmotifvsnon_SVM,Ntest_F, probability = TRUE)
#probability for each cytosin was regulated by NSUN2
attr(pred, 'probabilities')[,1]


#predict substrate of YBX1 based on the mature mRNA model
Pdata <- readRNAStringSet("/FA/mature_mRNA_positive_trainYBX1_35")
Ndata <- readRNAStringSet("/FA/mature_mRNA_negative_trainYBX1_35")

#feature encoding for training data
FTT <- PSNP(Pdata1,Ndata,NI=2)
Pdata_F <- FTT[1:length(Pdata),]
Ndata_F <- FTT[-c(1:length(Pdata)),]

Train_nunber <- c(1:length(Pdata))[c(1:round(length(Pdata)/10*10))]
label_train <- c(rep("motif",length(Train_nunber)),rep("non",length(Train_nunber)))
label_train <- factor(label_train,labels=c("motif","non"))
training <- as.data.frame(rbind(Pdata_F[Train_nunber,],Ndata_F[Train_nunber,]))

BIOmotifvsnontrain <- rbind(Pdata_F,Ndata_F)

#model training for YBX1 substrate prediction
BIOmotifvsnon_SVM <- svm(BIOmotifvsnontrain,label_train,cross=5, probability = TRUE)

####RNA sequences with 71bp flanking windows 
data_test <- as.character(DNAStringSet(Views(Hsapiens,m5C+30)))
data_test <- gsub("T","U",data_test)

##feature encoding for test data 
FTT <- PSNP(data_test,Ndata,NI=2)
data_test_F <- FTT[1:length(data_test),]

pred <- predict(BIOmotifvsnon_SVM,Ntest_F, probability = TRUE)
#probability for each cytosin was regulated by YBX1
attr(pred, 'probabilities')[,1]


#predict substrate of ALYREF based on the mature mRNA model
Pdata <- readRNAStringSet("/FA/mature_mRNA_positive_trainALYREF_25")
Ndata <- readRNAStringSet("/FA/mature_mRNA_negative_trainALYREF_25")

#feature encoding for training data
FTT <- PSNP(Pdata1,Ndata,NI=2)
Pdata_F <- FTT[1:length(Pdata),]
Ndata_F <- FTT[-c(1:length(Pdata)),]

Train_nunber <- c(1:length(Pdata))[c(1:round(length(Pdata)/10*10))]
label_train <- c(rep("motif",length(Train_nunber)),rep("non",length(Train_nunber)))
label_train <- factor(label_train,labels=c("motif","non"))
training <- as.data.frame(rbind(Pdata_F[Train_nunber,],Ndata_F[Train_nunber,]))

BIOmotifvsnontrain <- rbind(Pdata_F,Ndata_F)

#model training for ALYREF substrate prediction
BIOmotifvsnon_SVM <- svm(BIOmotifvsnontrain,label_train,cross=5, probability = TRUE)

####RNA sequences with 71bp flanking windows 
data_test <- as.character(DNAStringSet(Views(Hsapiens,m5C+30)))
data_test <- gsub("T","U",data_test)

##feature encoding for test data 
FTT <- PSNP(data_test,Ndata,NI=2)
data_test_F <- FTT[1:length(data_test),]

pred <- predict(BIOmotifvsnon_SVM,Ntest_F, probability = TRUE)
#probability for each cytosin was regulated by ALYREF
attr(pred, 'probabilities')[,1]