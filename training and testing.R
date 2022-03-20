##encoding methods
source("PSNP.R")

library(BSgenome)
library(dplyr)
library(e1071)
library(dplyr)
library(caret)
library(ROCR)
library(xgboost)
library(plyr)

###training data from supplementary material, YBX1 as example
Pdata <- readRNAStringSet("/FA/full_transcript_positive_trainYBX1_35")
Ndata <- readRNAStringSet("/FA/full_transcript_negative_trainYBX1_35")

Pdata_test <- readRNAStringSet("/FA/full_transcript_positive_testYBX1_35")
Ndata_test <- readRNAStringSet("/FA/full_transcript_negative_testYBX1_35")



FTT <- PSNP(as.data.frame(Pdata)$x,as.data.frame(Ndata)$x,NI=2)
Pdata_F <- FTT[1:length(Pdata),]
Ndata_F <- FTT[-c(1:length(Pdata)),]

FTT_test <- PSNP(as.data.frame(Pdata_test)$x,as.data.frame(Ndata_test)$x,NI=2)
Pdata_F_test <- FTT_test[1:length(Pdata_test),]
Ndata_F_test <- FTT_test[-c(1:length(Pdata_test)),]

Train_nunber <- length(Pdata)
test_nunber <- length(Pdata_test)

label_train <- c(rep("motif",Train_nunber),rep("non",Train_nunber))
label_train <- factor(label_train,labels=c("motif","non"))
testlabel <- c(rep(1,test_nunber),rep(0,test_nunber))

BIOmotifvsnontrain <- rbind(Pdata_F,Ndata_F)
BIOmotifvsnontest <- rbind(Pdata_F_test,Ndata_F_test)

BIOmotifvsnon_SVM <- svm(BIOmotifvsnontrain,label_train,cross=5, probability = TRUE)
pred <- predict(BIOmotifvsnon_SVM,BIOmotifvsnontest, probability = TRUE)
motif <- attr(pred, 'probabilities')[,1]
testppred <- prediction(motif,testlabel)
testpppred_auc <- performance(testppred,"auc")

#####performance
attr(testpppred_auc,"y.values")[[1]]

##########used input FA
# user_test <- readRNAStringSet("XXXXXXXXXXXXXXXXX.fa")
# FTT <- PSNP(as.data.frame(user_test)$x,as.data.frame(Ndata)$x,NI=2)
# user_test_F <- FTT[1:length(user_test),]
# user_test_pred <- predict(BIOmotifvsnon_SVM,user_test_F, probability = TRUE)
# user_test_result <- attr(user_test_pred, 'probabilities')[,1]
