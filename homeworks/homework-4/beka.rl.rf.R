library(randomForest)
library(RecordLinkage)
data(RLdata500)
source("C:/Users/sam/Google Drive/stats/code/calc.pcs.R")

true = identity.RLdata500
dat = cbind(RLdata500[,-c(2,4)], true)
dtf = calc.pcs(dat, type = c(rep("j", ncol(dat)-1), "e"))

rf = randomForest(true.comp~., data = dtf[,1:6])
probs = predict(rf, dtf)
preds = 1*(probs > 0.5)

table(dtf$true.comp, preds)

###  Results
#  124700 true negative matches
#  3 false negative matches
#  0 false positive matches
#  47 true positive matches
