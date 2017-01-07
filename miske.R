## made by blaz skrlj
## data exploration of mouse gene expression dataset





miske <- read.csv("../miske.csv")

## quick insight

summary(miske)

## description of class attribute>

# c-CS-m: control mice, stimulated to learn, injected with memantine (10 mice)
# c-SC-s: control mice, not stimulated to learn, injected with saline (9 mice)
# c-SC-m: control mice, not stimulated to learn, injected with memantine (10 mice)
# 
# t-CS-s: trisomy mice, stimulated to learn, injected with saline (7 mice)
# t-CS-m: trisomy mice, stimulated to learn, injected with memantine (9 mice)
# t-SC-s: trisomy mice, not stimulated to learn, injected with saline (9 mice)
# t-SC-m: trisomy mice, not stimulated to learn, injected with memantine (9 mice) 

## Definition of the problem> 
## main aim is to discriminate between the classes. This can be further splitted to two subproblems,
## in one  we are trying to classify trisomy and in other control mice. So first, we take this
## problem as a two class classificaiton problem. 

## Function definitions
recode2 <- function(val){
  tval <- ""
  if (substr(val,0,1) == "t"){
    tval <- "T"
  }else{
    tval = "C"
  }
  return (tval)
}

recode4 <- function(val){
  tval <- ""
  if (substr(val,0,1) == "t" && substr(val, 3,4) == "CS"){
    tval <- "Ts"
  }else if(substr(val,0,1) == "t" && substr(val, 3,4) == "SC"){
    tval <- "Tn"
  }else if(substr(val,0,1) == "c" && substr(val, 3,4) == "CS"){
    tval <- "Cs"
  }else if(substr(val,0,1) == "c" && substr(val, 3,4) == "SC"){
    tval <- "Cn"
  }else{
    tval <- NA
  }
  return (tval)
}

  
  to_means <- function(df){
    for (j in unique(df$class)){
      for (m in names(which(apply(df,2, function(x) any(is.na(x)))))){
        df[,m]<-sapply(df[,m],
               function(x) ifelse(is.na(x),
                                  x<-mean(df[df$class == j,m]),
                                  x<-x))

      }
    }
    return (df)
    #print(length(na.omit(df)))
  }

to_num <- function(vec){
  return (as.numeric(as.factor(vec)))
}

## create a new attribute, which will contain only two classes, ot 4 classes.
miske$two_class <- sapply(miske$class,function(x) recode2(x))
miske$four_class <- sapply(miske$class,function(x) recode4(x))

table(miske$class)/nrow(miske)
table(miske$two_class)/nrow(miske)
table(miske$four_class)/nrow(miske)

#this confirms attribute creation. 
## We know the distributions, if algorithm should perform better than
## 14%, 53% and 27% for specific classes.

## PREPROCESSING:
## assign numeric attributes for last 6 columns - this will be useful for learning. 
## assign numeric values to last 6 columns>


# miske[, ((ncol(miske)-5):ncol(miske))]<- 
#   sapply(miske[, ((ncol(miske)-5):ncol(miske))], function(x) to_num(x))

## this worked as expected. Moving forward to dealing with missing values:
## first let's see which columns include missing values.

missing_sums <- which(sapply(miske, function(x) sum(is.na(x)))>0)
length(missing_sums) ## here we can see, that 49 genes have missing values.

## From this point on, we have at least two options, first: impute, second: means, third: extrapolate
## we shall test all 3 of those options.

miske_imputed <- na.omit(miske)
miske_means <- to_means(miske) ## this is a more sophisticated approach.
missing_sums2 <- which(sapply(miske_means, function(x) sum(is.na(x)))>0)
length(missing_sums2) ## here we can see, that 49 genes have missing values.

colnames(miske_means) <- sapply(colnames(miske_means), function(x) gsub('_N','',x))
## two classes have missing values, those two couldn't be redefined.
# Exploratory data analysis is conducted with miske_means dataset


to_explore <- miske_means[,1:78]

to_explore %>% melt() %>% group_by(variable) %>% 
  filter(value<1.5) %>% ggplot(aes(value))+geom_density(aes(fill=variable),alpha=0.5)+
  theme(legend.position="bottom")+
  labs(variable='Gene name', value='Expression level', fill='gene name')+
  ggtitle("Gene expression distributions")+xlab("Expression level")

to_explore %>% melt() %>% group_by(variable) %>% 
  filter(value>1.5) %>% ggplot(aes(value))+geom_density(aes(fill=variable),alpha=0.5)+
  theme(legend.position="bottom")+
  labs(variable='Gene name', value='Expression level', fill='gene name')+
  ggtitle("Gene expression distributions")+xlab("Expression level")


to_explore %>% melt() %>% group_by(variable) %>%
  summarise(v_mean = mean(value), v_sd = sd(value))%>% select(v_mean) -> mns

qns <- quantile(mns$v_mean, na.rm=T)
## 83k rows
library(data.table)
library(dplyr)
library(ggplot2)
to_explore %>% melt() %>% group_by(variable) %>%
summarise(v_mean = mean(value), v_sd = sd(value)) %>%
arrange(v_mean) %>% 
ggplot(aes(variable, v_mean))+geom_bar(stat="identity", aes(fill=v_mean))+ 
  geom_point(aes(size=v_sd^2, color=v_sd))+
  geom_hline(yintercept=as.numeric(quantile(mns$v_mean, na.rm=T)[1]))+
  geom_hline(yintercept=as.numeric(quantile(mns$v_mean, na.rm=T)[2]))+
  geom_hline(yintercept=as.numeric(quantile(mns$v_mean, na.rm=T)[3]))+
  geom_hline(yintercept=as.numeric(quantile(mns$v_mean, na.rm=T)[4]))+
  geom_hline(yintercept=as.numeric(quantile(mns$v_mean, na.rm=T)[5]))+
  geom_text(aes(0,as.numeric(quantile(mns$v_mean, na.rm=T)[1]),label = '0%', hjust = -1, vjust=+1), color='black')+
  geom_text(aes(0,as.numeric(quantile(mns$v_mean, na.rm=T)[2]),label = '25%', hjust = -1),color='black')+
  geom_text(aes(0,as.numeric(quantile(mns$v_mean, na.rm=T)[3]),label = '50%', hjust = -1),color='black')+
  geom_text(aes(0,as.numeric(quantile(mns$v_mean, na.rm=T)[4]),label = '75%', hjust = -1),color='black')+
  geom_text(aes(0,as.numeric(quantile(mns$v_mean, na.rm=T)[5]),label = '100%', hjust=-1),color='black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  scale_colour_gradient(limits=c(0, 3), low="gray", high="yellow")+
  labs(fill='mean expression', 
       color='Variance of measurements', 
       size='deviation',
       title='Visualization of expression dataset',
       xlab='Gene',
       ylab='Expression level')+xlab('Gene name')+ylab('Expression mean')


## genes with the largest sd
to_explore %>% melt() %>% group_by(variable, std = sd(value)) %>% 
  filter(std > 0.78 && value>2) %>% 
  ggplot(aes(variable, value))+geom_boxplot(aes(fill=variable))+
  geom_hline(yintercept=2)+
  labs(title='Deviation of 8 most variable genes',
       xlab='',ylab='Expression')+
  theme(axis.title.x=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.ticks.x=element_blank())
  
  
to_explore$class <- miske_means$class

to_explore[,2:ncol(to_explore)] %>% melt(id.var = c('class')) %>% 
  group_by(class,variable) %>% summarise(m_val = mean(value)) %>% 
  ggplot(aes(variable, m_val))+
  geom_bar(stat="identity", aes(fill=class))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(fill='target attribute', 
       color='Target class',
       title='Visualization of expression dataset, vol.2',
       variable='Gene',
       m_val='Expression level')+ylab('Expression')+xlab('Gene name')

to_explore %>% select(class) %>% melt() %>% group_by(class) %>% summarise(sums = n()) %>%
  ggplot(aes(class,reorder(sums, -sums)))+geom_bar(stat="identity")+ 
  ylab('Count of measurements')+
  xlab('Target class')+ggtitle('Counts of target variable values')

############## T TESTS

t.test(to_explore$pCAMKII_N, to_explore$Ubiquitin_N)
library('scatterplot3d')

scatterplot3d(to_explore$pCAMKII_N,to_explore$Ubiquitin_N, to_explore$NR1_N, pch=16, highlight.3d=TRUE,
              type="h", main="3D Scatterplot of gene expression") 

###

## clustering analysis:

clusters <- hclust(dist(to_explore))
plot(clusters)
### small to big expression comparisons:

ggplot(data= to_explore, aes(pCAMKII_N, Ubiquitin_N))+geom_density_2d()+ggtitle("Comparison of two expression levels")
## MACHINE LEARNING PART - 2:78 are expression levels used for prediction.
library(caret)
library(mlbench)
library(doMC)
registerDoMC(cores=2)
control <- trainControl(method="repeatedcv", number=10, repeats=2)
seed <- 7
metric <- "Accuracy"
# set.seed(seed)
# fit.lda <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="lda", metric=metric, preProc=c("center", "scale"), trControl=control)
# # Logistic Regression
# set.seed(seed)
# fit.glm <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="glm", metric=metric, trControl=control)
# # GLMNET
# set.seed(seed)
# fit.glmnet <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="glmnet", metric=metric, preProc=c("center", "scale"), trControl=control)
# # SVM Radial
# set.seed(seed)
# fit.svmRadial <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="svmRadial", metric=metric, preProc=c("center", "scale"), trControl=control, fit=FALSE)
# # kNN
# set.seed(seed)
# fit.knn <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="knn", metric=metric, preProc=c("center", "scale"), trControl=control)
# # Naive Bayes
# set.seed(seed)
# fit.nb <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="nb", metric=metric, trControl=control)
# # CART
# set.seed(seed)
# fit.cart <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="rpart", metric=metric, trControl=control)
# # C5.0
# set.seed(seed)
# fit.c50 <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="C5.0", metric=metric, trControl=control)
# # Bagged CART
# set.seed(seed)
# fit.treebag <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="treebag", metric=metric, trControl=control)
# # Random Forest
# set.seed(seed)
# fit.rf <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="rf", metric=metric, trControl=control)
# # Stochastic Gradient Boosting (Generalized Boosted Modeling)
# set.seed(seed)
#fit.gbm <- train(as.factor(class)~., data=na.omit(miske[,1:82]), method="gbm", metric=metric, trControl=control, verbose=FALSE)


# results <- resamples(list(lda=fit.lda, logistic=fit.glm, glmnet=fit.glmnet,
#                           svm=fit.svmRadial, knn=fit.knn, nb=fit.nb, cart=fit.cart, c50=fit.c50,
#                           bagging=fit.treebag, rf=fit.rf, gbm=fit.gbm))

# Table comparison
train_set <- na.omit(miske_means[,-(79:81)]) ## with classes
fit.rf <- train(as.factor(class)~., data=train_set, method="rf", metric=metric, trControl=control)
set.seed(seed)
fit.gbm <-train(as.factor(class)~., data=train_set, method="gbm", metric=metric, trControl=control, verbose=FALSE)
set.seed(seed)
fit.nb <- train(as.factor(class)~., data=train_set, method="nb", metric=metric, trControl=control)
# # CART
set.seed(seed)
fit.svmRadial <- train(as.factor(class)~., data=train_set, method="svmRadial", metric=metric, preProc=c("center", "scale"), trControl=control, fit=FALSE)
set.seed(seed)
fit.c50 <- train(as.factor(class)~., data=train_set, method="C5.0", metric=metric, trControl=control)
set.seed(seed)
fit.cart <-train(as.factor(class)~., data=train_set, method="rpart", metric=metric, trControl=control)


# set.seed(seed)
fit.treebag <- train(as.factor(class)~., data=train_set, method="treebag", metric=metric, trControl=control)

results <- resamples(list( rf=fit.rf, 
                           gbm=fit.gbm, 
                           nb=fit.nb, 
                           svm=fit.svmRadial, 
                           c50=fit.c50,
                           tbg=fit.treebag,
                           cart = fit.cart))
summary(results)
# summarize the distributions
summary(results)
# boxplots of results
bwplot(results)
# dot plots of results
dotplot(results)

## NO ADDITIONAL VARIABLES!
train_set2 <- na.omit(miske_means[,-(79:81)]) ## with classes
train_set2 <- train_set2[,-c(80,81)] ## with classes
fit.rf1 <- train(as.factor(class)~., data=train_set2, method="rf", metric=metric, trControl=control)
set.seed(seed)
fit.gbm1 <-train(as.factor(class)~., data=train_set2, method="gbm", metric=metric, trControl=control, verbose=FALSE)
set.seed(seed)
fit.nb1 <- train(as.factor(class)~., data=train_set2, method="nb", metric=metric, trControl=control)
# # CART
set.seed(seed)
fit.svmRadial1 <- train(as.factor(class)~., data=train_set2, method="svmRadial", metric=metric, preProc=c("center", "scale"), trControl=control, fit=FALSE)
set.seed(seed)
fit.c501 <- train(as.factor(class)~., data=train_set2, method="C5.0", metric=metric, trControl=control)
set.seed(seed)
fit.cart1 <-train(as.factor(class)~., data=train_set2, method="rpart", metric=metric, trControl=control)


# set.seed(seed)
fit.treebag <- train(as.factor(class)~., data=train_set2, method="treebag", metric=metric, trControl=control)

results2 <- resamples(list( rf=fit.rf1, 
                           gbm=fit.gbm1, 
                           nb=fit.nb1, 
                           svm=fit.svmRadial1, 
                           c50=fit.c501,
                           cart = fit.cart1))
summary(results2)
# summarize the distributions
summary(results2)
# boxplots of results
bwplot(results2)
# dot plots of results
dotplot(results2)

miske_means %>% select(2:78, 83) %>% melt(id.var=c('two_class')) %>% 
  ggplot(aes(variable, value))+geom_bar(stat='identity', aes(color=two_class))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



miske_means %>% select(2:78, 83) %>% group_by(two_class) %>% summarise_each(funs(mean)) -> pairwise_summary

## identify genes with more than 5% expression in the control side. Those will be marked as +
csum<-apply(pairwise_summary,2,function(x) ifelse(as.numeric(x[1])> as.numeric(x[2])*1.05, x, NA))

## get that list with:
plusgenes<-names(which(csum != 'NA'))
plusframe <- data.frame(as.vector(plusgenes), rep('c', length(plusgenes)))
minusgenes <- setdiff(names(csum), plusgenes)
minusframe<- data.frame(as.vector(minusgenes), rep('t', length(minusgenes)))

colnames(minusframe) <- c('gene','state')
colnames(plusframe) <- c('gene','state')
finalframe <- rbind(plusframe, minusframe)
finalframe$gene <-gsub('_N','',finalframe$gene)
write.table(file='candidates.csv', finalframe,sep=',', quote=F, row.names=F)



# trans = preProcess(to_explore[,2:ncol(to_explore)-1],
#                    method=c("BoxCox", "center",
#                             "scale", "pca"))
# 
# 
# 
# PC = predict(trans, iris[,1:4])



## PCA
cc <-prcomp(na.omit(to_explore[,2:ncol(to_explore)]))
ggplot(data=as.data.frame(cc$x), aes(as.vector(cc$x[,2])))+geom_bar()
as.data.frame(cc$x) %>% melt() %>% ggplot(aes(variable,value))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("PCs")+
  ylab("PC Values")+ggtitle("PC plot")




### NMF

predictors <- miske[,2:78]

w <- matrix(1, nrow(predictors),ncol(predictors))
w [predictors[is.na(predictors)]] <- 0

## asign some random number!
predictors[is.na(predictors)] <- 99999
predictors[predictors < 0] <- 0
library(NMF)
res <- nmf(predictors, 3, 'ls-nmf', weight = w)
rx <- fitted(res)

nas <- which(is.na(as.matrix(predictors)))
t2 <- as.matrix(predictors)
t2[nas] <- rx[nas]
