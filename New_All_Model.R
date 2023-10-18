rm(list = ls())
library(ggplot2)
library(ggsci)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(devtools)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(snowfall) #+
load('Data/New-Model-prepare.rda')
mm <- lapply(ee,function(x){
  rownames(x) <- x[,1]    ###给与行名
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})
##################################
#### 准备工作 ####
##################################

result <- data.frame()
est_data <- mm$TCGA_LUAD   ###修改建模集
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',pre_var)]})
rm(mm)

rf_nodesize <- 3
seed <- 10  ##需要改Seed数

##################################
#### 1-1.RSF ####
##################################

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result,cc)


##################################
#### 1-2.rsf+enet ####
##################################

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
vi <- data.frame(imp=vimp.rfsrc(fit)$importance)
vi$imp <- (vi$imp-min(vi$imp))/(max(vi$imp)-min(vi$imp))
vi$ID <- rownames(vi)

ggplot(vi,aes(imp,reorder(ID,imp)))+
  geom_bar(stat = 'identity',fill='#FF9933',color='black',width=0.7)+
  geom_vline(xintercept = 0.4,color='grey50',linetype=2)+
  labs(x='Relative importance by Random Forest',y=NULL)+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.x = element_text(size = 11,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.title = element_text(size=13,color='black'),
        legend.text = element_text(size=12,color='black'),
        legend.title = element_text(size=13,color='black'))+
  scale_y_discrete(expand = c(0.03,0.03))+
  scale_x_continuous(expand = c(0.01,0.01))
ggsave(filename = 'rf-importance.pdf',width = 5.5,height = 8)

rid <- rownames(vi)[vi$imp>0.2]
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}


##################################
#### 1-3.rsf+stepcox ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

##################################
#### 1-4.rsf+CoxBoost ####
##################################

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
plot(fit)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result,cc)

##################################
#### 1-5.rsf+plsRcox ####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = F)

fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=3)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + plsRcox')
result <- rbind(result,cc)

##################################
#### 1-6.rsf+superpc ####
##################################

data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=3,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + SuperPC')
result <- rbind(result,cc)


##################################
#### 1-7.rsf+gbm ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
plot(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + GBM')
result <- rbind(result,cc)

##################################
#### 1-8.rsf+survivalsvm ####
##################################

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + survival-SVM')
result <- rbind(result,cc)

##################################
#### 2-1.Enet ####
##################################

x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}

##################################
#### 2-2.Lasso+RSF####
##################################

set.seed(seed)
fit = cv.glmnet(x1, x2,family = "cox",alpha=1,nfolds = 10)
coef.min = coef(fit, s = "lambda.min") 
active.min = which(as.numeric(coef.min)!=0)
rid <- colnames(x1)[active.min]

est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'Lasso + RSF'
result <- rbind(result,cc)

##################################
#### 2-3.Lasso+StepCox ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

##################################
#### 2-4.Lasso+CoxBoost ####
##################################

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result,cc)

##################################
#### 2-5.Lasso+plsRcox ####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = F)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + plsRcox')
result <- rbind(result,cc)

##################################
#### 2-6.Lasso+superpc ####
##################################

data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=3,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + SuperPC')
result <- rbind(result,cc)


##################################
#### 2-7.Lasso+gbm ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + GBM')
result <- rbind(result,cc)

##################################
#### 2-8.Lasso+survivalsvm ####
##################################

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + survival-SVM')
result <- rbind(result,cc)


##################################
#### 3-1.StepCox ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}

##################################
#### 3-2.StepCox+RSF ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
               ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + RSF')
  result <- rbind(result,cc)
}

##################################
#### 3-3.StepCox+Enet ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  x1 <- as.matrix(est_dd2[,rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))
  
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
    rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
    
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox','[',direction,']',' + Enet','[α=',alpha,']')
    result <- rbind(result,cc)
  }
}

##################################
#### 3-4.StepCox+CoxBoost ####
##################################

for (direction in c("backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit <- CoxBoost(est_dd2[,'OS.time'],est_dd2[,'OS'],as.matrix(est_dd2[,-c(1,2)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox',' + CoxBoost')
  result <- rbind(result,cc)
}

##################################
#### 3-5.StepCox+plsRcox ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  
  set.seed(seed)
  cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = F)
  fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + plsRcox')
  result <- rbind(result,cc)
}

##################################
#### 3-6.StepCox+superpc ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  
  data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                       n.fold = 10,
                       n.components=3,
                       min.features=3,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  rs <- lapply(val_dd_list2,function(w){
    test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2],RS=rr)
    return(rr2)
  })
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + SuperPC')
  result <- rbind(result,cc)
}

##################################
#### 3-7.StepCox+gbm ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + GBM')
  result <- rbind(result,cc)
}

##################################
#### 3-8.StepCox+gbm ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
  rid <- names(coef(fit))
  est_dd2 <- est_data[,c('OS.time','OS',rid)]
  val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})
  
  fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + survival-SVM')
  result <- rbind(result,cc)
}


##################################
#### 4-1.CoxBoost ####
##################################

set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)


##################################
#### 4-2.CoxBoost+Enet ####
##################################

rid <- names(coef(fit)[which(coef(fit)!=0)])
est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

x1 <- as.matrix(est_dd2[,rid])
x2 <- as.matrix(Surv(est_dd2$OS.time,est_dd2$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}


##################################
#### 4-3.CoxBoost+stepcox ####
##################################

for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,est_dd2),direction = direction)
  rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

##################################
#### 4-4.CoxBoost+RSF ####
##################################

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'CoxBoost + RSF'
result <- rbind(result,cc)

##################################
#### 4-5.rsf+plsRcox ####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd2[,rid],time=est_dd2$OS.time,status=est_dd2$OS),nt=10,verbose = F)
fit <- plsRcox(est_dd2[,rid],time=est_dd2$OS.time,event=est_dd2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + plsRcox')
result <- rbind(result,cc)

##################################
#### 4-6.CoxBoost+superpc ####
##################################

data <- list(x=t(est_dd2[,-c(1,2)]),y=est_dd2$OS.time,censoring.status=est_dd2$OS,featurenames=colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=3,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + SuperPC')
result <- rbind(result,cc)

##################################
#### 4-7.CoxBoost+gbm ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + GBM')
result <- rbind(result,cc)

##################################
#### 4-8.CoxBoost+survivalsvm ####
##################################

fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + survival-SVM')
result <- rbind(result,cc)

##################################
#### 5.plsRcox####
##################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=est_dd[,pre_var],time=est_dd$OS.time,status=est_dd$OS),nt=10,verbose = F)
fit <- plsRcox(est_dd[,pre_var],time=est_dd$OS.time,event=est_dd$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)

##################################
#### 6.superpc####
##################################

data <- list(x=t(est_dd[,-c(1,2)]),y=est_dd$OS.time,censoring.status=est_dd$OS,featurenames=colnames(est_dd)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=3,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)


##################################
#### 7.GBM ####
##################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = est_dd,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)


##################################
#### 8.survivalsvm ####
##################################
fit = survivalsvm(Surv(OS.time,OS)~., data= est_dd, gamma.mu = 1)

rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')
result <- rbind(result,cc)
result2 <- result
####over####
save(result2,file = '7.New_All_Model/New-all-model.rda') 

load('7.New_All_Model/New-all-model.rda')  
result2$Model <- gsub('α','a',result2$Model)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
range(result2$Cindex)
result2%>%filter(ID!='TCGA_LUAD')%>%    ###建模集
  ggplot(aes(Cindex,reorder(Model,Cindex)))+
  geom_bar(width = 0.7,stat = 'summary',fun='mean',fill='orange2')+
  theme_classic()+
  labs(y=NULL)

dd <- result2%>%
  filter(ID!='TCGA_LUAD')%>%   ###建模集
  group_by(Model)%>%
  summarise(Cindex=mean(Cindex))

dd%>%
  ggplot(aes(Cindex,reorder(Model,Cindex)))+
  geom_bar(width=0.7,stat = 'identity',fill='orange')
  #+scale_x_break(c(0.05,0.945),scales = 20)    ###scale_x_break（）修改


dd2 <- pivot_wider(result2,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()
dd2[,-1] <- apply(dd2[,-1],2,as.numeric)

table(result2$Model)

# -------------------------------------------------------------------------
###把选择的最佳模型复制粘贴此处。   Enet [a=0]
#### Enet [a=0] ####
##################################
x1 <- as.matrix(est_dd[,pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time,est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,family = "cox",alpha=0,nfolds = 10)
coef.apprx = coef(fit, s = fit$lambda.min)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
coef.apprx@x
save(fit,rs,file = 'Data/New-final-Enet[a=0].rda')

