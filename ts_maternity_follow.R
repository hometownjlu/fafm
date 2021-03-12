#time series of maternity follow

library(compareGroups)
library(Nbclust)
library(HDclassif)
library(sparcl)
library(forecast)

library(arm)
setwd("D:/R test")


tcga <-read.csv("ts-mas-wd-sd.csv")
tcga=tcga[,-1]

full=glm(num~wd+sd+bus+stu+jj+fes,family =poisson,data = tcga)
full2=glm(num~stu,family = poisson,data = tcga)
full
summary(full)
summary(full2)
plot(full$residuals)
plot(full$residuals)
plot(full$residuals)
plot(full$residuals)
plot(full$residuals)
plot(full)

mean(tcga$num)
var(tcga$num)
hist(tcga$num)
head(tcga)
par(tcga=c(2,2))
par(mfrow=c(1,1))
plot(tcga)
?par

attach(tcga)
aovr=aov(num~fes)
summary(aovr)
library(gplots)
plotmeans(num~wd)
install.packages("stats")
install.packages("jtools")
install.packages("broom")
install.packages("ggstance")
install.packages("interactions")

library(broom)
library(ggstance)
library(jtools)
library(interactions)
full=glm(num~wd+sd+bus+stu+jj+fes,family =poisson,data = tcga)
full2=glm(num~stu,family = poisson,data = tcga)
full
summary(full2)

plot_summs(full,scale=TRUE,exp=TRUE)

interact_plot(full,pred=wd,modx=fes)
cat_plot(full,pred=wd,modx=sd)

?cat_plot



coef1=coef(full)
se.coef1=se(coef1)

coef2=coef(full2)
se.coef2=se.coef(full2)

cmodel=cbind(coef1,coef2,exponent=exp(coef1))
cmodel

newdata=data.frame(wd=-18,sd=15,bus=0,stu=1,jj=2,fes=0)
predict(full,newdata = newdata,type="response")



r=summary(full)$r.sq
r
par(mfrow=c(2,2))
plot(full)

summary(full)
如何处理离群点
library(car)
install.packages(car)
car::influencePlot(full)
full
coef(full)
exp(coef(full))
?glm

cor(tcga)
pairs(tcga)



library(psych)
corr.test(tcga)
Mat_fol=read.csv('tst-mas-w.csv')
Mat_fol=Mat_fol[,-1]
Mat_folts=ts(Mat_fol,frequency =53)
Mat_folts=ts(Mat_fol)
plot.ts(Mat_folts)
Mat_folts

Mat_fol=read.csv('tst-mas-w1.csv')
Mat_fol=Mat_fol[,-1]
Mat_folts=ts(Mat_fol,frequency =4)
Mat_folts=ts(Mat_fol)
plot.ts(Mat_folts)
Mat_folts
Forecast=forecast(Mat_fol,h=12)
plot(Forecast)

Mat_fol=read.csv('ts-mas12.csv')
Mat_fol=read.csv('ts-mas12-bak.csv')
Mat_fol=Mat_fol[,-1]

summary(Mat_fol)
names(Mat_fol)
Mat_fol

Mat_folts=ts(Mat_fol,frequency =12)
Mat_folts=ts(Mat_fol,frequency =53)
Mat_folts=ts(Mat_fol)
plot.ts(Mat_folts)
Mat_folts
#Mat_folts=ts(Mat_fol)
#plot(Mat_fol)

mean(Mat_fol)
var(Mat_fol)

#decompose random seasonal trend
Dec_Fol=decompose(Mat_folts)
Dec_Fol$seasonal
plot(Dec_Fol)

#del seasonal and random factor
Dec_Foldel=Mat_fol-Dec_Fol$seasonal-Dec_Fol$random
plot(Dec_Foldel)
Dec_Fol$seasonal
Dec_Fol
plot(Dec_Fol$seasonal)
plot(Dec_Fol$random)
plot(Dec_Fol$x)

plot(Dec_Fol$trend)
plot(Dec_Fol$figure)
Dec_Fol$figure
plot(Dec_Fol$type)
Dec_Fol$type
Dec_Foldel=Mat_fol-Dec_Fol$seasonal-Dec_Fol$random
plot(Dec_Foldel)

#forecast

Forecast=forecast.HoltWinters(Mat_folts,h=12)
plot(Forecast)

tsdiag(Mat_fol)
##去平滑时间序???
library(TTR)
Mat_folma=SMA(Mat_folts,n=5)
plot.ts(Mat_folma)
#log for smoonth
logbirthsts=log(Mat_folts)
logbirthsts
plot.ts(logbirthsts)
Forecast=forecast(logbirthsts,h=12)
Forecast=forecast(Mat_folma,h=12)
plot(Forecast)

#“l.start”和“b.start”的参数指定水平和趋势的初始值
#常见的设定水平初始值为时间序列的第一个值（608）
#斜率的初始值则是其第二个值减去第一个值（9）
Mat_fol_forest=HoltWinters(Mat_folts,gamma = FALSE,l.start = 580,b.start = 10)
Mat_fol_forest=HoltWinters(Mat_folts,gamma = FALSE)
Mat_fol_forest
Mat_fol_forest$SSE
#sse 是误差
Mat_fol_forest$fitted
plot(Mat_fol_forest)
#
fcsthot=forecast(Mat_fol_forest,h=12)
Forecast=forecast:::forecast.HoltWinters(Mat_fol_forest,h=12)
plot(Forecast)
Forecast$residuals
Forecast

#log for smoonth
logMat_fol=log(Mat_folts)
Mat_fol_forest=HoltWinters(logMat_fol,gamma = FALSE,l.start = 580,b.start = 10)
Mat_fol_forest
Mat_fol_forest$SSE
#sse 是误差
Mat_fol_forest$fitted
plot(Mat_fol_forest)



#两个大、小阴影区间分别代表 80,95%的可能性
#lag.max 是预测误差延迟1-20阶的相关图
acf(Forecast$residuals,na.action = na.pass,lag.max = 20)
acf(Forecast,na.action = na.pass,lag.max = 20)
Forecast
plotForecastErrors(Forecast$residuals)
Box.test(Forecast$residuals,lag=20,type = "Ljung-Box")
#P值小???0.05 ，证??? 预测样本内预测误差在1-20阶内???0自相关？？？
#p??? 0.114
pacf(Forecast$residuals,lag.max = 20)
plot.ts(Forecast$residuals)
#residuals 正态分布，且均值为0，则模型不可再优化
hist(Forecast$residuals)
#用图形的方式显示正太分布


#两个大、小阴影区间分别代表 80,95%的可能性
#lag.max 是预测误差延迟1-20阶的相关图
acf(Forecast$residuals,na.action = na.pass,lag.max = 20)
acf(Forecast,na.action = na.pass,lag.max = 20)
Forecast
plotForecastErrors(Forecast$residuals)
Box.test(Forecast$residuals,lag=20,type = "Ljung-Box")
#P值小???0.05 ，证??? 预测样本内预测误差在1-20阶内???0自相关？？？
#p??? 0.114
pacf(Forecast$residuals,lag.max = 20)
plot.ts(Forecast$residuals)
#residuals 正态分布，且均值为0，则模型不可再优化
hist(Forecast$residuals)
#用图形的方式显示正太分布


#预测模型构建 & 机器学习（R语言进阶???
#time sequence 
# holt ，arima,model to smooth,and then forecast
#box.test or hist is applied to test performance


# diff 1阶差
bir_dif=diff(Mat_folts,differences = 1)
plot.ts(bir_dif)
#差分后，如不平稳再加1???
ski_dif=diff(Mat_folts,differences = 2)
plot.ts(ski_dif)

#holt 只预测已有时间内的数据，不预测未来数据。forecast 预测将来



auto.arima(Mat_folts)
ski_ari=arima(Mat_folts,order=c(0,1,0))
plot.ts(ski_ari)
ski_ari
ski_ari_fo=forecast(ski_ari)
plot(ski_ari_fo)
Box.test(ski_ari_fo$residuals,lag=20,type = 'Ljung-Box')
hist(ski_ari_fo$residuals)

tsdiag(ski_ari)
#Ljung p 值大于0.1时就是准确的。
