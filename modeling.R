library(BayesFactor)
library(ggplot2) 
library(cowplot) # for arranging ggplot objects in grid
library(pryr) # for saving base plots and later manipulating them

source("waldFunctions.R")
source("plottingFunctions.R")

rawdata <- read.csv("data.csv")

# remove errors and potential contaminant RTs (-3 < MAD < 6) (Leys et al., 2013)
trimmed <- subset(rawdata,subset=correct==1)
medRT <- median(trimmed$rt)
mad <- mad(trimmed$rt)
data <- subset(trimmed, subset = rt>medRT-3*mad & rt<medRT+6*mad)
data$subj <- as.factor(data$subj)

#----------------------------
# plot mean RTs
#----------------------------

aggRT=aggregate(rt~subj+problemSize+format+truth,data=data,FUN="mean")
graphSummaryRT <- summarySEwithin(aggRT, measurevar="rt", withinvars=c("problemSize","format","truth"), idvar="subj")

graphSummaryRT$problemSize <- factor(graphSummaryRT$problemSize, levels = c("Small","Large"))
graphSummaryRT$truth <- factor(graphSummaryRT$truth, levels = c("True","False"))

meanRTs=ggplot(graphSummaryRT,aes(x=problemSize,y=rt,shape=format))+geom_line(aes(group=format,linetype=format))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=rt-ci,ymax=rt+ci))+facet_grid(.~truth)+labs(x="Problem size",y="Mean RT (ms)")+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))
meanRTs

plot=ggplot(graphSummaryRT,aes(x=problemSize,y=rt,fill=format))
bars=geom_bar(position=position_dodge(0.9),color="black",stat="identity")
errors=geom_errorbar(position=position_dodge(0.9),width=0.25,aes(ymin=rt-ci,ymax=rt+ci))
labels=labs(x="Problem Size",y="Mean RT (msec)")
faceting=facet_grid(.~truth)
stripFormat=theme(strip.text=element_text(face="bold",size=rel(1.5)))
legendFormat=theme(legend.title=element_text(face="bold",size=rel(1.5)),legend.text=element_text(size=rel(1.5)))
axesFormat=theme(axis.text=element_text(face="bold",size=rel(1.3)))
axesTitleFormat=theme(axis.title=element_text(face="bold",size=rel(1.4)))


basePlot=plot+bars+errors+labels+faceting+stripFormat+stripFormat+legendFormat+axesFormat+axesTitleFormat
basePlot+labs(colour="format")+theme(legend.background=element_rect(fill="white",colour="black"))+scale_fill_manual(values=c("#BBBBBB","#FFFFFF"))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))


#----------------------------
# model RTs with single mean
#----------------------------

# true problems
d1 <- subset(data,subset=truth=="True")
aggRT=aggregate(rt~subj+problemSize+format,data=d1,FUN="mean") # RT performance data aggregated by subject
RT.aov=aov(rt~problemSize*format+Error(as.factor(subj)/(problemSize*format)),data=aggRT)
summary(RT.aov)
print(model.tables(RT.aov,"means"),digits=3)

bf1 <- anovaBF(rt~problemSize*format+subj,data=d1,whichRandom="subj")
bf1
bf1[3]/bf1[4]  # favors model with no interaction by factor of 15.5

# false problems
d2 <- subset(data,subset=truth=="False")
aggRT=aggregate(rt~subj+problemSize+format,data=d2,FUN="mean") # RT performance data aggregated by subject
RT.aov=aov(rt~problemSize*format+Error(as.factor(subj)/(problemSize*format)),data=aggRT)
summary(RT.aov)
print(model.tables(RT.aov,"means"),digits=3)

bf1 <- anovaBF(rt~problemSize*format+subj,data=d2,whichRandom="subj")
bf1
bf1[3]/bf1[4]  # favors model with no interaction by factor of 13.4



#----------------------------
# model RTs with shifted Wald
#----------------------------

# True trials
# fit shifted Wald model
fit1 <- swapply(dat=d1, obsvar="rt", facs = c("problemSize", "format", "subj")) # fit the data

# model fit plots

## QQ plot (panel A)
mxqnts <- apply(fit1$xqnts,2,mean)
mpxqnts <- apply(fit1$pxqnts,2,mean)
residraw <- abs(fit1$xqnts-fit1$pxqnts); stdev <- apply(residraw,2,sd)

QQ.pryr %<a-% {
op <- par(cex.main = 1.5, mar = c(5, 6, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(as.matrix(fit1$xqnts),as.matrix(fit1$pxqnts), col="light grey",cex = 1.3,
     xlim = c(0,5000), ylim = c(0,5000), ylab = "", xlab = "", axes = FALSE)
points(mxqnts,mpxqnts,pch=20,cex=1.3)
arrows(mxqnts,mpxqnts-stdev,mxqnts,mpxqnts+stdev,angle=90,length=.025,code=3)
abline(a=0,b=1)
axis(1); mtext("RT (observed)",side=1,line=2.5,cex=1.3,font=2)
axis(2)
par(las=0)
mtext("RT (predicted)",side=2,line=4,cex=1.3,font=2)
}

## decile plot (panel B)
decile.pryr %<a-% {
ndat <- dim(fit1$vars)[1]
xse <- pxse <- sapply(1:ndat, function(ind)  sqrt(fit1$vars$alpha[ind]/(fit1$vars$gamma[ind]^3)))
resid <- abs((fit1$xqnts/xse)-(fit1$pxqnts/pxse))
dens <- lapply(1:9,function(i) density(resid[,i]))
xmax <- max(sapply(1:9,function(i) quantile(resid[,i],p=.5)  ))
ymax <- max(sapply(1:9,function(i) max(density(resid[,i])$y)  ))
op <- par(cex.main = 1.5, mar = c(5, 3, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(NA,ylim=c(0,ymax),xlim=c(0,xmax),las=1,main="",xlab="",ylab="")
invisible(sapply(1:9,function(i){
  points(dens[[i]],type="l",lwd=1.5)
  text(dens[[i]]$x[which(dens[[i]]$y==max(dens[[i]]$y))], .9*max(dens[[i]]$y),labels=paste("",i,sep=""),cex=1.3) 
}))
mtext("Std. Residual (per Decile)",side=1,line=2.5,cex=1.3,font=2)
par(las=0)
mtext("Density",side=2,line=2.5,cex=1.3,font=2)
}

# cell plot

cell.pryr %<a-% {
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
cresid <- apply(abs(resid),1,sum)
op <- par(cex.main = 1.5, mar = c(5, 4, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(NA, xlab="",ylab="",main="",las=1,ylim=c(0,1.25*max(cresid)),xlim=c(0,ndat),col=c("dark grey"))
mtext("Design cell",side=1,line=2.5,cex=1.3,font=2)
par(las=0)
mtext("Std. Residual",side=2,line=3,cex=1.3,font=2)
segments(x0=1:length(cresid),x1=1:length(cresid),y0=rep(0,times=length(cresid)),y1=cresid,col="dark grey", lwd=2)
abline(b=0,a=quantile(cresid,c(.05,.95))[2],col="black",lty=2,lwd=2)
abline(b=0,a=quantile(cresid,c(.05,.95))[1],col="black",lty=2,lwd=2)
text(ndat/2,quantile(cresid,c(.01,.978))[2],labels="95%",cex=1.3,font=2)
text(ndat/2,quantile(cresid,c(.001,.978))[1],labels="5%",cex=1.3,font=2)
points(x=ndat/2,y=mean(cresid),col="black",bg="black",pch=19,cex=1.5)
arrows(ndat/2,mean(cresid)-stderr(cresid),ndat/2,mean(cresid)+stderr(cresid),angle=90,length=.3,code=3,lwd=2)
legend(x=65,y=3.2,legend=sapply(c(bquote(bar(Delta) == .(round(mean(cresid),2))),
                                  bquote(bar(sigma)[X] == .(round(mean(xse),0))), 
                                  bquote(rho[Delta][sigma] == .(round(cor(cresid,xse),2)))
),as.expression) ,bty="n" ,xjust=0.5,cex=1.3,text.font=2)
}






# combine diagnostic plots into one plot
split.screen(c(1,3))
screen(1)
QQ.pryr

screen(2)
decile.pryr

screen(3)
cell.pryr

close.screen(all=TRUE)


# analyze drift rate (gamma)
library(dplyr)
gammas <- fit1$vars %>%
  select(gamma,problemSize,format,subj)
detach(package:dplyr)

gamma.aov=aov(gamma~problemSize*format+Error(as.factor(subj)/(problemSize*format)),data=gammas)
summary(gamma.aov)
print(model.tables(gamma.aov,"means"),digits=3)

bf <- anovaBF(gamma~problemSize*format+subj,data=gammas,whichRandom="subj")
bf
bf[4]/bf[3] # 

# analyze response threshold (alpha)
library(dplyr)
alphas <- fit1$vars %>%
  select(alpha,problemSize,format,subj)
detach(package:dplyr)

alpha.aov=aov(alpha~problemSize*format+Error(as.factor(subj)/(problemSize*format)),data=alphas)
summary(alpha.aov)
print(model.tables(alpha.aov,"means"),digits=3)

bf <- anovaBF(alpha~problemSize*format+subj,data=alphas,whichRandom="subj")
bf
bf[2]/bf[3] # 


# analyze nondecision times (theta)
library(dplyr)
thetas <- fit1$vars %>%
  select(theta,problemSize,format,subj)
detach(package:dplyr)

theta.aov=aov(theta~problemSize*format+Error(as.factor(subj)/(problemSize*format)),data=thetas)
summary(theta.aov)
print(model.tables(theta.aov,"means"),digits=3)

bf <- anovaBF(theta~problemSize*format+subj,data=thetas,whichRandom="subj")
bf
bf[2]/bf[3] # 


#############################################
## plot of Wald parameters for True problems
#############################################

# drift rates
graphSummaryGamma <- summarySEwithin(gammas, measurevar="gamma", withinvars=c("problemSize","format"), idvar="subj")
graphSummaryGamma$problemSize <- factor(graphSummaryGamma$problemSize, levels = c("Small","Large"))

driftTrue=ggplot(graphSummaryGamma,aes(x=problemSize,y=gamma,shape=format))+geom_line(aes(group=format,linetype=format))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=gamma-ci,ymax=gamma+ci))+labs(x="Problem size",y=expression(paste("Drift rate ",gamma)))+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))+theme(legend.position=c(0.75,0.8))
driftTrue

# response thresholds
graphSummaryAlpha <- summarySEwithin(alphas, measurevar="alpha", withinvars=c("problemSize","format"), idvar="subj")
graphSummaryAlpha$problemSize <- factor(graphSummaryAlpha$problemSize, levels = c("Small","Large"))

alphaTrue=ggplot(graphSummaryAlpha,aes(x=problemSize,y=alpha,shape=format))+geom_line(aes(group=format,linetype=format))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=alpha-ci,ymax=alpha+ci))+labs(x="Problem size",y=expression(paste("Response threshold ",alpha)))+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))+theme(legend.position="none")
alphaTrue

# nondecision times
graphSummaryTheta <- summarySEwithin(thetas, measurevar="theta", withinvars=c("problemSize","format"), idvar="subj")
graphSummaryTheta$problemSize <- factor(graphSummaryTheta$problemSize, levels = c("Small","Large"))

thetaTrue=ggplot(graphSummaryTheta,aes(x=problemSize,y=theta,shape=format))+geom_line(aes(group=format,linetype=format))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=theta-ci,ymax=theta+ci))+labs(x="Problem size",y=expression(paste("Nondecision time ",theta)))+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))+theme(legend.position="none")
thetaTrue

plot_grid(driftTrue,alphaTrue,thetaTrue,nrow=1,ncol=3,labels="AUTO")





#----------------------------
# shifted Wald for FALSE trials
#----------------------------

# fit shifted Wald model
fit2 <- swapply(dat=d2, obsvar="rt", facs = c("problemSize", "format", "subj")) # fit the data

# model fit plots

## QQ plot (panel A)
mxqnts <- apply(fit2$xqnts,2,mean)
mpxqnts <- apply(fit2$pxqnts,2,mean)
residraw <- abs(fit2$xqnts-fit1$pxqnts); stdev <- apply(residraw,2,sd)

QQ.pryr %<a-% {
  op <- par(cex.main = 1.5, mar = c(5, 6, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
  plot(as.matrix(fit2$xqnts),as.matrix(fit2$pxqnts), col="light grey",cex = 1.3,
       xlim = c(0,5000), ylim = c(0,5000), ylab = "", xlab = "", axes = FALSE)
  points(mxqnts,mpxqnts,pch=20,cex=1.3)
  arrows(mxqnts,mpxqnts-stdev,mxqnts,mpxqnts+stdev,angle=90,length=.025,code=3)
  abline(a=0,b=1)
  axis(1); mtext("RT (observed)",side=1,line=2.5,cex=1.3,font=2)
  axis(2)
  par(las=0)
  mtext("RT (predicted)",side=2,line=4,cex=1.3,font=2)
}

## decile plot (panel B)
decile.pryr %<a-% {
  ndat <- dim(fit2$vars)[1]
  xse <- pxse <- sapply(1:ndat, function(ind)  sqrt(fit2$vars$alpha[ind]/(fit2$vars$gamma[ind]^3)))
  resid <- abs((fit2$xqnts/xse)-(fit2$pxqnts/pxse))
  dens <- lapply(1:9,function(i) density(resid[,i]))
  xmax <- max(sapply(1:9,function(i) quantile(resid[,i],p=.5)  ))
  ymax <- max(sapply(1:9,function(i) max(density(resid[,i])$y)  ))
  op <- par(cex.main = 1.5, mar = c(5, 3, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
  plot(NA,ylim=c(0,ymax),xlim=c(0,xmax),las=1,main="",xlab="",ylab="")
  invisible(sapply(1:9,function(i){
    points(dens[[i]],type="l",lwd=1.5)
    text(dens[[i]]$x[which(dens[[i]]$y==max(dens[[i]]$y))], .9*max(dens[[i]]$y),labels=paste("",i,sep=""),cex=1.3) 
  }))
  mtext("Std. Residual (per Decile)",side=1,line=2.5,cex=1.3,font=2)
  par(las=0)
  mtext("Density",side=2,line=2.5,cex=1.3,font=2)
}

# cell plot
cell.pryr %<a-% {
  stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
  resid <- abs((fit2$xqnts/xse)-(fit2$pxqnts/pxse))
  cresid <- apply(abs(resid),1,sum)
  op <- par(cex.main = 1.5, mar = c(5, 4, 1, 1) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
  plot(NA, xlab="",ylab="",main="",las=1,ylim=c(0,3.2),xlim=c(0,ndat),col=c("dark grey"))
  mtext("Design cell",side=1,line=2.5,cex=1.3,font=2)
  par(las=0)
  mtext("Std. Residual",side=2,line=3,cex=1.3,font=2)
  segments(x0=1:length(cresid),x1=1:length(cresid),y0=rep(0,times=length(cresid)),y1=cresid,col="dark grey", lwd=2)
  abline(b=0,a=quantile(cresid,c(.05,.95))[2],col="black",lty=2,lwd=2)
  abline(b=0,a=quantile(cresid,c(.05,.95))[1],col="black",lty=2,lwd=2)
  text(ndat/2,quantile(cresid,c(.01,.978))[2],labels="95%",cex=1.3,font=2)
  text(ndat/2,quantile(cresid,c(.001,.978))[1],labels="5%",cex=1.3,font=2)
  points(x=ndat/2,y=mean(cresid),col="black",bg="black",pch=19,cex=1.5)
  arrows(ndat/2,mean(cresid)-stderr(cresid),ndat/2,mean(cresid)+stderr(cresid),angle=90,length=.3,code=3,lwd=2)
  legend(x=40,y=3.2,legend=sapply(c(bquote(bar(Delta) == .(round(mean(cresid),2))),
                                    bquote(bar(sigma)[X] == .(round(mean(xse),0))), 
                                    bquote(rho[Delta][sigma] == .(round(cor(cresid,xse),2)))
  ),as.expression) ,bty="n" ,xjust=0.5,cex=1.3,text.font=2)
}






# combine diagnostic plots into one plot
split.screen(c(1,3))
screen(1)
QQ.pryr

screen(2)
decile.pryr

screen(3)
cell.pryr

close.screen(all=TRUE)


# analyze drift rate (gamma)
library(dplyr)
gammas <- fit2$vars %>%
  select(gamma,problemSize,format,subj)
detach(package:dplyr)

gamma.aov=aov(gamma~problemSize*format+Error(as.factor(subj)/(problemSize*format)),data=gammas)
summary(gamma.aov)
print(model.tables(gamma.aov,"means"),digits=3)

bf <- anovaBF(gamma~problemSize*format+subj,data=gammas,whichRandom="subj")
bf
bf[4]/bf[3] # 

# analyze response threshold (alpha)
library(dplyr)
alphas <- fit2$vars %>%
  select(alpha,problemSize,format,subj)
detach(package:dplyr)

alpha.aov=aov(alpha~problemSize*format+Error(as.factor(subj)/(problemSize*format)),data=alphas)
summary(alpha.aov)
print(model.tables(alpha.aov,"means"),digits=3)

bf <- anovaBF(alpha~problemSize*format+subj,data=alphas,whichRandom="subj")
bf
bf[3]/bf[4] # 


# analyze nondecision times (theta)
library(dplyr)
thetas <- fit2$vars %>%
  select(theta,problemSize,format,subj)
detach(package:dplyr)

theta.aov=aov(theta~problemSize*format+Error(as.factor(subj)/(problemSize*format)),data=thetas)
summary(theta.aov)
print(model.tables(theta.aov,"means"),digits=3)

bf <- anovaBF(theta~problemSize*format+subj,data=thetas,whichRandom="subj")
bf
bf[2]/bf[3] # 


#############################################
## plot of Wald parameters for False problems
#############################################

# drift rates
graphSummaryGamma <- summarySEwithin(gammas, measurevar="gamma", withinvars=c("problemSize","format"), idvar="subj")
graphSummaryGamma$problemSize <- factor(graphSummaryGamma$problemSize, levels = c("Small","Large"))

driftFalse=ggplot(graphSummaryGamma,aes(x=problemSize,y=gamma,shape=format))+geom_line(aes(group=format,linetype=format))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=gamma-ci,ymax=gamma+ci))+labs(x="Problem size",y=expression(paste("Drift rate ",gamma)))+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))+theme(legend.position=c(0.75,0.8))


# response thresholds
graphSummaryAlpha <- summarySEwithin(alphas, measurevar="alpha", withinvars=c("problemSize","format"), idvar="subj")
graphSummaryAlpha$problemSize <- factor(graphSummaryAlpha$problemSize, levels = c("Small","Large"))

alphaFalse=ggplot(graphSummaryAlpha,aes(x=problemSize,y=alpha,shape=format))+geom_line(aes(group=format,linetype=format))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=alpha-ci,ymax=alpha+ci))+labs(x="Problem size",y=expression(paste("Response threshold ",alpha)))+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))+theme(legend.position="none")

# nondecision times
graphSummaryTheta <- summarySEwithin(thetas, measurevar="theta", withinvars=c("problemSize","format"), idvar="subj")
graphSummaryTheta$problemSize <- factor(graphSummaryTheta$problemSize, levels = c("Small","Large"))

thetaFalse=ggplot(graphSummaryTheta,aes(x=problemSize,y=theta,shape=format))+geom_line(aes(group=format,linetype=format))+geom_point(size=4)+geom_errorbar(width=0.1,aes(ymin=theta-ci,ymax=theta+ci))+labs(x="Problem size",y=expression(paste("Nondecision time ",theta)))+theme(legend.title=element_text(face="bold",size=rel(1.3)),legend.text=element_text(size=rel(1.3)))+theme(axis.title=element_text(face="bold",size=rel(1.3)))+theme(axis.text.x=element_text(size=rel(1.3)))+theme(axis.text.y=element_text(size=rel(1.3)))+theme_classic(20)+theme(axis.line.x=element_line(color="black",size=0.5,linetype="solid"),axis.line.y=element_line(color="black",size=0.5,linetype="solid"))+theme(legend.position="none")

plot_grid(driftFalse,alphaFalse,thetaFalse,nrow=1,ncol=3,labels="AUTO")




