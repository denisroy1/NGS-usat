#Update_scores.R 
#reads scored data from autoscr_locus.R (filename.rda) and reproduces 
#scanning plots for individuals. Allows individuals scores to be 
#adjusted by using the mouse to click on likely peaks.
#Data in the saved score file is updated. 

#written by dr March 2016

rm(list=ls())

#library(doParallel)
library(plyr)
library(dplyr)
library(data.table)

starttime<-Sys.time()

###################  attfr   ###################################################
attfr<-function (x) {
  wti<-paste("fr",x,sep="")
# adds "fr" to x (likely the individual id)
}
################### end attfr  #################################################

################### setyincr  #################################################
setyincr<-function(x) {
  yincr<-50
  # default yincr value
  
  if (x>=1000) yincr<-100
  if (x<=150 && x) yincr<-25
  if (x<=100 && x) yincr<-20
  if (x<=50 && x) yincr<-10
  if (x<=20 && x) yincr<-5
  if (x<=10) yincr<-2
  
  return(yincr)
  # sets the y-increment for plotting
}
################## end setyincr ##############################################

rdafile<-file.choose()
fd<-dirname(rdafile)
setwd(fd)
load(rdafile)
# asks user for the .rda file saved from previous automatic scoring

scores_file<-paste(fd,"/",mrkr,".scores",sep="")
prvscr<-read.table(scores_file,header=T)
# accesses the '.scores' file for updating 

redolist<-file.choose()
#redos<-scan(redolist)
redos<-scan(redolist, character(), blank.lines.skip = T)
rnum<-length(redos)

# locates, scans, and counts the redo list 
# (list of individuals with problematic scores)

fnum<-length(cfd)

frredos<-attfr(redos)
# reads in list of individuals to rescore

pti<-paste(mrkr,"-",redos,sep="")
# builds individual redo names 

cfdindx<-unlist(dimnames(cfd)[2])
# makes cfd indices searchable

subdir<-paste(mrkr,'-traces',sep='')
setwd(file.path(fd,subdir))
# locate and use subdir to re-store traces

rin<-0
rscr<-0
# initialises rin to 0

for (a in 1:rnum) rin[a]<-which(cfdindx==frredos[a])
for (b in 1:rnum) rscr[b]<-which(prvscr$indiv==redos[b])
#minor loop to find subset of individials to re-plot

subcfd<-cfd[rin]
cmy<-pmxys[rin]
als<-cbind(as.numeric(prvscr$al1[rscr]),as.numeric(prvscr$al2[rscr]))
alfr<-cbind(as.numeric(prvscr$as1[rscr]),as.numeric(prvscr$as2[rscr]))
sista<-as.numeric(prvscr$si[rscr])
mph<-minpeakheight[rscr]
# setting appropriate variables from loaded data 

for (b in 1:rnum) {

  trname<-paste(pti[b],"-trace",".jpeg",sep="")
  # sets name of the file to be saved in trace dir
  
  par(oma=c(0,0,0,0))    
  par(mar=c(3.7,3.7,0.5,0.5))  
  
  # sets plot window as in automatic scoring
  
  jtfy<-.5
  # default justification for title
  jtfy<-.5
  #justify title placement
  
  if (length(xntv)<=20) { 
    locgenran<-which(als[b,]==tail(xntv, 6))
  } else {
    locgenran<-which(als[b,]==tail(xntv, 16))
  }
  
  if (length(locgenran)>=1) jtfy<-.1 else jtfy<-.9
  #tests if alleles called within 16 repeat of end, if yes makes adjustment to jtfy
  # adjusts justification based on where alleles are
  
  yincr<-setyincr(cmy[b])
  # sets y increment
  
  numb<-max(alfr[b,])
  # finds max frequency in individual redo b
  
  cals<-als[b,]
  # sets which individual's alleles to be changed
  
  plot(cfd$alleles,cfd$mxfr, type='l', lwd=0.5, xlab="", ylab="", col='green',xaxt='n', yaxt='n', 
       frame.plot=F, ylim=c(0,(cmy[b]+10)))
  axis(side = 1, at = xntv, lwd=2, cex.axis=1.2, cex.lab=2)
  axis(side = 2, at = seq(0,cmy[b]+10,by=yincr), lwd=2,cex.axis=1.2)
  title(xlab='# of repeats', line=2.3, cex.lab=1.4)
  title(ylab='frequency', line=2.5, cex.lab=1.4)
  for (p in 2:fnum) lines(cfd$alleles,cfd[,p], lwd=0.5, lty=3, col='gray')
  lines(cfd$alleles,subcfd[,b], type ="b", col='dodgerblue4', lwd=3, pch=20, cex=1)
  title(outer=F,adj=jtfy, main=(pti[b]), cex.main=1.4, font=12, col.main='darkgreen',line=-1)
  text(cals,(numb*1.07),as.character(cals), cex=1.1, col="lightsteelblue")
  arrows(x0=0, y0=mph[b], x1=tail(xntv,1), y1=mph[b], col='orange',lty=2, code=0,lwd=2)
  # plots the frequency distribution of individual relative to that of others 
  # including the minpeakheight
  
  posal<-locator(2, "n")
  cals<-round(posal$x)
  numb<-round(posal$y)
  # get's the new values from the users interaction with the graph
  
  arrows(x0=cals, y0=max(numb),x1=cals,y1=-30, col='red', code=0,lwd=3)
  text(cals,(max(numb)*1.07),as.character(cals), cex=1.4, col="darkred")
  # plots the new scores on the graph 
  
  if(cals[1]==cals[2]) hetsta<-0 else hetsta<-1
  sista[b]<-sista[b]+1
  # updates whether individual has changed to homo or heterozygote
  # and how many times it's been re-scored
  
  prvscr[rin[b]-1,]$al1<-as.character(cals[1])
  prvscr[rin[b]-1,]$al2<-as.character(cals[2])
  
  prvscr[rin[b]-1,]$as1<-as.character(numb[1])
  prvscr[rin[b]-1,]$as2<-as.character(numb[2])
  
  prvscr[rin[b]-1,]$xev1<-as.character(round(numb[1]/mph[b],2))
  prvscr[rin[b]-1,]$xev2<-as.character(round(numb[2]/mph[b],2))
  
  prvscr[rin[b]-1,]$het<-as.character(hetsta)
  prvscr[rin[b]-1,]$si<-as.character(sista[b])
  # above 8 commands updates the scores
  
  dev.copy(jpeg,trname)
  dev.off()
  # save the newly produced jpeg of the individual's scores
}

setwd(fd)
write.table(prvscr, file=paste(mrkr,".scores",sep=""), append=F, quote=F, sep="\t", eol="\n", col.names=T, row.names=F)
#updating scores to the '.scores' file

Sys.time() - starttime
#time stamp
