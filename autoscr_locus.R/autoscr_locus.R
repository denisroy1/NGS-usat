#autoscr_locus.R
#Scanning frequency distribution created by prepare_data.R 
#and to find most likely alleles

#written by dr Feb 2016

rm(list=ls())

#library(doParallel)
library(foreach)
library(plyr)
library(dplyr)
library(data.table)
library(Biostrings)

###################### fid #####################################################
fid<-function(x) {
  if (as.numeric(x)>3) x<-as.numeric(substr(x,nchar(x)-2,nchar(x)))
  #fixes the id number of individuals to be represented byu 3 digits
  
  return(x)
}
###################### end fid #################################################
#applies to 'ionxpressXX' files

extid<-function(x) {
  dash<-as.numeric(gregexpr2("-",x))
  take<-substr(x,1,dash-1)
}

#################### corf #####################################################
corf<-function(y, mph=0) {
  nmxy<-which(y==max(y))
  #finds the max value(s) of y 
  
  if (length(nmxy) > 1) y[min(nmxy)] <- y[min(nmxy)]+10
  #if more than 1 max adds 10 to the first occurance to break ties
 
   return(y)  
}
################### end corf ##################################################

###################### sdf #####################################################
sdf<-function(x) {
  sdl<-paste(as.character(sign(diff(x))),collapse = "")
  #converts fqdt to +1 and -1s
  
  sdl <- gsub("1", "+", gsub("-1", "-", sdl))
  #converts 1s to +,- and 0s
  
  return(sdl)
}
##################### end sdf #################################################

##################### alst ####################################################
alst<-function(x, pat, y) {
  
  y<-y[,]
  #vectorizes y
  
  rc<-gregexpr(pat,x)[[1]]
  #looks for '+-' and continues along until + or 0 encountered
  #in each ccfd
  
  if (rc[1] < 0) {
    al<-c(0,0,0,0)
  }else{
  y1<-rc
  #finds where pattern begins and how long it goes
  
  y2<-rc + attr(rc,"match.length")
  #adds the length of the peak to its beginning
  
  attributes(y1)<-attributes(y2)<-NULL
  #removes attributes  
  
  n<-length(y1)
  yv <- yp <- numeric(n)
  
  for (i in 1:n){
    frg<-y[y1[i]:y2[i]]
    yp[i]<-which.max(frg)+y1[i]-1
    #finds position of max in 'frg' and -1 to account for initial shift
  
    yv[i]<-y[yp[i]]
    #returns the value of y at yp  
    
}
  al<-cbind(yp,yv,y1,y2)
  #saves the recovered peaks as possible alleles
}
  return(al)
}
#################### end alst #################################################

##################   rfp #######################################################
rfp<-function (x) {
  pmx<-0 
  if (x > 350) pmx<-round_any(x, 100, f=ceiling)
  if (x <= 350 && x > 100) pmx<-round_any(x, 50, f=ceiling)
  if (x <= 100 && x > 50) pmx<-round_any(x, 20, f=ceiling) 
  if (x <= 50 && x > 20) pmx<-round_any(x, 10, f=ceiling)
  if (x <= 20 && x > 10) pmx<-20
  if (x <= 10) pmx<-10
  
  if (pmx-x < 20 && pmx > 350) pmx<-pmx + 100
  if (pmx-x < 10 && length(!which(pmx==seq(101,350)))) pmx<-pmx + 50
  if (pmx-x < 10 && length(!which(pmx==seq(51,100)))) pmx<-pmx + 20
  if (pmx-x < 10 && length(!which(pmx==seq(21,50)))) pmx<-pmx + 10
  if (pmx-x < 5 && length(!which(pmx==seq(11,20)))) pmx<-pmx + 5
  if (pmx-x < 2 && length(!which(pmx==seq(0,10)))) pmx<-pmx + 2
  
  return (pmx)
  #rounds and pads the max frequency values for plotting so numbers
  #generally stay on the graph
}
###################  end rfp  ################################################

################### setyincr  #################################################
setyincr<-function(x) {
  yincr<-50
  #default yincr value
  
  if (x>=1000) yincr<-200
  if (x<=150 && x) yincr<-25
  if (x<=100 && x) yincr<-20
  if (x<=50 && x) yincr<-10
  if (x<=20 && x) yincr<-5
  if (x<=10) yincr<-2
  
  return(yincr)
}
################## end setyincr ##############################################

###################  alviz  ###################################################
alviz<-function(xs,x,mxs,xntv,cmy,yincr,pti,mph,AL,mrkr,fdd) {

  subdir<-paste(mrkr,'-traces',sep='')
  #sets subdir off main to store traces
  
  if (file.exists(subdir)) {
    setwd(file.path(fdd,subdir))
  } else {
    dir.create(file.path(fdd,subdir),showWarnings = F)
    setwd(file.path(fdd,subdir))
  }
  #is subdir exists then use it, otherwise create it and setwd to it
  
  trname<-paste(gsub(" ","",pti),"-trace",".jpeg",sep="")
  #sets name of the file to be saved in trace dir
  
  jpeg(trname, width=10, height=10, units="cm", pointsize=12, bg="white", quality=90, res=250)
  #set type of file with specific to save in traces dir
  
  #xs<-cfd
  #x<-cfd[k]
  #mxs<-mxfr
  #cmy<-pmxys[k]
  #mph<-minpeakheight[k]
  #above lines not read - for testing purposes
  
  x<-x[,]
  #vectorizes x
  
  xsl<-length(xs)-1
  #finds length of of dataframe passed (i.e., cfd)
  
  als<-xntv[AL[,1]]
  #gets saved alleles
  
  numb<-max(AL[,2])
  #gets max freq
  
  jtfy<-.5
  #justify title placement
  
  if (length(xntv)<=20) { 
    locgenran<-which(als==tail(xntv, 6))
    } else {
    locgenran<-which(als==tail(xntv, 16))
    }
  
  if (length(locgenran)>=1) jtfy<-.1 else jtfy<-.9
  #tests if alleles called within 16 repeat of end, if yes makes adjustment to jtfy
  # adjusts justification based on where alleles are
  
  if (cmy >= 1000) {
    cmy<-cmy+100
    yincr<-200
  } else {
    cmy<-cmy+10
  }
  #add empty space padding to max values for plotting
  
  par(oma=c(0,0,0,0))    
  par(mar=c(3.7,3.7,0.5,0.5))  
  #sets plotting areas
  
  plot(xs$alleles, mxs, type='l', lwd=0.5, xlab='', ylab='',col='green', xaxt='n', yaxt='n', 
       frame.plot=F, ylim=c(0,(cmy+10)))
  title(xlab='# of repeats', line=2.3, cex.lab=1.4)
  title(ylab='frequency', line=2.5, cex.lab=1.4)
  axis(side = 1, at = xntv, lwd=2, cex.axis=1.2, cex.lab=2)
  axis(side = 2, at = seq(0,cmy+10,by=yincr), lwd=2,cex.axis=1.2)
  title(outer=F,adj=jtfy, main=(pti), cex.main=1.4, font=12, col='black',line=-1)
  #sets defaults for plot with axes and titles using max values
  
  for (p in 2:xsl) lines(xs$alleles,xs[,p], lwd=0.5, lty=3, col='gray')
  #plots all freq. dist. for comparisons
  
  if (substr(pti,(nchar(pti)-3),nchar(pti))!="mxfr") {
    lines(xs$alleles,x, type ="b", col='dodgerblue4', lwd=3, pch=20, cex=1)
    arrows(x0=als, y0=numb[1],x1=als,y1=-30, col='red', code=0,lwd=3)
    text(als,(numb[1]*1.07),as.character(als), cex=1.1, col="darkred")
  } else {
    
    lines(xs$alleles,x, type ="l", col='green', lwd=2, cex=1) 
  }
  arrows(x0=0, y0=mph, x1=tail(xntv,1), y1=mph, col='orange',lty=2, code=0,lwd=2)
  #does actual plotting with special conditions for the mxfr plot
  #also shows minpeakheight from which x - eveneness is calculated
  
  dev.off()
  #shuts jpeg creation off
}
###################  end alviz  ################################################

#cl<-makeCluster(6)
#registerDoParallel(cores=6)

starttime<-Sys.time()

fdfile<-file.choose()
fdd<-dirname(fdfile)
setwd(fdd)
#get the fd list file generated from prepdat

fdl<-read.table(fdfile, header=F)
attach(fdl)
fnum<-nrow(fdl)
#reads in fd list file and finds the number of entries

scrfile<-paste(substr(basename(fdfile),1,nchar(basename(fdfile))-7),".scores", sep='')
#generate the name of the score file where scores will be dumped

cfdfile<-paste(substr(basename(fdfile),1,nchar(basename(fdfile))-7),".cfds", sep='')
#generate the name of the cfdfile where cfd data frame will be dumped

pdst<-apply(fdl,1,as.character)
prid<-lapply(pdst,basename)
iid<-lapply(prid,extid)
#id<-sapply(as.numeric(gsub("\\D","",prid)),fid,USE.NAMES=F)

#mrkrl<-head(regexpr('-',prid[1]))
#mrkr<-substr(prid[[1]],1,mrkrl-1)
mrkr<-substr(basename(fdfile),1,nchar(basename(fdfile))-7)

foreach (n = 1:fnum, .combine=c, .export="cfd" ) %do% { 
#cfd<- foreach(s = 1:fnum, .combine='cbind', .inorder=T) %dopar%   
#above comment is for parallel implementation which cannot work
#when output depend on each other - no easy way to combine data
  
  f<-read.table(pdst[n], header=T)
  attach (f)
  if(max(f$Freq)==head(f$Freq,1)) f<-rbind(c(head(f$alleles,1)-1,0),f)
  colnames(f)[2]<-paste('fr', as.character(iid[n]),sep='')
  if(n==1) { 
    cfd<-f 
  } else {
  cfd<-merge(cfd,f, by='alleles', all=T)
  }
}  
#loop combining all frequency dist. into 1 giant dataframe (df)

cfd[is.na(cfd)]<-0   
#removing NAs from cfd

mna<-min(cfd$alleles)
mxa<-max(cfd$alleles)
xntv<-seq(mna,mxa,by=1)
if (mxa%%2!=0) xntv<-c(xntv,mxa+1)
#finding the population allele span

if (!length(xntv)==length(cfd$alleles)) {
  misal<-xntv[!(xntv %in% cfd$alleles)]
  mtyal<-data.frame(matrix(0,ncol=ncol(cfd),nrow=length(misal)))
  mtyal[,1]<-misal
  colnames(mtyal)<-names(cfd)
  cfd<-arrange(rbind(cfd,mtyal),alleles)
}
#condition to make alleles span gaps in allele dist.  

mxfr<-as.integer(apply(cfd[,-1], 1, max))
#finding the max ferquency at all alleles 

cfd<-cbind(cfd,mxfr)
#putting mxfr at end of cfd

nups<-ndowns<-1
peakpat=NULL
npeaks<-2
dth<-0.2
AL<-matrix(NA, nrow=npeaks, ncol=4, dimnames=list(c('al1','al2'),c('yp','yv','y1','y2')))
#df to temporarily store called alleles

scores<-as.data.frame(matrix(NA,nrow=(fnum+1),ncol=10,byrow=F,
dimnames=list(c(),c("mrkr","indiv","al1","al2","as1","as2","xev1","xev2","het","si"))))
#df that will from basis for individual scores and results - exportable

si<-1
#variable outlining score iteration

het<-1
#variable outlining hetero (=1) or homozygote (=0); default is 1

#id<-c(0,id)
lsid<-as.list(c(iid,'mxfr'))
#adding to and making list out of individual ids

mxys<-apply(cfd, 2, max)
#vector of each inds max allele freq

sumys<-apply(cfd, 2, sum)
#vector of each inds sum freq

minpeakheight<-sumys %/% length(xntv)
minpeakheight<-replace(minpeakheight,which(minpeakheight<2),2)
#vector of each inds minpeakheight

ccfd<-apply(cfd,2,corf)
ccfd<-apply(ccfd,2,sdf)
#using corf and sdf functions above to genetate compressed corrected freq dists

mxindx<-length(ccfd)

if (is.null(peakpat)) peakpat <- sprintf("[+]{%s,}[-]{%s,}", nups, ndowns)  
#sets peak pattern to recognize

pmxys<-sapply(mxys,rfp)
#calls rfp to fix max values for plotting

for (k in 2:mxindx) {
  tmp<-cfd[k]
  ctmp<-ccfd[k]
  mph<-minpeakheight[k]
  mxy<-mxys[k]
  for (a in 1:npeaks){
    pals<-alst(ctmp, peakpat,tmp)
    if (is.null(nrow(pals))) ca<-c(0,0,0,0) else ca<-pals[which.max(pals[,2]),]
    #if  pals is empty then fill cal with 0s otherwise use largest of pals
    
    if (a==2){
      if (ca[1] >= AL[1,1]) dth<-dth*0.5 else dth<-dth*2.5
    #  mph<-sum(tmp[-(AL[1,3]:AL[1,4]),]) %/% length(tmp[-(AL[1,3]:AL[1,4]),])
    }
    
    #default threshold adjusted based on 2nd allele's post- or 
    #pre- position. Pre- threshold set higher, while post- threshold 
    #set lower
    
    #if (mph<2) mph<-2
    
    thd<-as.integer(mxy*dth)
    #this is where the threshold gets calculated
  
    evalthd<-ca[2]-pmax(tmp[ca[3],],tmp[ca[4],])
    if (length(evalthd)<=0) evalthd<-0
    #this is where threshold is used to evaluate alleles
     
    if (ca[2] >= mph && evalthd >= thd) {
      AL[a,]<-ca
    } else {
      AL[a,]<-0
    }
  #allows only alleles above minpeakheight and passed threshold
  #(established by the difference between the peak and largest margin)
  #to be called. 
    
    if (a==1 && nrow(pals)==1) {
      tmp<-as.data.frame(rep(0,length(cfd$alleles)))
      ctmp<-sdf(corf(tmp))
    }
    #handles unlikely condition of only one possible allele passed
    
    if (a==1 && nrow(pals) > 1) {
      pt<-tmp[ca[3]-1,]
      if(length(pt)<=0) pt<-0
      pk2<-max(pals[,2][-which.max(pals[,2])])
      if (pt >= pk2 && abs(pt-tmp[ca[3],]) <= 0.2*(ca[2]-tmp[ca[3],])) {
        x<-ca[3]
        while(tmp[x,]>2) {x<-x-1; if(length(tmp[x,])==0 || tmp[x,]<pk2) break; print(tmp[x,])}
        ca[3]<-x
      }  
      #scans pre peak values for studder and if smaller than 20% of actual peak
      #its ignored and slides the start peak position further down until it reaches 
      #>2 or a larger peak is found (likely somewhere else).
      
      tmp[ca[3]:ca[4],]<-0
      ctmp<-sdf(corf(tmp[,],mph))
      #After largest allele called, replaces it with 
      #minpeakheight values
    }
    dth<-0.2
    #resets the dth for next allele
  }
  
  AL<-AL[order(AL[,1],decreasing=F),] 
  #reorders alleles by position
  
  #print(c(id[k],xntv[AL[,1]]))
  #print(AL)
  #above 2 lines prints to screen the pass thru data
    
  pti<-paste(mrkr,'-',as.character(lsid[k-1]))
  #setting pti to pass to alviz
  
  yincr<-setyincr(pmxys[k]) 
  #conditons for changing yincr for plotting
  
  alviz(cfd,cfd[k],mxfr,xntv,pmxys[k],yincr,pti,mph,AL,mrkr,fdd)
  #calls plotting function alviz

  indiv<-as.character(lsid[k-1])
  
  if(AL[1,1]==0) {
    al1<-xntv[AL[2,1]]
    as1<-AL[2,2]
    het<-0
    } else {
    al1<-xntv[AL[1,1]]
    as1<-AL[1,2]
    xev1<-round(as1/minpeakheight[k], digits=2)
    het<-1
    }
  #setting conditions for genotype reporting (i.e., hetero vs homozygotes)
  
  al2<-xntv[AL[2,1]]
  as2<-AL[2,2]
  xev1<-round(as1/minpeakheight[k], digits=2)
  xev2<-round(as2/minpeakheight[k], digits=2)
  scores[(k-1),]<-cbind(mrkr,indiv,al1,al2,as1,as2,xev1,xev2,het,si)
  #initializing scores and binding them in scores df
  }
#above is main loop running through cfd df

scores<-as.data.table(scores)
#converts scores df into data table for easier storage and manipulation

setwd(fdd)

write.table(scores, file=scrfile, append=F, quote=F, sep="\t", eol="\n", col.names=T, row.names=F)
#write.table(cfd, file=cfdfile,append=F, quote=F, sep="\t", eol="\n", col.names=T, row.names=F)

rdafile<-paste(mrkr,".rda",sep="")
save(list=c("xntv","pmxys","mrkr", "mxys","cfd","minpeakheight"), file=rdafile)

print(Sys.time() - starttime)

#time stamp