#Prepare_data.R 
#bioinformatics framework using R and R-based bioconductor scripts and tools.
#Written to prepare the raw data for usat scoring by opening all listed files,
#checking that files can be used and preparing a freq.dist to be processed in
#the next autoscr_locus.r script. 

#raw reads ought to be de-multiplexed by individual and locus
#this script recognises either fasta or fastq formats

#Use shift & right click on mouse to copy files to text in textpad
#list should start with "ind"
#dr Aug 2013: updated Jul 2015, updated Feb 2016, Oct 2019, Jul 2022

rm(list=ls())
#re-setting all instances

t1<-Sys.time()

## Only run once to update Bio Conductor packages
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("BiocGenerics")
# BiocManager::install("IRanges")
# BiocManager::install("XVector")
# BiocManager::install("Biostrings")

#library(Defaults)

#loading needed libraries
{ 
	library(zoo)
	library(xts)
	library(TTR)
	library(pracma)
	library(parallel)
	library(BiocGenerics)
	library(IRanges)
	library(XVector)
	library(S4Vectors)
	library(Biostrings)
	library(tcltk)
}

############################ is.0 function #####################################
is.num0 <- function (x) {
return(identical(x,numeric(0)))
}
############################ end is.0  ########################################

############################ fchk -file format check and conversion ###########
fchk<-function(tind) {    #, drc){

trmdr<-dirname(tind)
setwd(trmdr)

exts<-substr(tind,nchar(tind)-4,nchar(tind))
posex<-"fastq"
posex2<-"fasta"

selfile<-readLines(tind)
fixfile<-tind

if(exts==posex) {
    tfile1<-selfile[-seq(4,length(selfile),4)]
  	tfile<-tfile1[-seq(3,length(tfile1),3)]
    tfile<-sub("@",">",tfile,fixed=T)
    tmpn<-basename(tind)
    fixfile<-paste(substr(tmpn,1,nchar(tmpn)-5),"fasta",sep='')
    sink(file=fixfile, append=F)
    writeLines(tfile)
    sink()
 # this is where the file ought to be re-read into the program to start the rest.
  	} else {
	    if(exts!=posex2) warning("only fasta and fastq files are acceptable!")
    } 
    return(fixfile)
 }
############################ end fchk  ########################################

############################ scanal - scanning data to retrieve frequency distribution of alleles ##
scanal<-function(tramp, sat, ampnum)
{
als<-vector(mode='numeric', length=ampnum)
smd<-nchar(sat)

scusat<-gregexpr2(sat,tramp)
difrep<-lapply(scusat,diff)
if (!is.na(locprox)) { 
  screp <-lapply(difrep,rle)
  } else {
cdifrep<-lapply(difrep, function (x) {
      if(is.numeric(x))
      x[x>smd & x<=smd+smd-1]<-smd
      x
      }
      )
screp<-lapply(cdifrep,rle)
}
for (b in 1:ampnum) {
    if(is.num0(as.numeric(screp[[b]]$length))) als[b]<-0 else als[b]<-(max(screp[[b]]$length))+1
}
#als<-als[als > 2]
return(als)
}
############################ end scanal ########################################

lf<-file.choose() 
drc<-dirname(lf)
setwd(drc)
# ask user for list file

indlist<-read.table(lf, header=T)   
attach(indlist)                           
gnum<-length(ind)
# read list into table

pf<-file.choose()
params<-readLines(pf)
# ask user for paramerter file

mrkr<-params[1]
sat<-params[2]
ufr<-params[3]
dfr<-params[4]
locprox<-params[5]

# set parameters from parameter file

for (c in 1:gnum) {                     
  tind<-as.character(ind[c])
  if(!file.exists(tind)) next 
  tind<-fchk(tind)           
  myseqs<-readDNAStringSet(tind, format='fasta')       
  length(myseqs)->seqnum                                           
  if (seqnum <= 3) {     
  #migth be an issue -> if seqnum < 1 then goes to next ind. 
  next 
  } else {
  idn<-as.numeric(gsub("\\D","",basename(tind)))
  idnm<-substr(basename(tind),1,nchar(basename(tind))-6)
  if (idn < 10) idn<-paste(0,0,idn,sep='')
  if (idn < 100) idn<-paste(0,idn,sep='') 
 # frqdt<-matrix(nrow=seqnum,ncol=2)     
  
  tmpseq<-trimLRPatterns(Lpattern=ufr, Rpattern=dfr, subject=myseqs, max.Lmismatch=5, max.Rmismatch=0.25, with.Lindels=T, with.Rindels=T)
 # above removes the flanking regions allowing for some errors
  
  frqdt<-scanal(tmpseq, sat, seqnum)
  # calls scanal function above
  
  if(max(frqdt)<3) next else frqdt<-subset(frqdt,frqdt >= 3)
  
  nfd<-paste(mrkr,'-tpsq',sep='')
  if (file.exists(nfd)) 
     {setwd(file.path(drc,nfd))
     } else {
     dir.create(file.path(drc,nfd), showWarnings = FALSE)
     setwd(file.path(drc,nfd))
     }
  nf<-paste(idnm,"-",mrkr,'.tpsq',sep='')
  writeXStringSet(tmpseq, nf , append=F, format="fasta")   
  setwd(drc)
  # makes a '.tpsq' file storing the truncated amplicons for each individual processed

  frfr<-as.data.frame(table(frqdt))
  # makes the freq dist. of counted repeats and
  # converts to data frame for easier use
  
  names(frfr)[names(frfr)=="frqdt"]<-'alleles'
  #renames frqdt to alleles
  
  fdd<-paste(mrkr,'-fqdt',sep='')

  if (file.exists(fdd)) 
     {setwd(file.path(drc,fdd))
     } else {
     dir.create(file.path(drc,fdd), showWarnings = FALSE)
     setwd(file.path(drc,fdd))
     }
  fd<-paste(idnm,"-",mrkr,'.fqdt',sep='')
  write.table(frfr, file=fd, sep="\t", append=F, quote=F)
  setwd(drc)
  # above cmds writes frequency dist. to an '.fqdt' file
  
  if (max(frfr[,2])<3) next
  
  fdpath<-paste(drc,fdd,fd,sep='/')
  plist<-paste(mrkr,'.fdlist',sep='')
  cat(fdpath,file=plist,append=T, sep='\n')
   
}
}

print(Sys.time()-t1)
# End initial data processing 
