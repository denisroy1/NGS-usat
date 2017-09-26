# NGS-usat 1.0

### An R based bioinformatics framework for scoring microsatellites generated from next generation sequencing platforms.

## 1 – System requirements & suggestions: (what has worked for us during development)

**1.1 –** NGS-usat scripts were written using R version 3.4.0, but have since been updated to run with latest R (3.4.1). See instructions from the CRAN website (https://cran.r-project.org/) on how to download and install R, depending on your platform. All packages and dependencies for the scripts have also been updated to their most recent versions and these can be accessed through your favourite CRAN mirror (we suggest the one at Dal). Instruction on how to download and install packages from the CRAN mirrors, and from the BIOCONDUCTOR sources are also available online (see https://www.r-bloggers.com/installing-r-packages/ and https://www.bioconductor.org/install/, respectively).

**1.2 –** Much of the processing and manipulation for NGS-usat is accomplished rather easily using the RStudio framework (https://www.rstudio.com/products/rstudio/download/). This is not essential, but it does make visualisation and updating the scores much easier. The version of RStudio used here is 1.0.136. 

**1.3 –** Mac OSX users need to have the latest XQuartz application loaded properly in their root system. Otherwise updating the scores will not load on your machines. Instructions on how to do this, and the application itself, are available at the link provided (https://www.xquartz.org/). 

* PLEASE NOTE! RStudio works on Mac/Windows/Unix/Linux platforms, but the latest Mac and Windows version is 64bit only. This means the JAVA version on Windows has to be the 64bit application which is available here (https://java.com/en/download/faq/java_win64bit.xml#Java%20for%2064-bit). On Macs the newest OSX (Sierra 10.12.2) still doesn’t incorporate the XQuartz.app. So, even if you end up using RStudio to run the NGS-usat scripts, the updated XQuartz app is still needed.

**1.4 –** For list creation and manipulation, we recommend using Textpad® in Windows. However, any basic text editor will do, so long as it doesn’t introduce silent/unseen characters (MS WORD, not a good candidate). Mac OSX users can find the equivalent (or better) text manipulation using TextWrangler®. Both are free and easily installed.

## 2 – NGS-usat Introduction and Setup:

**2.1 –** NGS-usat is a series of scripts written in R that will allow automatic and fast genotyping of individuals sequenced at a series of microsatellites using NGS techniques. The scripts are essentially platform independent and so, ought to work with both Ion Torrent and Illumina sequence reads, as long as the raw data for each individual and locus can be demultiplexed first. Data demultiplexing can be performed using other freely available software packages (we recommend the RDpipeline available here: http://rdp.cme.msu.edu/).

The easiest way to run the scripts is through the RStudio interface. Once RStudio is loaded and running, open the scripts in the ‘sourcing’ tab (Figure 1). Highlight the script in its entirety and press the enter/return key while holding down the CTRL/CMD key (in Windows or Mac, respectively). This will activate the script(s) and will prompt users for appropriate files. 

PLEASE NOTE! Most scripts prompt the user for at least one file. So, if it looks as if R has stalled, it may very likely be that there is a dialog box waiting for user input. It might be wise to leave part of the desktop available to see ‘file choose’ dialog boxes. 





