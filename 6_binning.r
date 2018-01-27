# <<<--------------------------------------------------------- >>>
# <<<--- 
# <<<--- 	R-script to process geneious-files containing
# <<<---	microsatellite scorings with 3 multiplexes 
# <<<--- 
# <<<--- 	6. binning fragment length into integers
# <<<--- 
# <<<--------------------------------------------------------- >>>

{# <<<--- pre-requisites
	# read data
		filein="mplexeges_allPOPs_NONBinnedPeaks.csv"
		workingdir=paste0(maindir,"/02_processed_data/")
		mplexe=read.csv(paste(workingdir, filein, sep="/"), stringsAsFactors=FALSE)
		str(mplexe)
		# change names
		names(mplexe)=gsub("local2p2w_","",names(mplexe))
}

{# <<<--- code for processing
	# function
		splitandsort = function(x)
		{
			if(x!="NA/NA")
			{
				allels=as.numeric(unlist(strsplit(x, "/")))
				for(alleli in unique(allels))
				{
					allabs=abs(start-alleli)
					mini=which(allabs==min(allabs))
					allabs[mini]=NA
					mini2=which(na.omit(allabs)==min(na.omit(allabs)))
					
					allels[which(allels==alleli)]=namen_bins[min(c(mini,mini2))]
				}
				
				return(paste(allels,collapse="/"))
			} else
			{
				return("NA/NA")
			}
		}
			

	# bin all peaks and export a diagnostic plot		
		pdf(paste(workingdir, "msatall_NONBinnedPeaks_Name_ALLELs.pdf", sep="/"),width=30,height=10)
			
			# parameter for binning
			maxbreak=600
			
			# K241
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K241"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(140,210))
				start=seq(140,210,2)
					abline(v=start, col="red")
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(start)[-1])
						namen_bins=(namen_bins+start[-1])/2
					# 2. bin peaks
						mplexe[,"K241"]=sapply(mplexe[,"K241"], splitandsort)
					
			# K066
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K066"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K066"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(150,185))
					startbreak=152.7
					binsraw=seq(startbreak,startbreak+200,2)
					abline(v=binsraw)
					start=binsraw+(binsraw-startbreak)*0.05
					abline(v=start, col="red")
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"K066"]=sapply(mplexe[,"K066"], splitandsort)

			# K253
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K253"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K253"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(205,250))
					startbreak=215.9
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					start=binsraw-(binsraw-startbreak)*0.02
					abline(v=start, col="red")
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"K253"]=sapply(mplexe[,"K253"], splitandsort)


			# K211
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K211"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K211"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(190,250))
				start=seq(101,401,2)
				abline(v=start, col="red")
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(start)[-1])
						namen_bins=(namen_bins+start[-1])/2
					# 2. bin peaks
						mplexe[,"K211"]=sapply(mplexe[,"K211"], splitandsort)
					

			# d101
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"d101"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"d101"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(190,234))
					startbreak=200.3
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					start=binsraw-(binsraw-startbreak)*0.04
					abline(v=start, col="red")
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"d101"]=sapply(mplexe[,"d101"], splitandsort)


			# K260
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K260"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K260"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(100,160))
					startbreak=107.5
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.07
					factorstretch_down=0.0175
					start=c(seq(startbreak-200,startbreak-2,2)+(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)+(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"K260"]=sapply(mplexe[,"K260"], splitandsort)

						
			# K056
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K056"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K056"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(150,255))
					startbreak=162.1
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.02
					factorstretch_down=0.02
					start=c(seq(startbreak-200,startbreak-2,2)+(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)-(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"K056"]=sapply(mplexe[,"K056"], splitandsort)
					
			# K224
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K224"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K224"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(140,160))
					startbreak=138.8
					start=seq(startbreak,startbreak+200,2)
					abline(v=start, col="red")
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(start)[-1])
						namen_bins=round((namen_bins+start[-1])/2)
					# 2. bin peaks
						mplexe[,"K224"]=sapply(mplexe[,"K224"], splitandsort)
					
			# Ld45
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"Ld45"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"Ld45"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(210,260))
					startbreak=219.0
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.02
					factorstretch_down=0.02
					start=c(seq(startbreak-200,startbreak-2,2)-(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)-(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"Ld45"]=sapply(mplexe[,"Ld45"], splitandsort)
					
			# K263
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K263"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K263"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(195,280))
					startbreak=200.2
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.022
					factorstretch_down=0.022
					start=c(seq(startbreak-200,startbreak-2,2)+(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)+(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"K263"]=sapply(mplexe[,"K263"], splitandsort)
						
						
			# K228
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K228"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K228"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(195,280))
					startbreak=201.6
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.03
					factorstretch_down=0.03
					start=c(seq(startbreak-200,startbreak-2,2)-(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)-(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"K228"]=sapply(mplexe[,"K228"], splitandsort)
					
			# K225
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K225"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K225"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(170,225))
					startbreak=170
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.01
					factorstretch_down=0.01
					start=c(seq(startbreak-200,startbreak-2,2)-(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)-(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=ceiling(namen_bins)
					# 2. bin peaks
						mplexe[,"K225"]=sapply(mplexe[,"K225"], splitandsort)
					
			# K235
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K235"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K235"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(135,275))
					startbreak=212.09
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.025
					factorstretch_down=0.01
					start=c(seq(startbreak-200,startbreak-2,2)+(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)+((seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up+c(0, .5, .4, rep(0,length((seq(startbreak,startbreak+200,2)-startbreak))-3))))
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=floor(namen_bins)
					# 2. bin peaks
						mplexe[,"K235"]=sapply(mplexe[,"K235"], splitandsort)
					
			# K189
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"K189"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"K189"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(150,245))
					startbreak=165
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.01
					factorstretch_down=0
					start=c(seq(startbreak-200,startbreak-2,2)+(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)-(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=floor(namen_bins)
					# 2. bin peaks
						mplexe[,"K189"]=sapply(mplexe[,"K189"], splitandsort)
					
			# Ld56
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"Ld56"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"Ld56"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(225,280))
					startbreak=260.5
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.0
					factorstretch_down=0
					start=c(seq(startbreak-200,startbreak-2,2)+(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)-(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=floor(namen_bins)
					# 2. bin peaks
						mplexe[,"Ld56"]=sapply(mplexe[,"Ld56"], splitandsort)
				
			# Ld42
			range(na.omit(as.numeric(unlist(strsplit(mplexe[,"Ld42"], "/")))))
				plot(hist(na.omit(as.numeric(unlist(strsplit(mplexe[,"Ld42"], "/")))), breaks=c(seq(0,maxbreak,.1)), plot=FALSE), xlim=c(170,205))
					startbreak=190.5
					binsraw=seq(startbreak-200,startbreak+200,2)
					abline(v=binsraw)
					factorstretch_up=0.0
					factorstretch_down=0.05
					start=c(seq(startbreak-200,startbreak-2,2)-(seq(startbreak+200,startbreak+2,-2)-startbreak)*factorstretch_down,seq(startbreak,startbreak+200,2)-(seq(startbreak,startbreak+200,2)-startbreak)*factorstretch_up)
					abline(v=start, col="red");abline(v=startbreak, lwd=3, col="red", lty=2)
					# 1. determine names for the bins as the mean of the two limits of the bin
						namen_bins=rev(rev(binsraw)[-1])
						namen_bins=(namen_bins+binsraw[-1])/2
						namen_bins=floor(namen_bins)
					# 2. bin peaks
						mplexe[,"Ld42"]=sapply(mplexe[,"Ld42"], splitandsort)
		dev.off()
			
	# export data
		write.csv(mplexe, paste(workingdir, "mplexeges_allPOPs_ageheight_BINNEDPeakslocal2p2w_.csv", sep="/"),row.names=FALSE)
		write.csv2(mplexe, paste(workingdir, "mplexeges_allPOPs_ageheight_BINNEDPeakslocal2p2w_csv2.csv", sep="/"),row.names=FALSE)
}		
	
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
	