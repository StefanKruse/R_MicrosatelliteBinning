# <<<--------------------------------------------------------- >>>
# <<<--- 
# <<<--- 	R-script to process geneious-files containing
# <<<---	microsatellite scorings with 3 multiplexes 
# <<<--- 
# <<<--- 	5. add additional information 
# <<<--- 	(original individual TreeID and SiteID)
# <<<--- 
# <<<--------------------------------------------------------- >>>


{# <<<--- pre-requisites
	# read data 
		workingdir=paste0(maindir,"/02_processed_data/")
		mplexeinputdatei="msatall_NONBinnedPeaks_Name.csv"
		mplexe=read.csv2(paste(workingdir, mplexeinputdatei, sep="/"), stringsAsFactors=FALSE)
		# change loci names
		names(mplexe)=gsub("local2p2w_","",names(mplexe))

	# read additional data to search in	
		indf=read.csv2(paste0(maindir,"/Treeinformation.csv"), stringsAsFactor=FALSE)
}
	
	
{# <<<--- code for processing
	mplexe$ind=NA
	mplexe$pop=NA
	for(rowi in 1:dim(mplexe)[1])
	{
		mplexe[rowi,]$ind=indf[which(indf$ExtractionID==mplexe[rowi,]$Name),]$TreeID
		mplexe[rowi,]$pop=indf[which(indf$ExtractionID==mplexe[rowi,]$Name),]$SiteID
	}

	# export data
		locis= c("K241" ,"K066","K253","K211","d101" , "K260", "K056" , "K224","Ld45" , "K263" , "K228"  ,"K225", "K235", "K189", "Ld56", "Ld42" )
		mplexeges=mplexe[,c("ind","pop",locis)]
		table(mplexeges$pop)
		write.csv(mplexeges, paste(workingdir,"/mplexeges_allPOPs_",unlist(strsplit(mplexeinputdatei,"_"))[2], ".csv",sep=""), row.names=FALSE)
}

	
if(FALSE)
{# <<<--- code for processing --- additional approach with height and age information for strata separation
	# all separated by population identifier
	# 	# 	agesformat			... 3 numbers giving the age; missing == 999
	#	#						... separator="_"
	#	#	heightformat		... 4 numbers giving the height; missing == 9999
	#	#						... separator="_"
	#	#	ageclass			... 1 number giving the age class (1,2,3,4,5 => <20, <55, <100, <150, >=150 years); missing == 9
	#	#						... separator="_"
	#	#	heightclass			... 4 characters giving the height class (<40 cm == "Seedl", >40 & <200 cm == "Sapling", >200 cm == "Tree"); missing == NNNN
	#	#						... separator="_"
	#	#	heightmeterclasses	... 2 numbers giving the height class in steps by one meter (<1, <2, ...); missing == 99


	# functions for class separation
		getageclass=function(x)
		{
			agei=x
			altersgruppei=NA
			if(is.na(agei)==FALSE)
			{
				if(agei<20)
				{
					altersgruppei=1
				}
				if(agei>=20 & agei<55)
				{
					altersgruppei=2
				}
				if(agei>=55 & agei<100)
				{
					altersgruppei=3
				}
				if(agei>=100 & agei<150)
				{
					altersgruppei=4
				}
				if(agei>=150)
				{
					altersgruppei=5
				}
			} else
			{
				altersgruppei=9
			}
			return(altersgruppei)
		}	
		
		getheightclass=function(x)
		{
			heighti=x
			heightgruppei=NA
			if(is.na(heighti)==FALSE)
			{
				if(heighti<40)
				{
					heightgruppei="Seed"
				}
				if(heighti>=40 & heighti<200)
				{
					heightgruppei="Sapl"
				}
				if(heighti>=200)
				{
					heightgruppei="Tree"
				}
			} else
			{
				heightgruppei="NNNN"
			}
			return(heightgruppei)
		}	
		
		getheightclassinmetersteps=function(x)
		{
			heighti=x
			heightgruppei=NA
			if(is.na(heighti)==FALSE)
			{
				heightgruppei=floor(heighti/100)
			} else
			{
				heightgruppei=99
			}
			return(heightgruppei)
		}	
		
	# reformat contents and reformat height and age into classes
		mplexeges=mplexe[,c("ind","pop",locis)]
		table(mplexeges$pop)

		mplexe$age[is.na(mplexe$age)]=999
		agesformat=formatC(round(mplexe$age,0),width=3, flag=0)
		agesformat=gsub(" NA","999",agesformat)
		
		heightformat=formatC(round(mplexe$height),width=4, flag=0)
		heightformat=gsub("  NA","9999",heightformat)
		
		ageclass=formatC(unlist(lapply(mplexe$age,FUN = function(x) {sapply(x,FUN=getageclass)})),width=1, flag=0)
		heightclass=formatC(unlist(lapply(mplexe$height,FUN = function(x) {sapply(x,FUN=getheightclass)})),width=4, flag=0)
		heightmeterclasses=formatC(unlist(lapply(mplexe$height,FUN = function(x) {sapply(x,FUN=getheightclassinmetersteps)})),width=2, flag=0)

		mplexeges$pop=paste(agesformat,heightformat,ageclass,heightclass,heightmeterclasses,mplexeges$pop,sep="_")

		# export data
		write.csv(mplexeges, paste(workingdir,"/mplexeges_allPOPs_ageheight_",unlist(strsplit(mplexeinputdatei,"_"))[2], ".csv",sep=""), row.names=FALSE)
}

				
				
				
				
				
				
				
				
				
				
				
					
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				