# <<<--------------------------------------------------------- >>>
# <<<--- 
# <<<--- 	R-script to process geneious-files containing
# <<<---	microsatellite scorings with 3 multiplexes 
# <<<--- 
# <<<--- 	4. reformat alleles in two alleles per loci
# <<<--- 
# <<<--------------------------------------------------------- >>>

# check the console output for >2 alleles "not2buffer"

{# <<<--- pre-requisites
	# read data
		workingdir=paste0(maindir,"/02_processed_data/")
		msatdf_geneious=read.csv2(paste(workingdir,"msatdf_geneious_realFragSizes.csv",sep=""),stringsAsFactors=FALSE)
		str(msatdf_geneious)

		# create subset containing name, microsatellite loci
		source(paste0(maindir,"/3_split_colors_to_loci_RSOURCEFILE.r"))
		msatdf_geneious_msats_name=msatdf_geneious[,c("name",paste(as.character(unlist((as.data.frame(strsplit(msats_colnames$Colname,"_")))[4,])),msats_colnames$MSat_ID, sep="_"))]
		str(msatdf_geneious_msats_name)

		# ... subset for fragment length determination method of "local2p2w"
		msatdf_geneious_msats_name=msatdf_geneious_msats_name[,c(1,grep("local2p2w",names(msatdf_geneious_msats_name)))]
		
		# ... export data
		write.csv2(msatdf_geneious_msats_name, paste(workingdir, "msatall_RAWPeaks_Name.csv", sep=""), row.names=FALSE)

		# read the final data frame
		msatdf_geneious_msats_name=read.csv2(paste(workingdir, "msatall_RAWPeaks_Name.csv", sep=""), stringsAsFactors=FALSE)
		str(msatdf_geneious_msats_name)
}

{# <<<--- code for processing	

	# function for reformatting the alleles of each cell in the data frame
	# ... split alleles by "/"
	# ... assure that only two elements per cell are present, if not, copy!
	# ... generally all peaks <50 bp are excluded prior processing
	reformatpeaksinCells <- function(x) 
	{
		lengthsi=unlist(strsplit(x,split="/"))
	
		# change to NUMBER/NUMBER or NA/NA
		if(lengthsi[1]!="NA")
		{
			# check if below the threshold
			numslengthsi=as.numeric(lengthsi)
			lengthsi=numslengthsi[!numslengthsi<50]
			if(length(lengthsi)==3)
			{
				return(paste(lengthsi,collapse="/"))
			}
			if(length(lengthsi)==2)
			{
				return(paste(lengthsi,collapse="/"))
			}
			if(length(lengthsi)==1)
			{
				return(paste(c(lengthsi,lengthsi),collapse="/"))
			}	
			if(length(lengthsi)==0)
			{
				return("NA/NA")
			}				
		} else
		{
			return("NA/NA")
		}
	}

	
	# reformat peaks in each cell
		msatdf_geneious_msats_name_reformatted=as.data.frame(lapply(msatdf_geneious_msats_name[,-1],FUN = function(x) {sapply(x,FUN=reformatpeaksinCells)}), stringsAsFactors=FALSE)
		msatdf_geneious_msats_name_reformatted$Name=msatdf_geneious_msats_name$name
		str(msatdf_geneious_msats_name_reformatted)
	
	# check if >2 allels are present in any of the microsatellite loci
		columnbuffer=not2buffer=NULL
		for(spaltei in 1:16)
		{
			print(names(msatdf_geneious_msats_name_reformatted)[spaltei])
			nicht2=which(sapply(strsplit(msatdf_geneious_msats_name_reformatted[,spaltei],split="/"), length)!=2)
			if(length(nicht2)>0)
			{
				# print(paste(spaltei, nicht2, msatdf_geneious_msats_name_reformatted[spaltei,nicht2]))
				print(paste(spaltei, nicht2, msatdf_geneious_msats_name_reformatted[nicht2,spaltei]))
				columnbuffer=c(columnbuffer, rep(spaltei, length(nicht2)))
				not2buffer=c(not2buffer, nicht2)
			}
		}
		columnbuffer
			names(msatdf_geneious_msats_name_reformatted)[columnbuffer]
		not2buffer
			msatdf_geneious_msats_name_reformatted[not2buffer,"Name"]
		cbind(msatdf_geneious_msats_name_reformatted[not2buffer,"Name"],msatdf_geneious_msats_name_reformatted[not2buffer,columnbuffer])
		msatdf_geneious_msats_name_reformatted[not2buffer,columnbuffer+1]
		# ... if here individuals are printed in the console, please check the given sample for >2 peaks and rerun the whole scripts again
		
	
	# export data	
		write.csv2(msatdf_geneious_msats_name_reformatted, paste(workingdir, "msatall_NONBinnedPeaks_Name.csv", sep=""), row.names=FALSE)

}
				
				
				
				
























