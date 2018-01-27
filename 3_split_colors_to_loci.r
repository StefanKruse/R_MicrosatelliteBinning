# <<<--------------------------------------------------------- >>>
# <<<--- 
# <<<--- 	R-script to process geneious-files containing
# <<<---	microsatellite scorings with three multiplexes 
# <<<--- 
# <<<--- 	3. read raw fragment lengths and separate into
# <<<--- 	individual microsatellite loci
# <<<--- 
# <<<--------------------------------------------------------- >>>


{# <<<--- pre-requisites
	
	# read additional information on microsatellite loci separation
		source(paste0(maindir,"/3_split_colors_to_loci_RSOURCEFILE.r"))
		str(msats_colnames)

	
	# read data
		workingdir=paste0(maindir,"/02_processed_data/")
		dateinamen=list.files(workingdir, pattern="123.csvprocessed.csv")

		for(dateinameni in dateinamen)
		{
			print(paste(which(dateinamen==dateinameni),"/",length(dateinamen), " == ",dateinameni,sep=""))
			
				# read data
				msatdf_geneious_NNN=read.csv2(paste(workingdir,dateinameni,sep=""),stringsAsFactors=FALSE)
				
				# merge all data
				if(!exists("msatdf_geneious"))
				{
					msatdf_geneious=msatdf_geneious_NNN
				} else 
				{
					msatdf_geneious=rbind(msatdf_geneious,msatdf_geneious_NNN)
				}	
		}
		str(msatdf_geneious)

}

{# <<<--- code for processing	

	# reformat peaks sorted by fluorescent dyes to individual microsatellite loci
		for(colnamei in msats_colnames$Colname)
		{
			boundariesi=msats_colnames[msats_colnames$Colname==colnamei,"Boundaries"]
				
			outputcolname=paste(as.character(unlist((as.data.frame(strsplit(msats_colnames[msats_colnames$Colname==colnamei,]$Colname,"_")))[4,])),msats_colnames[msats_colnames$Colname==colnamei,]$MSat_ID, sep="_")
			
			# copy contents if only one microsatellite loci of actual fluorescent dye
			if(is.na(boundariesi)==TRUE)
			{
				msatdf_geneious[,outputcolname]=msatdf_geneious[,colnamei]
			} else
			{
				fragsi=msatdf_geneious[,colnamei]
				kleinerfragsi=grosserfragsi=NULL
				for(fragsii in fragsi)
				{
					if(fragsii=="NA/NA" | fragsii=="NA/NA/NA/NA")
					{
						kleinerfragsi=c(kleinerfragsi,fragsii)
						grosserfragsi=c(grosserfragsi,fragsii)
					} else 
					{
						fragsii_split=as.numeric(unlist(strsplit(fragsii,"/")))
						# for peaks smaller than the boundary
						wertkleiner=paste(fragsii_split[which(fragsii_split < as.numeric(unique(boundariesi)))],collapse="/")
						if(wertkleiner=="") {wertkleiner="NA/NA"}
						kleinerfragsi=c(kleinerfragsi, wertkleiner)
						# for peaks greater than the boundary
						wertgroeszer=paste(fragsii_split[which(fragsii_split >= as.numeric(unique(boundariesi)))],collapse="/")
						if(wertgroeszer=="") {wertgroeszer="NA/NA"}
						grosserfragsi=c(grosserfragsi, wertgroeszer)
					}
				}
				msatdf_geneious[,outputcolname[1]]=kleinerfragsi
				msatdf_geneious[,outputcolname[2]]=grosserfragsi
			}
		}
	



	# manually convert microsatellite loci for two primer of which products enter the range of others
	# ... before that, check in Geneious for 
	# ... ... A multiplex 3 in fluorescent dye NED (loci Ld56 into bcLK189)
	# ... ... B multiplex 2 in fluorescent dye 6-FAM (loci bcLK056 into bcLK260)
	{# A
		# read list with names and remaining peaks in the upper microsatellite loci
			manedit=read.csv2(paste0(maindir,"/mp3_ned_bcLK189_Ld56.csv"), stringsAsFactors=FALSE)
			str(manedit)
			submanedit=manedit[manedit$bp>226.5,]
			
		# move individual peaks
			for(zubearbeitenzeilei in 1:dim(submanedit)[1])
			{
				zeilennum=which(msatdf_geneious$name==submanedit[zubearbeitenzeilei,"eskname"])
				remaininLd56=submanedit[zubearbeitenzeilei,"peaksremainLd56"]
				
				# delete in Ld56 and add to bcLK189
				buf0=msatdf_geneious[zeilennum,grep("Ld56",names(msatdf_geneious))]
				# buf0
				buf1=msatdf_geneious[zeilennum,grep("K189",names(msatdf_geneious))]
				# buf1
				
				if(dim(buf0)[1]==1) # proceed only if information are available
				{
					for(veki in 1:5)
					{
						alleallele=c(as.character(unlist(strsplit(buf1[,veki],"/"))),as.character(unlist(strsplit(buf0[,veki],"/"))))
						
						# print before and after move of peaks to console
						print(msatdf_geneious[zeilennum,"name"])
						alleallele=alleallele[alleallele!="NA"]
						print(alleallele)
						
						# redistribute peaks to loci
						alK189=alleallele[1:(length(alleallele)-remaininLd56)]
						if(remaininLd56!=0)
						{
							alLd56=alleallele[(length(alleallele)-remaininLd56+1):length(alleallele)]
						}else
						{
							alLd56=c("NA","NA")
						}
						
						buf0[,veki]=paste(alLd56,collapse="/")
						buf1[,veki]=paste(alK189,collapse="/")
					}
					
					# change data in original data frame
					msatdf_geneious[zeilennum,grep("Ld56",names(msatdf_geneious))]=buf0
					msatdf_geneious[zeilennum,grep("K189",names(msatdf_geneious))]=buf1
				} else 
				{
					print(paste(submanedit[zubearbeitenzeilei,"eskname"]," not found in data frame"))
				}

				# delete individual buffer container
				rm(buf0)
				rm(buf1)
			}
	}# 

	{# B
		# read list with names and remaining peaks in the upper microsatellite loci
			manedit=read.csv2(paste0(maindir,"/mp2_fam_bcLK260_bcLK056.csv"), stringsAsFactors=FALSE)
			str(manedit)
			submanedit=manedit[manedit$bp>150,]
			
		# move individual peaks
			for(zubearbeitenzeilei in 1:dim(submanedit)[1])
			{
				zeilennum=which(msatdf_geneious$name==submanedit[zubearbeitenzeilei,"eskname"])
				remaininK056=submanedit[zubearbeitenzeilei,"peaksremainK056"]
				
				# delete in K056 and add to bcLK260
				buf0=msatdf_geneious[zeilennum,grep("K056",names(msatdf_geneious))]
				# buf0
				buf1=msatdf_geneious[zeilennum,grep("K260",names(msatdf_geneious))]
				# buf1
				
				if(dim(buf0)[1]==1) # proceed only if information are available
				{
					for(veki in 1:5)
					{
						alleallele=c(as.character(unlist(strsplit(buf1[,veki],"/"))),as.character(unlist(strsplit(buf0[,veki],"/"))))
						
						# print before and after move of peaks to console
						print(msatdf_geneious[zeilennum,"name"])
						alleallele=alleallele[alleallele!="NA"]
						print(alleallele)
						
						# redistribute peaks to loci
						alK260=alleallele[1:(length(alleallele)-remaininK056)]
						if(remaininK056!=0)
						{
							alK056=alleallele[(length(alleallele)-remaininK056+1):length(alleallele)]
						}else
						{
							alK056=c("NA","NA")
						}
						
						buf0[,veki]=paste(alK056,collapse="/")
						buf1[,veki]=paste(alK260,collapse="/")

					}
					
					# change data in original data frame
					msatdf_geneious[zeilennum,grep("K056",names(msatdf_geneious))]=buf0
					msatdf_geneious[zeilennum,grep("K260",names(msatdf_geneious))]=buf1
				} else 
				{
					print(paste(submanedit[zubearbeitenzeilei,"eskname"]," not found in data frame"))
				}
				
				# delete individual buffer container
				rm(buf0)
				rm(buf1)
			}
	}
	
	# export data
		write.csv2(msatdf_geneious, paste(workingdir, "msatdf_geneious_realFragSizes.csv", sep="/"),row.names=FALSE)
	
}	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	