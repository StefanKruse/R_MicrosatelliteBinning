# <<<--------------------------------------------------------- >>>
# <<<--- 
# <<<--- 	R-script to process geneious-files containing
# <<<---	microsatellite scorings with 3 multiplexes 
# <<<--- 
# <<<--- 	1. read scored geneious XML-files and transform
# <<<--- 	to data frame
# <<<--- 
# <<<--------------------------------------------------------- >>>


{# <<<--- pre-requisites
	# export all files from Geneious to one file for each multiplex
	# ... e.g. file input name for each multiplex MP1.geneious
	# sample names must meet the format: "COORDINATE_NAME_MULTIPLEX..."
	# ... first, before split "_" coordinate in plate, 
	# ... second, name of sample, DNA-extraction
	# ... third, the multiplex number after a split character "_" of choice
	# ... set the split character between sample name and multiplex ID
	samplenamempsplit_name_mplex="-"

	# unzip to folder ./MP1/
	# ... content must contain all files "fileData.[0 running until max file number]" and one XML file MP1.geneious

	# set the working directory folder path
	workingdir_in=paste0(maindir,"/01_geneious_data/")

	# set the file base name
	# ... single folders ending with "1", "2" and "3" with this base name will be read in
	dateinamen=c(paste0("MP")) 
}

{# <<<--- code for processing
	library(XML)

	# declaration of function
	readGeneiousXML=function(dateinamevollstaendig, samplenamempsplit="_") 
	{
		xmldata=xmlParse(dateinamevollstaendig)
		str(xmldata)

		anzahltraces=xmlSize(getNodeSet(xmldata, "//geneiousDocument"))

		ladderpeakscollapsed =NULL
		for(i in 1:anzahltraces)
		{
			anzahlladderpeaks=xmlSize(getNodeSet(xmldata, "//TracesData/ladderTrace/peaks")[[i]])
			ladderpeaks=NULL
			for(ladpeaki in 1: anzahlladderpeaks)
			{
				ladderpeaks =c(ladderpeaks, xmlValue(getNodeSet(xmldata, "//TracesData/ladderTrace/peaks")[[i]][[ladpeaki]]))
			}

			ladderpeakscollapsed=c(ladderpeakscollapsed, paste(as.character(unlist(ladderpeaks)), collapse="/"))
		}

		tracedateiname=as.character(xmldata["//originalSourceFile/@xmlFileData"])
		FAM=VIC=NED=ROX=NULL
		tubename=NULL
		for(i in 1:anzahltraces)
		{
			if(is.na(as.numeric(xmlValue(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]][[5]])))==FALSE)
			{
				tbnamei=xmlValue(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]][[2]])
				if(nchar(tbnamei)==0)
				{
					tbnamei=xmlValue(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]][[6]])
				}
				if(nchar(tbnamei)==0)
				{
					print(paste("Error => no name for location ", i, " in plate ", dateinamevollstaendig, " found"))
				}
			} else
			{
				if(xmlSize(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]])==6)
				{
					tbnamei=xmlValue(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]][[5]])
				}
				if(xmlSize(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]])==7)
				{
					tbnamei=xmlValue(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]][[6]])
				}
				if(xmlSize(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]])==8)
				{
					tbnamei=xmlValue(getNodeSet(xmldata, "//geneiousDocument")[[i]][[1]][[6]])
				}
			}

			tubename=c(tubename, tbnamei)

			fam=vic=ned=rox=NULL
			for(coli in 1:4)
			{
				speicherfragmente=NULL
				
				anzahlmsattraces=xmlSize(getNodeSet(xmldata, "//TracesData/microsatelliteTraces/microsatelliteTrace")[[(i*4-4)+coli]][[1]])
				for(peaki in 1:anzahlmsattraces)
				{
					if(peaki==0)
					{
						speicherfragmente=c(speicherfragmente,NA)
					} else
					{
						speicherfragmente=c(speicherfragmente,xmlValue(getNodeSet(xmldata, "//TracesData/microsatelliteTraces/microsatelliteTrace")[[(i*4-4)+coli]][[1]][[peaki]]))
					}
				}
				if(coli==1)
				{
					fam=speicherfragmente
				}
				if(coli==2)
				{
					vic=speicherfragmente
				}
				if(coli==3)
				{
					ned=speicherfragmente
				}
				if(coli==4)
				{
					rox=speicherfragmente
				}
				
			}
			FAM=c(FAM,paste(fam, collapse="/"))
			VIC=c(VIC,paste(vic, collapse ="/"))
			NED=c(NED,paste(ned, collapse ="/"))
			ROX=c(ROX,paste(rox, collapse ="/"))
		}

		msatdf_geneious=data.frame(tubename,tracedateiname,FAM,VIC,NED,ROX,ladderpeakscollapsed, stringsAsFactors=FALSE)
		
		msatdf_geneious$name=unlist(strsplit(split=samplenamempsplit,as.character(t(as.data.frame(strsplit(as.character(msatdf_geneious$tubename),"_")))[,2])))[c(TRUE,FALSE)]
		msatdf_geneious$coordinates=as.character(t(as.data.frame(strsplit(as.character(msatdf_geneious$tubename),"_")))[,1])
		
		foldername=rev(unlist(strsplit(dateinamevollstaendig,"/")))[2]
		mplexi=substr(foldername,start=nchar(foldername),stop=nchar(foldername))
		
		names(msatdf_geneious)[-which(names(msatdf_geneious)=="name")]=paste(names(msatdf_geneious)[-which(names(msatdf_geneious)=="name")],mplexi,sep="_")
		
		return(msatdf_geneious)
	}



	# analysis of all multiplexes
	for(dateinameni in dateinamen)
	{
		# print progress in console
			print(paste(which(dateinamen==dateinameni),"/",length(dateinamen), " == ",dateinameni,sep=""))
			Sys.sleep(.5)

		# read and transform data
			dateinamevollstaendig_alle=paste(paste(paste(workingdir_in,paste(dateinameni,1:3,sep=""), sep="/"),paste(dateinameni,1:3,sep=""),sep="/"),".geneious",sep="")
			m1=readGeneiousXML(dateinamevollstaendig_alle[1], samplenamempsplit=samplenamempsplit_name_mplex)
			m2=readGeneiousXML(dateinamevollstaendig_alle[2], samplenamempsplit=samplenamempsplit_name_mplex)
			m3=readGeneiousXML(dateinamevollstaendig_alle[3], samplenamempsplit=samplenamempsplit_name_mplex)	
			merge123=merge(merge(m1, m2, all=TRUE), m3, all=TRUE)
		
		# export data
			write.csv2(merge123, paste(workingdir_in, "/msatdf_geneious_",dateinameni,"123.csv",sep=""), row.names=FALSE)
	}


}

