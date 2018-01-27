# <<<--------------------------------------------------------- >>>
# <<<--- 
# <<<--- 	R-script to process geneious-files containing
# <<<---	microsatellite scorings with 3 multiplexes 
# <<<--- 
# <<<--- 	2. read processed geneious data frame and determine
# <<<--- 	fragment lengths
# <<<--- 
# <<<--------------------------------------------------------- >>>

# check ...ERRORLOG_not36ladderpeaks text file in 02_processed_data/ folder and check if ladder peaks were correctly assigned

{# <<<--- pre-requisites
	# ladder peaks of LIZ600
	liz600peaks_all=c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600 )

	# read in table with raw data and determine fragment lengths
	workingdir_in=paste0(maindir,"/01_geneious_data/")
	workingdir_out=paste0(maindir,"/02_processed_data/")
	dateinamen=paste("msatdf_geneious_",c(paste0("MP")),"123.csv",sep="")

}

{# <<<--- code for processing
	for(dateinamei in dateinamen)
	{
		# progress in console
			print(paste(which(dateinamen==dateinamei),"/",length(dateinamen), " == ",dateinamei,sep=""))
			
		# read data
			msatdf_geneious=read.csv2(paste(workingdir_in,dateinamei,sep=""),stringsAsFactors=FALSE)
			msatdf_geneious[,paste(paste(c("FAM_real","VIC_real","NED_real","ROX_real"),sort(rep(1:3,4)),sep="_"),sort(rep(c("local2","local2p1","local2p1w","local2p2w","global"),4*3)), sep="_")]=paste(c(NA,NA),collapse="/")
			# ... for each sample
			for(indii in 1:dim(msatdf_geneious)[1])
			{
				# ... for each multiplex
				for(mplexi in 1:3)
				{
					xladder =as.numeric(unlist(strsplit(as.character(msatdf_geneious[indii,paste("ladderpeakscollapsed",mplexi,sep="_")]),"/")))
					# print error log if more or less than 36 ladder peaks are present
					nxladder=length(xladder)
					if(nxladder!=36)
					{
						sink(paste(workingdir_out,dateinamei,"ERRORLOG_not36ladderpeaks.txt",sep=""),append=TRUE)
							print(paste("mplexno=",mplexi,"nxladder=",nxladder,msatdf_geneious[indii,]))
						sink()
					}
					if(is.na(xladder)==FALSE)
					{
						yladder=rev(rev(liz600peaks_all)[1:nxladder])
						
						# calculate peak lengths with different strategies
						# ... local2
							lmsteiglocal2=NULL
							lmynulllocal2=NULL
							for(i in 1:(nxladder-1))
							{
								lmlocal2=lm(yladder[i:(i+1)]~xladder[i:(i+1)])
								# points(yladder[i:(i+1)]~xladder[i:(i+1)], col=colvec[i])
								# abline(lmlocal2, col=colvec[i], lwd=3)
								lmynulllocal2=c(lmynulllocal2, as.numeric(lmlocal2$coefficients[1]))
								lmsteiglocal2=c(lmsteiglocal2, as.numeric(lmlocal2$coefficients[2]))
							}
						# ... local2+1
							lmsteiglocal2p1=NULL
							lmynulllocal2p1=NULL
							for(i in 2:(nxladder-1))
							{
								lmlocal2p1=lm(yladder[c(i-1,i,i+1)]~xladder[c(i-1,i,i+1)])
								# points(yladder[c(i-1,i,i+1)]~xladder[c(i-1,i,i+1)], col=colvec[i])
								# abline(lmlocal2p1, col=colvec[i], lty=2)
								lmynulllocal2p1=c(lmynulllocal2p1, as.numeric(lmlocal2p1$coefficients[1]))
								lmsteiglocal2p1=c(lmsteiglocal2p1, as.numeric(lmlocal2p1$coefficients[2]))
							}
						# ... local2+1 + weights
							lmsteiglocal2p1w=NULL
							lmynulllocal2p1w=NULL
							for(i in 2:(nxladder-1))
							{
								lmlocal2p1w=lm(yladder[c(i-1,i,i+1)]~xladder[c(i-1,i,i+1)], weights=c(0.5,1,0.5))
								# points(yladder[c(i-1,i,i+1)]~xladder[c(i-1,i,i+1)], col=colvec[i])
								# abline(lmlocal2p1w, col=colvec[i], lty=3)
								lmynulllocal2p1w=c(lmynulllocal2p1w, as.numeric(lmlocal2p1w$coefficients[1]))
								lmsteiglocal2p1w=c(lmsteiglocal2p1w, as.numeric(lmlocal2p1w$coefficients[2]))
							}
						# ... local2+2 + weights
							lmsteiglocal2p2w=NULL
							lmynulllocal2p2w=NULL
							for(i in (1+2):(nxladder-2))
							{
								lmlocal2p2w=lm(yladder[c(i-2,i-1,i,i+1,i+2)]~xladder[c(i-2,i-1,i,i+1,i+2)], weights=c(0.25,0.5,1,0.5,0.25))
								lmynulllocal2p2w=c(lmynulllocal2p2w, as.numeric(lmlocal2p2w$coefficients[1]))
								lmsteiglocal2p2w=c(lmsteiglocal2p2w, as.numeric(lmlocal2p2w$coefficients[2]))
							}
						# ... global
							lmglob=lm(yladder~xladder)
							# abline(lmglob)
							lmglobnull=as.numeric(lmglob$coefficients[1])
							lmglobsteig=as.numeric(lmglob$coefficients[2])

						
						# add new columns for each length determination
						for(interpolarti in c("local2","local2p1","local2p1w","local2p2w","global"))
						{
							for(farbei in c("FAM", "VIC", "NED", "ROX"))
							{
								if(msatdf_geneious[indii,paste(farbei,mplexi,sep="_")]!="NA/NA")
								{
									toestimate=as.numeric(unlist(strsplit(as.character(msatdf_geneious[indii,paste(farbei,mplexi,sep="_")]),"/")))
									estiloc2=estiloc2p1=estiloc2p1w=estiloc2p2w=estiglob=NULL
									for(toestimatei in toestimate)
									{
										nearestval=which.min(abs(xladder-toestimatei))
										#local2p1 + 2p1w
										# ... nearest
										estiloc2p1=c(estiloc2p1, toestimatei*lmsteiglocal2p1[nearestval]+lmynulllocal2p1[nearestval])
										estiloc2p1w=c(estiloc2p1w, toestimatei*lmsteiglocal2p1w[nearestval]+lmynulllocal2p1w[nearestval])
										estiloc2p2w=c(estiloc2p2w, toestimatei*lmsteiglocal2p2w[nearestval]+lmynulllocal2p2w[nearestval])

										#local2
										# ... next value
										if((xladder[nearestval]-toestimatei)>0)
										{
											nearestval=nearestval-1
										}
										estiloc2=c(estiloc2, toestimatei*lmsteiglocal2[nearestval]+lmynulllocal2[nearestval])
									}
									estiglob=c(estiglob, toestimate*lmglobsteig+lmglobnull)

									if(interpolarti=="local2")
									{
										vals=estiloc2
									} else if(interpolarti=="local2p1")
									{
										vals=estiloc2p1
									} else if(interpolarti=="local2p1w")
									{
										vals=estiloc2p1w
									} else if(interpolarti=="local2p2w")
									{
										vals=estiloc2p2w
									} else if(interpolarti=="global")
									{
										vals=estiglob
									}
								
									msatdf_geneious[indii,paste(farbei,"real",mplexi,interpolarti,sep="_")]=paste(vals	,collapse="/")
								}
							}
						}
					}
				}
			}
		
		# export data
			write.csv2(msatdf_geneious,paste(workingdir_out,dateinamei,"processed.csv",sep=""),row.names=FALSE)
	}
}

	

	
	
	
	
	
	
	

