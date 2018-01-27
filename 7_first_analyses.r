# <<<--------------------------------------------------------- >>>
# <<<--- 
# <<<--- 	R-script to process geneious-files containing
# <<<---	microsatellite scorings with 3 multiplexes 
# <<<--- 
# <<<--- 	7. first example analyses
# <<<--- 
# <<<--------------------------------------------------------- >>>

{# <<<--- pre-requisites
	library(poppr)

	# read data
		mplexe=read.csv2(paste0(maindir,"/02_processed_data/", "mplexeges_allPOPs_ageheight_BINNEDPeakslocal2p2w_csv2.csv"),stringsAsFactors=FALSE)
		# order by pop
		mplexe=mplexe[order(mplexe$pop),]
		str(mplexe)
		
	# .. transform to genind-object
		dfi=df2genind(mplexe[,3:18], ploidy=2, sep="/", pop=mplexe$pop, ind=make.unique(mplexe$ind), NA.char = "NA")
		dfi
		summary(dfi)
		
	# subset of only selected loci (see details in Kruse et al. submitted to Tree Genes & Genomes 2017)
		loc_8=locNames(dfi)[which(!locNames(dfi)%in%c("K241", "Ld56", "K225", "K224", "Ld45", "K235", "K260", "K066"))]
		dfi_lessmissing=missingno(dfi[loc=loc_8], type="genotype", cutoff=0.2)
		dfi_lessmissing
}

{# <<<--- analyse data
	# number of samples
	png(paste0(maindir,"/03_analyses/samplenumber.png"))
		par(mar=c(4,15,2,2), las=1)
		barplot(table(pop(dfi_lessmissing)), horiz=TRUE, xlab="Number of samples", main="")
		abline(v=10,lty=2)
	dev.off()

	png(paste0(maindir,"/03_analyses/genotype_curve.png"))
		genotype_curve(dfi_lessmissing, sample = 1000, quiet = TRUE)	
	dev.off()
	
	# ordination
	png(paste0(maindir,"/03_analyses/PCA.png"))
		dpc=dudi.pca(missingno(dfi_lessmissing, type="mean"), scannf=FALSE)
		s.class(dpc$li, pop(dfi_lessmissing), clab=1, col=rainbow(length(levels(pop(dfi_lessmissing)))))
	dev.off()

	# DAPC
	png(paste0(maindir,"/03_analyses/DAPC.png"))
		dapc1=dapc(dfi_lessmissing, n.pca=30,n.da=2)
		scatter(dapc1)
	dev.off()
	
	png(paste0(maindir,"/03_analyses/DAPC_compoplot.png"), width=3000, height=400)
		compoplot(dapc1)
	dev.off()
}
	
