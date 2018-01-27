# table with information on loci
msats_colnames=data.frame(rbind(
	c("FAM_real_1","1"	,"bcLK241", NA),
	c("VIC_real_1","1"	,"bcLK066", 190),
	c("VIC_real_1","2"	,"bcLK253", 190),
	c("NED_real_1","1"	,"bcLK211", NA),
	c("ROX_real_1","1"	,"Ld101", NA),
	
	c("FAM_real_2","1"	,"bcLK260", 150),
	c("FAM_real_2","2"	,"bcLK056", 150),
	c("VIC_real_2","1"	,"bcLK224", 205),
	c("VIC_real_2","2"	,"Ld45", 205),
	c("NED_real_2","1"	,"bcLK263", NA),
	c("ROX_real_2","1"	,"bcLK228", NA),
	
	c("FAM_real_3","1"	,"bcLK225", NA),
	c("VIC_real_3","1"	,"bcLK235", NA),
	c("NED_real_3","1"	,"bcLK189", 226.5),	
	c("NED_real_3","2"	,"Ld56", 226.5),
	c("ROX_real_3","1"	,"Ld42", NA)
), stringsAsFactors=FALSE)
names(msats_colnames)=c("Colname","Bin_No","MSat_Name","Boundaries")
msats_colnames$MSat_ID=	substr(msats_colnames$MSat_Name, nchar(msats_colnames$MSat_Name)-3,nchar(msats_colnames$MSat_Name))

# reformat length estimators
	msats_colnameslocal2=msats_colnames
	msats_colnameslocal2$Colname=paste(msats_colnameslocal2$Colname,"local2",sep="_")
	msats_colnameslocal2p1=msats_colnames
	msats_colnameslocal2p1$Colname=paste(msats_colnameslocal2p1$Colname,"local2p1",sep="_")
	msats_colnameslocal2p1w=msats_colnames
	msats_colnameslocal2p1w$Colname=paste(msats_colnameslocal2p1w$Colname,"local2p1w",sep="_")
	msats_colnameslocal2p2w=msats_colnames
	msats_colnameslocal2p2w$Colname=paste(msats_colnameslocal2p2w$Colname,"local2p2w",sep="_")
	msats_colnamesglobal=msats_colnames
	msats_colnamesglobal$Colname=paste(msats_colnamesglobal$Colname,"global",sep="_")
	msats_colnames=rbind(rbind(rbind(msats_colnameslocal2,msats_colnameslocal2p1),rbind(msats_colnameslocal2p1w,msats_colnamesglobal)),msats_colnameslocal2p2w)


