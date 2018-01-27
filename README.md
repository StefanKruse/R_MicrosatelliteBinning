# R Microsatellite Binning
Functions to process scored microsatellite track data exported from Geneious.

Individual steps after preparing the samples in Geneious (check ladder peaks and score alleles in Geneious/microsatellite plugin; export all files from Geneious to one file for each multiplex; extract the individual multiplex files):
1. read scored Geneious-XML-files and transform to data frame
2. read processed geneious data frame and determine fragment lengths
3. read raw fragment lengths and separate into individual microsatellite loci
4. reformat alleles in two alleles per loci
5. add additional information (original individual TreeID and SiteID)
6. binning fragment length into integers
7. import the binned data and convert to genind-object and run first example analyses
Expample data is provided and further details on each step can be found in the R-script.

Currently supported:
- size ladder LIZ600
- microsatellite assay consisting of three multiplexes
  1. bcLK241, bcLK066, bcLK253, bcLK211, Ld101 (N=5)
  2. bcLK260, bcLK056, bcLK224, Ld45, bcLK263, bcLK228 (N=6)
  3. bcLK225, bcLK235, bcLK189, Ld56, Ld42 (N=5)

The multiplexes were developed and applied, and accordingly, these functions compiled/developed during the analyses of microsatellite data for the manuscript: "High gene flow and complex treeline dynamics of Larix Mill. stands on the Taymyr Peninsula (north-central Siberia) revealed by nuclear microsatellites" by Kruse, S., Epp, L. S., Wieczorek, M., Pestryakova, L. A., Stoof-Leichsenring, K. R., Herzschuh, U., accepted for publication in Tree Genetics & Genomes, doi:10.1007/s11295-018-1235-3
