library(tidyr)
library(dplyr)
Pvalue <- 0.01
FoldChange <- 2
PvaluePearson <- 0.05
PearsonCorrelation <- 0.3
BetaDifference <- 0.1

unificati <-
	as.data.frame(
		read.csv(
			file.choose(),
			stringsAsFactors = F,
			sep = ";",
			dec = '.',
			header = T,
			na.strings = c("", "NA"),
			row.names = 1
		)
	)

unificati <- na.omit(unificati)

unificati <- tibble::rownames_to_column(unificati, "sample")

unificati_new <-
	filter(
		unificati,
		(medianUP > medianMedium &
		 	medianMedium > medianDown) |
			(medianUP < medianMedium & medianMedium < medianDown)
	)

unificati_new <-
	filter(unificati_new,
		   (
		   	abs(bd_UpvsMID) >= BetaDifference &
		   		abs(bd_UpvsDOWN) >= BetaDifference &
		   		abs(bd_MIDvsDOWN) >= BetaDifference
		   ))

unificati_new <-
	filter(
		unificati_new,
		(
			pvalue_UpvsMID <= Pvalue &
				pvalue_UpvsDOWN <= Pvalue & pvalue_MIDvsDOWN <= Pvalue
		)
	)

unificati_new <-
	filter(
		unificati_new,
		(
			fc_UpvsMID.gene. >= FoldChange &
				fc_UpvsDOWN.gene. >= FoldChange &
				fc_MIDvsDOWN.gene. >= FoldChange
		)
	)


unificati_new <-
	filter(unificati_new, (pvalue_pearson_correlation <= PvaluePearson))

unificati_new <-
	filter(
		unificati_new,
		(
			pearson_correlation >= PearsonCorrelation |
				pearson_correlation <= -PearsonCorrelation
		)
	)

unificati_new <-
	filter(
		unificati_new,
		(
			pvalue_UpvsMID.gene. <= Pvalue &
				pvalue_UpvsDOWN.gene. <= Pvalue &
				pvalue_MIDvsDOWN.gene. <= Pvalue
		)
	)

write.csv(unificati_new,
		  row.names = FALSE,
		  file = 'CG_of_a_genes(filtrato solo valori).csv')

gene_file_filter <-
	read.csv(
		"/media/giuseppe/DATA/Workspace R/gepia 08.10.18.csv",
		sep = "\t",
		stringsAsFactors = F,
		header = T,
		colClasses = c(NA, "NULL", "NULL")
	)

#filtering of genes
unificati_new <- unificati_new %>%
	separate(sample, c("test", "gene"), "~")

unificati_new <-
	unificati_new[unificati_new$gene %in% gene_file_filter$Gene.Symbol, , drop = FALSE]
unificati_new <-
	unite(unificati_new, "sample", c("test", "gene"), sep = '~')

write.csv(unificati_new,
		  row.names = FALSE,
		  file = 'CG_of_a_genes(filtrato valori e geni).csv')
