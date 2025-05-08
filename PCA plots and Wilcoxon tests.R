# PCA plots and Wilcoxon tests
# Author: Yue Xing

pca_sig_plot = function(d1, ct, ra, wt, e_type, c_type) {

	library(ggplot2)
	library(vegan)
	library(rgl)
	library(ggbiplot)

	colnames(d1)=c("ID","Type")
	d1=d1[d1[,1]%in%colnames(ct),]
	ct=ct[,colnames(ct)%in%d1[,1]]
	ra=ra[,colnames(ra)%in%d1[,1]]
	a=which(rowSums(ra)==0)
	if (length(a)!=0) {
		ra=ra[-a,]
		ct=ct[-a,]
	}
	ra=ra[,d1[,1]]
	ct=ct[,d1[,1]]

	rownames(ct)=gsub(".*s__","",rownames(ct))
	rownames(ra)=gsub(".*s__","",rownames(ra))

	# only use the taxa with p<0.05
	wt=wt[wt$p_value<0.05,]
	wt=wt[!is.na(wt$p_value),]
	wt=wt$Taxon
	wt=wt[!duplicated(wt)]
	print(length(wt))
	if (length(wt)==1) {
		return("taxa length = 1!")
	}
	ct=ct[wt,]
	ra=ra[wt,]

	print(all(colnames(ra)==colnames(ct)))
	print(all(rownames(ra)==rownames(ct)))
	print(all(colnames(ct)==d1[,1]))

	coll=d1$Type
	coll=gsub("HV","green",coll)
	coll=gsub("IBS_D","blue",coll)
	coll=gsub("IBS_C","red",coll)

	# pca
	tryCatch(
		expr={

			ct=t(ct)
	
			ct.pca <- prcomp(ct, center = TRUE,scale. = TRUE)
			#summary(ct.pca)
	
			png(paste0("Plots/sig_taxa_only030925/",e_type,".",c_type,".pca.png"),width=2000,height=2000,res=300)

			print(ggbiplot(ct.pca,ellipse=TRUE,var.axes=FALSE,labels=rownames(ct),groups=as.factor(d1$Type))+theme_bw()+ scale_color_manual(values=c("green","blue","red")))

			dev.off()
	
			if (ncol(ct)>500) {
				ct2=ct
				rv <- rowMeans(ra)
				ntop=500
				select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
				ct2=ct2[,select]

				ct2.pca <- prcomp(ct2, center = TRUE,scale. = TRUE)
		
				png(paste0("Plots/sig_taxa_only030925/",e_type,".",c_type,".top500RA.pca.png"),width=2000,height=2000,res=300)

				print(ggbiplot(ct2.pca,ellipse=TRUE,var.axes=FALSE,labels=rownames(ct2),groups=as.factor(d1$Type))+theme_bw()+ scale_color_manual(values=c("green","blue","red")))

				dev.off()

				print(dim(ct))
				print(dim(ct2))
			}
		},
		error = function(e){print(paste("pca",e_type,c_type))}
	)
}
	
options(stringsAsFactors=FALSE)

d10=read.csv("All_on_list.csv")
ct0=read.csv("MP3_genus_count_all.ALR.csv",row.names=1)
ra0=read.csv("MP3_genus_ra_all.csv",row.names=1)
wt0=read.csv("MP3_genus.pairwise.wilcox.csv")
wt0$p_value=as.numeric(as.character(wt0$p_value))
pca_sig_plot(d10,ct0,ra0,wt0,nm,nm2)
