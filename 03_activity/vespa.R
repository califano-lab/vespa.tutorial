library(data.table)
library(plyr)
library(vespa)
library(viper)
library(pheatmap)

set.seed(21)

get_vpmx<-function(osw, spregulon, aregulon, library) {
  # VESPA (substrate-protein)
  spvpl<-viperize(osw, spregulon, xc=TRUE, normalize=TRUE)
  spvpmx<-reshape2::dcast(spvpl,protein_id~tag,value.var="nes",fill=0)
  spvpmx_proteins<-spvpmx$protein_id
  spvpmx$protein_id<-NULL
  spvpmx<-as.matrix(spvpmx)
  row.names(spvpmx)<-spvpmx_proteins
  
  # VESPA (substrate-protein)
  spvpl2<-viperize(osw, spregulon, xc=TRUE, normalize=FALSE)
  spvpmx2<-reshape2::dcast(spvpl2,protein_id~tag,value.var="nes",fill=0)
  spvpmx2_proteins<-spvpmx2$protein_id
  spvpmx2$protein_id<-NULL
  spvpmx2<-as.matrix(spvpmx2)
  row.names(spvpmx2)<-spvpmx2_proteins
  
  # VESPA (activity-protein)
  apvpl<-viperize(vmx2pv(spvpmx2, fasta=library), aregulon, xc=TRUE, normalize=TRUE)
  apvpmx<-reshape2::dcast(apvpl,protein_id~tag,value.var="nes",fill=0)
  apvpmx_proteins<-apvpmx$protein_id
  apvpmx$protein_id<-NULL
  apvpmx<-as.matrix(apvpmx)
  row.names(apvpmx)<-apvpmx_proteins
  
  # integrate substrate and activity-levels
  ipvpl<-ddply(rbind(melt(spvpmx),melt(apvpmx)),.(Var1,Var2),function(X){data.frame("value"=sum(X$value)/sqrt(length(X$value)))})
  ipvpmx<-dcast(ipvpl, Var1~Var2, value.var="value")
  ipvpmx_proteins<-ipvpmx$Var1
  ipvpmx<-as.matrix(ipvpmx[,-1])
  rownames(ipvpmx)<-ipvpmx_proteins
  
  ipvpmx<-ipvpmx[,colnames(apvpmx)]
  
  return(list("substrate"=spvpmx, "activity"=apvpmx, "integrated"=ipvpmx))
}

viperize<-function(osw, regulon, min_size = 10, topn = 500, xc = TRUE, normalize = TRUE){
  # generate quantitative matrix
  qmn<-export2mx(osw, fillvalues = "colmin")
  
  if (normalize) {
    # rank based normalization of the matrix
    qmn.d1 <- t(t(apply(qmn, 2, rank, na.last = "keep"))/(colSums(!is.na(qmn)) + 1))
    
    # rank based estimation of single sample gene expression signature across the matrix
    qmn.norm <- t(apply(qmn.d1, 1, rank, na.last = "keep"))/(rowSums(!is.na(qmn.d1)) + 1)
  } else {
    qmn.norm <- qmn
  }
  
  # run VIPER
  vpres<-viper(qmn.norm, vespa::pruneRegulon(vespa::subsetRegulon(regulon, rownames(qmn.norm), min_size=min_size), cutoff=topn), minsize=min_size, pleiotropy = xc, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, targets = 5, penalty = 20, method = "adaptive"), cores=6)
  
  # transform viper results
  vpd<-as.data.frame(vpres)
  vpd$protein_id<-row.names(vpd)
  vpl<-reshape2::melt(vpd, id=c("protein_id"), variable.name = "tag", value.name = "nes")

  return(vpl)
}

# load quantitative data
cptac<-readRDS("../01_import/CPTAC_S046_LUAD_phospho.rds")

# get FASTA library
library<-"../01_import/library.fasta"

# import substrate-protein regulon
spregulon<-readRDS("../02_network/stdpimeta_substrate_protein_regulon.rds")

# import activity-protein regulon
aregulon<-readRDS("../02_network/dpimeta_activity_protein_regulon.rds")

# run VESPA
vpmx<-get_vpmx(cptac, spregulon, aregulon, library)

# generate heatmap plots
legend_max<-max(abs(c(max(vpmx$substrate), min(vpmx$substrate), max(vpmx$activity), min(vpmx$activity), max(vpmx$integrated), min(vpmx$integrated))))
legend_breaks<-c(seq(-legend_max, -2.5, length.out=10),seq(-2,2,0.05),seq(2.5, legend_max, length.out=10))

pheatmap(vpmx$substrate, breaks=legend_breaks, show_rownames=FALSE, show_colnames=FALSE, filename = "vespa_substrate.pdf", width=10, height=10)
pheatmap(vpmx$activity, breaks=legend_breaks, show_rownames=FALSE, show_colnames=FALSE, filename = "vespa_activity.pdf", width=10, height=10)
pheatmap(vpmx$integrated, breaks=legend_breaks, show_rownames=FALSE, show_colnames=FALSE, filename = "vespa_integrated.pdf", width=10, height=10)

# store VESPA results in matrix
write.csv(vpmx$substrate, file="vespa_substrate.csv")
write.csv(vpmx$activity, file="vespa_activity.csv")
write.csv(vpmx$integrated, file="vespa_integrated.csv")
