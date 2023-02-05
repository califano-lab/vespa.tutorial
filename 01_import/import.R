library(vespa)
library(data.table)

# import CCT data
phospho<-importPhosphoCPTAC("CPTAC3_Lung_Adeno_Carcinoma_Phosphoproteome.phosphopeptide.tmt10.tsv.gz", fasta="library.fasta", cores=6)
proteo<-importProteoCPTAC("CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv.gz", fasta="library.fasta")

# exclude any non-tumor data
phospho_tumor<-subset(phospho, run_id %in% fread("sample_ids.txt")$sample_id,sep="")
proteo_tumor<-subset(proteo, run_id %in% fread("sample_ids.txt")$sample_id,sep="")

# save processed data in RDS file
saveRDS(subset(phospho_tumor, run_id %in% proteo_tumor$run_id), file="CPTAC_S046_LUAD_phospho.rds")
saveRDS(subset(proteo_tumor, run_id %in% phospho_tumor$run_id), file="CPTAC_S046_LUAD_proteo.rds")
