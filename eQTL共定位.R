if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("SummarizedExperiment") # For windows or linux
BiocManager::install("PupillometryR") # For windows or linux




setwd('C:/Users/dr_li/Desktop/HF')
# Install from github to get the latest version.
if(!require("devtools")){install.packages("devtools")}
devtools::install_github("dingruofan/xQTLbiolinks")

library(data.table)
library(xQTLbiolinks)
library(stringr)
library(coloc)


data=fread("G:/GWAS/心衰/MTAG/result/mtag后-MHC.txt")
data <- data[, c("SNP", "CHR","BP","mtag_pval","FRQ","mtag_beta","mtag_se")]
colnames(data) <- c("rsid", "chrom","position","pValue","AF","beta","se")

sentinelSnpDF <- xQTLanalyze_getSentinelSnp(data, centerRange=1e6,
                                            genomeVersion="grch37", grch37To38=TRUE)
tissueAll <- xQTLquery_tissue()
tissueSiteDetail="Skin - Not Sun Exposed (Suprapubic)"

#Artery - Aorta
#Artery - Coronary
#Artery - Tibial
#Uterus
#Lung
#Heart - Atrial Appendage
#Heart - Left Ventricle
#Skin - Not Sun Exposed (Suprapubic)
#Skin - Sun Exposed (Lower leg)

traitsAll <- xQTLanalyze_getTraits(sentinelSnpDF, detectRange=1e6, tissueSiteDetail=tissueSiteDetail)
#output <- xQTLanalyze_coloc(data, "ENSG00000130203", tissueSiteDetail=tissueSiteDetail) # using gene symbol
gencodeId<-c('ENSG00000130203','ENSG00000065361')
outputs <- rbindlist(lapply( unique(gencodeId), function(x){ # using gencode ID.
  xQTLanalyze_coloc(data, x, tissueSiteDetail=tissueSiteDetail, method = "Both")$colocOut }))
write.csv(outputs, file = "Skin - Not Sun Exposed (Suprapubic).csv", row.names = FALSE)









#xQTLvisual_eqtl("ENSG00000065361")

eqtlAsso <- xQTLdownload_eqtlAllAsso(gene="ENSG00000130203", 
                                     tissueLabel = tissueSiteDetail, data_source = "liLab")
gwasEqtldata <- merge(data, eqtlAsso, by="rsid", suffixes = c(".gwas",".eqtl"))
xQTLvisual_locusCompare(gwasEqtldata[,.(rsid, pValue.eqtl)], 
                        gwasEqtldata[,.(rsid, pValue.gwas)], legend_position = "bottomright")



