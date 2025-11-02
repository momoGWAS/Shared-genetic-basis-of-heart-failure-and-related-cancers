
install.packages(remotes)
library(remotes)
remotes::install_github("YuLab-SMU/createKEGGdb")
createKEGGdb::create_kegg_db("hsa")
library("KEGG.db")
library(KEGG.db)
library(clusterProfiler)
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1,use_internal_data =T)
kk_res <- kk@result
dotplot(kk, showCategory = 20)

