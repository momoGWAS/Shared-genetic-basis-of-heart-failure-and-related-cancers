#安装Y叔的包，
#安装创建KEGG数据库的包的包
install.packages(remotes)
library(remotes)
remotes::install_github("YuLab-SMU/createKEGGdb")
#创建自己的物种的包create_kegg_db，会自动创建名称为KEGG.db_1.0.tar,gz的包。物种名称的简写，在
createKEGGdb::create_kegg_db("hsa")
library("KEGG.db")
library(KEGG.db)
library(clusterProfiler)
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1,use_internal_data =T) #设置为T就是本地数据
kk_res <- kk@result
dotplot(kk, showCategory = 20)
