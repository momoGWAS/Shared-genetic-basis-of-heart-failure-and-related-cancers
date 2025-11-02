library(data.table)
library(tidyverse)
all=fread("MGI_all.txt")
all=all[,1:5]
colnames(all)=c("GENE","ENTREZ","mouse mark symbol","MGI","phenotype")
result <- all %>%
  group_by(GENE) %>%
  summarize(
    MGI = paste(unique(MGI), collapse = ","),
    phenotype = paste(unique(phenotype), collapse = ",")
  )
##mgi数据库匹配
gene=fread("gene.txt")
gene=gene[!duplicated(gene$ID),]
colnames(gene)="GENE"
phe=left_join(gene,result,by="GENE")
##提取有表型的gene
sig=filter(phe,phe$phenotype!="")
data_split <- sig %>%
  separate_rows(phenotype, sep = ",")
data_split =filter(data_split ,data_split$phenotype!="")
data_split$phenotype=trimws(data_split$phenotype)##去掉空格
bb=as.data.frame(table(data_split$phenotype))##显著性组表型数
###非显著性组表型数
all_split <- result %>%
  separate_rows(phenotype, sep = ",")
all_split$phenotype= trimws(all_split$phenotype)
final_split=filter(all_split,all_split$phenotype!="")
aa=as.data.frame(table(final_split$phenotype))

merge=left_join(bb,aa,by="Var1")
###创建四格表，fisher检验
merge=merge |> mutate(fre10=43-Freq.x,
                      Freq.y.1=Freq.y-Freq.x,
                      fre00=19325-43-Freq.y.1)
merge=merge[,c(1,2,5,4,6)]
out <- data.frame()
for (i in 1:nrow(merge)){
  t <- fisher.test(matrix(as.vector(t(merge[i, 2:5])), ncol=2),alternative="greater")
  
  d <- merge[i, ]
  d$p.value <- t$p.value
  out <- rbind(out, d)
}
colnames(out)=c("phenotype","11","01","10","00","p")

##把基因嵌套在后面
res_gene <- data_split %>%
  group_by(phenotype) %>%
  summarize(
    GENE = paste(unique(GENE), collapse = ",")
  )
final=merge(out,res_gene,by="phenotype")

write.csv(final,file="MGI_fisher.csv")
