setwd("G:/GWAS/心衰/COLOC")
library(data.table)
library(dplyr)
library(coloc)

# 读取 SNP.txt 文件
snp_data <- fread("SNP.txt")

# 获取唯一的 Group 列表
unique_groups <- unique(snp_data$Group)

# 遍历每个 Group
for (group_val in unique_groups) {
  # 提取 Group 列的前后缀部分
  parts <- unlist(strsplit(group_val, "-"))
  group1 <- parts[1]
  group2 <- parts[2]
  
  # 读取 GWAS 数据集
  data1 <- fread(paste0("G:/GWAS/心衰/COLOC/GWAS/", group1, ".txt"))
  data2 <- fread(paste0("G:/GWAS/心衰/COLOC/GWAS/", group2, ".txt"))
  data1$VAR <- data1$SE^2
  data2$VAR <- data2$SE^2
  
  # 筛选属于当前 Group 的 SNP
  group_snp_data <- snp_data %>% filter(Group == group_val)
  
  # 初始化结果列表
  result_list <- list()
  result_res_list <- list()
  
  # 遍历当前 Group 中的每一行
  for (i in 1:nrow(group_snp_data)) {
    # 获取当前行的 CHR 和 BP 值
    chr_val <- group_snp_data$CHR[i]
    start_val <- group_snp_data$BP[i]
    
    # 在当前 CHR 和 BP 值的范围内过滤数据
    data3 <- data1 %>% filter(CHR == chr_val, BP >= (start_val - 500000), BP <= (start_val + 500000))
    data4 <- data2 %>% filter(CHR == chr_val, BP >= (start_val - 500000), BP <= (start_val + 500000))
    
    # 检查是否有重叠的 SNP
    if (nrow(data3) == 0 || nrow(data4) == 0) {
      cat(paste("No overlapping SNPs for CHR", chr_val, "Group", group_val, "\n"))
      # 在结果列表中添加 "NA" 元素
      result_list[[i]] <- data.frame(CHR = chr_val, BP = NA, SNP = NA, coloc.post.prob = NA)
      result_res_list[[i]] <- data.frame(SNP = NA, SNP.PP.H4 = NA)
      next  # 跳过当前循环，进行下一个 SNP 的处理
    }
    
    # 合并过滤后的数据
    data <- merge(data3, data4, by = "SNP")
    data <- data[!duplicated(data$SNP), ]
    data <- data %>% mutate(BETA.y = ifelse(A1.x == A1.y, BETA.y, -BETA.y))
    
    data3 <- data[, c("BETA.x", "VAR.x", "SNP")]
    data4 <- data[, c("BETA.y", "VAR.y", "SNP")]
    colnames(data3) <- c("beta", "varbeta", "snp")
    colnames(data4) <- c("beta", "varbeta", "snp")
    data3 <- as.list(data3)
    data4 <- as.list(data4)
    data3$type <- "cc"
    data4$type <- "cc"
    
    # 运行 coloc.abf
    res <- coloc.abf(data3, data4, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    
    # 提取结果的 summary
    result_list[[i]] <- res$summary
    
    # 提取 SNP.PP.H4 最大的行，随机选择一个
    max_h4_rows <- res$results %>% filter(SNP.PP.H4 == max(SNP.PP.H4, na.rm = TRUE))
    max_h4_row <- max_h4_rows[sample(nrow(max_h4_rows), 1), ]
    result_res_list[[i]] <- max_h4_row
  }
  
  # 将结果合并为一个数据框
  result_df <- bind_rows(result_list)
  result_res_df <- bind_rows(result_res_list)
  
  # 查看结果
  print(result_df)
  print(result_res_df)
  
  # 保存结果，文件名以 Group 名称为前缀
  write.csv(result_df, paste0("结果/", group_val, "_summary.csv"), row.names = FALSE)
  write.csv(result_res_df, paste0("SNP/", group_val, "_max_SNP.PP.H4.csv"), row.names = FALSE)
}
