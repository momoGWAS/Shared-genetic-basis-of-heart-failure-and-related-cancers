#!/bin/bash

# 指定脑功能文件夹下的所有txt文件
input_folder="./GWAS/"
input_files=(${input_folder}*.txt)

# 循环处理每个txt文件
for input_file in "${input_files[@]}"; do
    # 获取当前txt文件的文件名（不含扩展名）
    filename=$(basename -- "${input_file}")
    filename_no_ext="${filename%.*}"

    # 构建输出文件名
    output_file="./结果/${filename_no_ext}"

    # 执行smr命令
    ./smr-1.3.1 --bfile ./g1000_eur/g1000_eur --gwas-summary "${input_file}" --beqtl-summary ./eQTL_besd_lite/Adipose_Subcutaneous.lite --diff-freq-prop 0.90 --out "${output_file}"_Adipose_Subcutaneous.lite --thread-num 30
done
