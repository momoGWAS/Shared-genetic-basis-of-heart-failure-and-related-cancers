#!/bin/bash

# FUMA文件夹的路径
FUMA_DIR="./FUMA"

# 结果文件夹的路径
OUT_DIR="./结果"

# 创建结果文件夹如果它不存在的话
mkdir -p "$OUT_DIR"

# 遍历FUMA文件夹下的每个子文件夹
for SUBDIR in "$FUMA_DIR"/*; do
  if [ -d "$SUBDIR" ]; then
    # 获取子文件夹的名字
    SUBDIR_NAME=$(basename "$SUBDIR")
    # 构建输出前缀
    OUT_PREFIX="$OUT_DIR/$SUBDIR_NAME"
    # 运行pops.py脚本
    python pops.py --gene_annot_path ./gene_all.txt \
                   --feature_mat_prefix ./features_munged/pops_features \
                   --num_feature_chunks 2 \
                   --magma_prefix "$SUBDIR/magma" \
                   --control_features_path ./features_jul17_control.txt \
                   --out_prefix "$OUT_PREFIX"
  fi
done
