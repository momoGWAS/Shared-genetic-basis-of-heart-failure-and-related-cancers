#!/bin/bash


FUMA_DIR="./FUMA"
OUT_DIR="./result"
mkdir -p "$OUT_DIR"


for SUBDIR in "$FUMA_DIR"/*; do
  if [ -d "$SUBDIR" ]; then
   
    SUBDIR_NAME=$(basename "$SUBDIR")
 
    OUT_PREFIX="$OUT_DIR/$SUBDIR_NAME"
 
    python pops.py --gene_annot_path ./gene_all.txt \
                   --feature_mat_prefix ./features_munged/pops_features \
                   --num_feature_chunks 2 \
                   --magma_prefix "$SUBDIR/magma" \
                   --control_features_path ./features_jul17_control.txt \
                   --out_prefix "$OUT_PREFIX"
  fi
done
