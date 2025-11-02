
input_folder="./GWAS/"
input_files=(${input_folder}*.txt)


for input_file in "${input_files[@]}"; do

    filename=$(basename -- "${input_file}")
    filename_no_ext="${filename%.*}"

    output_file="./result/${filename_no_ext}"

    ./smr-1.3.1 --bfile ./g1000_eur/g1000_eur --gwas-summary "${input_file}" --beqtl-summary ./eQTL_besd_lite/Adipose_Subcutaneous.lite --diff-freq-prop 0.90 --out "${output_file}"_Adipose_Subcutaneous.lite --thread-num 30
done
