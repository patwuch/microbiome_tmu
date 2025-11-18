#!/bin/bash
# Created on 04/10/2025
# @author: Joseph, PetSci, petsci.tw
# version: EN2.0
    
ROOT_PATH="/home/patwuch/projects/microbiome/experiments/邱"
DATE_OF_FOLDER="251011"
OUTPUT_PATH="${ROOT_PATH}/${DATE_OF_FOLDER}/04-FunctionAndBiomarker"
COLUMN_OF_GROUP_NAME=("Group")
COLORS=(
  "#19456e"  # Deep Blue
  "#8f3339"  # Dark Red
  "#3e7a4d"  # Forest Green
  "#f2a900"  # Golden Yellow
  "#6c5ce7"  # Vibrant Purple
  "#e17055"  # Coral Orange
  "#00b894"  # Mint Green
  "#fd79a8"  # Pink
  "#d63031"  # Bright Red
  "#0984e3"  # Sky Blue
)



mkdir -p "${OUTPUT_PATH}"

python ${ROOT_PATH}/lefse/plugin_lefse_input.py \
    --table_file "${ROOT_PATH}/qiime/table.qza" \
    --taxonomy_file "${ROOT_PATH}/qiime/taxonomy.qza" \
    --metadata_file "${ROOT_PATH}/metadata.tsv" \
    --output_file "${ROOT_PATH}/8_lefse_input_table.tsv" \
    --class_col "${COLUMN_OF_GROUP_NAME}" 

python ${ROOT_PATH}/lefse/plugin_lefse_format.py \
    "${ROOT_PATH}/8_lefse_input_table.tsv" \
    "${ROOT_PATH}/8_lefse_input_table.in" \
    -c 1 -o 1000000    

python ${ROOT_PATH}/lefse/plugin_lefse_run.py \
   "${ROOT_PATH}/8_lefse_input_table.in" \
   "${ROOT_PATH}/8_lefse_input_table.res" \
    -l 2.0

python ${ROOT_PATH}/lefse/plugin_lefse_barplot.py \
    "${ROOT_PATH}/8_lefse_input_table.res" \
    "${OUTPUT_PATH}/LEfSe LDA.png" \
    --format png \
    --dpi 300 \
    --colors "${COLORS[@]}"

python ${ROOT_PATH}/lefse/plugin_lefse_treeplot.py \
    "${ROOT_PATH}/8_lefse_input_table.res" \
    "${OUTPUT_PATH}/LEfSe Cladogram.svg" \
    --format svg  \
    --dpi 300 \
    --right_space_prop 0.95 \
    --left_space_prop 0.05 \
    --colors "${COLORS[@]}"