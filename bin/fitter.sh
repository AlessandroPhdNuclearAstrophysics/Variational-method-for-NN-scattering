#!/bin/bash

if [ $# -ne 1 ] || [ ! -d "$1" ] || ! ls "$1"/k2_kcotd_*.dat 1>/dev/null 2>&1; then
  echo "Usage: $0 <input_directory_with_k2_kcotd_*.dat_files>"
  exit 1
fi

# Remove trailing slash from input_dir if present
input_dir="${1%/}"
output_dir="fits/$(basename "$input_dir")"
echo "Input directory: $input_dir"
echo "Output directory: $output_dir"

mkdir -p "$output_dir"

for file in "$input_dir"/k2_kcotd_*.dat; do
  [ -e "$file" ] || continue
  base_name=$(basename "$file")
  # Extract everything after k2_kcotd_ prefix
  output_name="${base_name#k2_kcotd_}"
  ./build/main_evaluate_low_energy_scattering_observables.x "$file" > "$output_dir/$output_name"
done