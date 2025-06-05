#!/bin/bash

for y in {1..8}; do
  echo "========================================"
  for x in {1..4}; do
    file1="fort.${y}0${x}"
    file2="fort.${y}${y}${x}"

    if [[ -f "$file1" && -f "$file2" ]]; then
      echo "----------------------------------------"
      echo "Comparing: $file1  <->  $file2"
      echo "----------------------------------------"

      vals1=($(cat "$file1"))
      vals2=($(cat "$file2"))

      for i in {0..3}; do
        val1=${vals1[$i]}
        val2=${vals2[$i]}
        # Compare val1 and val2 as numbers, not as strings
        if (( $(echo "$val1 == 0" | bc -l) && $(echo "$val2 == 0" | bc -l) )); then
          continue
        fi

        if [[ $val1 != 0 ]]; then
          awk -v idx=$((i+1)) -v v1="$val1" -v v2="$val2" -v f1="$file1" -v f2="$file2" '
            function abs(x){return x<0?-x:x}
            BEGIN {
              diff = (v1==0 ? 0 : 100*abs(v2-v1)/abs(v1));
              color_start = (diff > 0.001) ? "\033[31m" : "\033[32m";
              color_end = "\033[0m";
              printf("  [%d] = %s\t%s : %s%.18f %% difference%s\n", idx, v1, v2, color_start, diff, color_end);
            }
          '
        else
          printf "  [%d] %s has zero value, cannot compute percentage difference.\n" $((i+1)) "$file1"
        fi
      done
      echo
    else
      echo "Missing file: $file1 or $file2"
    fi
  done
done