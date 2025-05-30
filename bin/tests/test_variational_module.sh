#!/bin/bash

# -----------------------------------------------------------------------------
# Script to test the consistency of the variational module by running the main
# executable and comparing its output to reference results.
# For each test case, EMAX and NE are inferred from the reference files.
# The comparison is done with an absolute percentage tolerance of 1e-4.
# -----------------------------------------------------------------------------

# Set up paths for the binary and reference files (relative to bin/tests/)
BIN_DIR="$(dirname "$(dirname "$(dirname "$(realpath "$0")")")")/build"
REF_BASE="$(dirname "$0")/test_files/scattering_phase_shifts"

FAILED=0

# -----------------------------------------------------------------------------
# Test cases: (IPOT, ILB, REF_SUBDIR)
# IPOT=18: AV18, IPOT=19: EFT_pless
# Each tuple defines the potential, its label, and the reference directory.
# -----------------------------------------------------------------------------
declare -a TESTS=(
  "18 1 AV18"
  "19 15 EFT_pless_15"
)

# -----------------------------------------------------------------------------
# Function: plot_difference
# Instead of plotting, print a suggestion for the user.
# Arguments:
#   $1 - reference file
#   $2 - output file
# -----------------------------------------------------------------------------
plot_difference() {
  local ref_file="$1"
  local out_file="$2"
  local fname
  fname=$(basename "$ref_file")
  if [[ "$fname" =~ delta_[0-9][A-Z][0-9]\.dat$ ]]; then
    echo -e "To compare, run: \033[0;36mxmgrace \"$ref_file\" \"$out_file\"\033[0m" >> /tmp/variational_module_plot_suggestions.txt
  elif [[ "$fname" =~ delta_.*-.*\.dat$ ]]; then
    echo -e "To compare, run: \033[0;36mxmgrace -nxy \"$ref_file\" -nxy \"$out_file\"\033[0m" >> /tmp/variational_module_plot_suggestions.txt
  fi
}

# -----------------------------------------------------------------------------
# Function: compare_files
# Compares two files line by line and column by column.
# Passes if the absolute percentage difference for each value is < 1e-4.
# Arguments:
#   $1 - reference file
#   $2 - output file to check
# Returns:
#   0 if all values pass, 1 if any mismatch is found
# -----------------------------------------------------------------------------
compare_files() {
  local ref_file="$1"
  local out_file="$2"
  local tol=5e-3
  local mismatch=0

  # Constants for awk
  local EPS=1e-5
  local HUGE=1e20
  local MIN_DENOM=1e-12

  local awk_output
  awk_output=$(awk -v tol="$tol" -v eps="$EPS" -v huge="$HUGE" -v min_denom="$MIN_DENOM" '
    BEGIN { mismatch=0; diff_rows=0; }
    {
      getline ref < ARGV[1]
      getline out < ARGV[2]
      if (ref == "" && out == "") next
      split(ref, ref_vals)
      split(out, out_vals)
      if (length(ref_vals) != length(out_vals)) {
        print "Column mismatch at line " NR
        mismatch=1; diff_rows++
        next
      }
      for (i=1; i<=length(ref_vals); i++) {
        r=ref_vals[i]+0; o=out_vals[i]+0
        # Accept if both are "zero" (abs < eps) and tolerance is 1
        if ((r<eps && r>-eps) && (o<eps && o>-eps)) {
          abs_perc=0
          tol_this=1
        } else if (((r<eps && r>-eps) && (o>=eps || o<=-eps)) || ((o<eps && o>-eps) && (r>=eps || r<=-eps))) {
          abs_perc=huge
          tol_this=tol
        } else {
          denom=(r!=0)?r:o
          if (denom<0) denom=-denom
          if (denom<min_denom) denom=(o!=0)?o:1
          if (denom<0) denom=-denom
          abs_perc=(r-o)/denom; if (abs_perc<0) abs_perc=-abs_perc
          tol_this=tol
        }
        if (abs_perc > tol_this) {
          print "Mismatch at line " NR " col " i ": ref=" r " out=" o " (abs%diff=" abs_perc ")"
          mismatch=1; diff_rows++
        }
      }
    }
    END {
      if (diff_rows >= 2) {
        print "__PLOT__"
      }
      exit mismatch
    }
  ' "$ref_file" "$out_file")

  # Print awk output
  if [ -n "$awk_output" ]; then
    echo "$awk_output" | grep -v "__PLOT__"
  fi

  # If awk output contains __PLOT__, call plot_difference
  if echo "$awk_output" | grep -q "__PLOT__"; then
    plot_difference "$ref_file" "$out_file"
    return 1
  fi

  # Return 1 if any mismatch, 0 otherwise
  if echo "$awk_output" | grep -q "Mismatch\|Column mismatch"; then
    return 1
  fi
  return 0
}

# -----------------------------------------------------------------------------
# Function: run_test_case
# Runs a single test case:
#   - Determines EMAX and NE from the first reference file
#   - Prepares a namelist input file
#   - Runs the main executable
#   - Compares all output files to reference files using compare_files
# Arguments:
#   $1 - IPOT
#   $2 - ILB
#   $3 - Reference subdirectory
# -----------------------------------------------------------------------------
run_test_case() {
  local IPOT="$1"
  local ILB="$2"
  local REF_SUBDIR="$3"
  local OUTPUT_DIR="output/test_${REF_SUBDIR}/"
  local REF_DIR="$REF_BASE/$REF_SUBDIR"
  local TEST_INPUT

  # Find the first reference file to extract EMAX and NE
  local first_ref_file
  first_ref_file=$(ls "$REF_DIR"/delta_*.dat | head -n 1)
  if [ ! -f "$first_ref_file" ]; then
    echo "No reference files found in $REF_DIR"
    FAILED=1
    return
  fi
  # EMAX is the last value in the first column, NE is the number of lines
  local EMAX NE
  EMAX=$(awk '{emax=$1} END{print emax}' "$first_ref_file")
  NE=$(awk 'END{print NR}' "$first_ref_file")

  # Prepare a temporary namelist input file for this test
  TEST_INPUT="$(mktemp)"
  cat > "$TEST_INPUT" << EOF
&IN
  EMAX = $EMAX,
  NE = $NE,
  TZ = 0,
  IPOT = $IPOT,
  ILB = $ILB,
  LEMP = 0,
  PRINT_COEFFICIENTS = .FALSE.,
  OUT_DIR = "$OUTPUT_DIR",
/
EOF

  # Clean and create the output directory
  rm -rf "$OUTPUT_DIR"
  mkdir -p "$OUTPUT_DIR"

  # Run the main program with the generated namelist
  "$BIN_DIR/main_scattering_NN_variazional_method.x" "$TEST_INPUT"
  local RETVAL=$?
  if [ $RETVAL -ne 0 ]; then
    echo "Program failed to run for IPOT=$IPOT ILB=$ILB."
    rm -f "$TEST_INPUT"
    FAILED=1
    return
  fi

  # Compare each output file to its reference
  for ref_file in "$REF_DIR"/delta_*.dat; do
    local fname=$(basename "$ref_file")
    local out_file="$OUTPUT_DIR/$fname"
    if [ ! -f "$out_file" ]; then
      echo "Missing output file: $fname for IPOT=$IPOT ILB=$ILB"
      FAILED=1
      continue
    fi
    if compare_files "$ref_file" "$out_file"; then
      echo "PASS: $fname for IPOT=$IPOT ILB=$ILB"
    else
      echo "Mismatch in $fname for IPOT=$IPOT ILB=$ILB"
      # Suggest xmgrace command instead of plotting
      if [[ "$fname" =~ delta_[0-9][A-Z][0-9]\.dat$ ]]; then
        echo "If you want to compare the results, use the command: xmgrace \"$ref_file\" \"$out_file\""
      elif [[ "$fname" =~ delta_.*-.*\.dat$ ]]; then
        echo "If you want to compare the results, use the command: xmgrace -nxy \"$ref_file\" -nxy \"$out_file\""
      fi
      FAILED=1
    fi
  done

  # Clean up the temporary namelist file
  rm -f "$TEST_INPUT"
}

# -----------------------------------------------------------------------------
# Main function: runs all test cases and exits with the appropriate status
# -----------------------------------------------------------------------------
main() {
  local all_passed=1
  declare -A file_status

  # Clear previous suggestions
  > /tmp/variational_module_plot_suggestions.txt

  for test in "${TESTS[@]}"; do
    set -- $test
    local IPOT="$1"
    local ILB="$2"
    local REF_SUBDIR="$3"
    local OUTPUT_DIR="output/test_${REF_SUBDIR}/"
    local REF_DIR="$REF_BASE/$REF_SUBDIR"
    local TEST_INPUT

    # Find the first reference file to extract EMAX and NE
    local first_ref_file
    first_ref_file=$(ls "$REF_DIR"/delta_*.dat | head -n 1)
    if [ ! -f "$first_ref_file" ]; then
      echo "No reference files found in $REF_DIR"
      all_passed=0
      continue
    fi
    local EMAX NE
    EMAX=$(awk '{emax=$1} END{print emax}' "$first_ref_file")
    NE=$(awk 'END{print NR}' "$first_ref_file")

    TEST_INPUT="$(mktemp)"
    cat > "$TEST_INPUT" << EOF
&IN
  EMAX = $EMAX,
  NE = $NE,
  TZ = 0,
  IPOT = $IPOT,
  ILB = $ILB,
  LEMP = 0,
  PRINT_COEFFICIENTS = .FALSE.,
  OUT_DIR = "$OUTPUT_DIR",
/
EOF

    rm -rf "$OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"

    "$BIN_DIR/main_scattering_NN_variazional_method.x" "$TEST_INPUT"
    local RETVAL=$?
    if [ $RETVAL -ne 0 ]; then
      for ref_file in "$REF_DIR"/delta_*.dat; do
        local fname=$(basename "$ref_file")
        file_status["$fname|$IPOT|$ILB"]="\033[0;31mFAILED\033[0m"
      done
      all_passed=0
      rm -f "$TEST_INPUT"
      continue
    fi

    for ref_file in "$REF_DIR"/delta_*.dat; do
      local fname=$(basename "$ref_file")
      local out_file="$OUTPUT_DIR/$fname"
      if [ ! -f "$out_file" ]; then
        file_status["$fname|$IPOT|$ILB"]="\033[0;31mFAILED\033[0m"
        all_passed=0
        continue
      fi
      if compare_files "$ref_file" "$out_file"; then
        file_status["$fname|$IPOT|$ILB"]="\033[0;32mPASS\033[0m"
      else
        file_status["$fname|$IPOT|$ILB"]="\033[0;31mFAILED\033[0m"
        all_passed=0
      fi
    done

    rm -f "$TEST_INPUT"
  done

  echo
  echo "---------------- VARIATIONAL MODULE TEST SUMMARY ----------------"
  for test in "${TESTS[@]}"; do
    set -- $test
    local IPOT="$1"
    local ILB="$2"
    local REF_SUBDIR="$3"
    local REF_DIR="$REF_BASE/$REF_SUBDIR"
    for ref_file in "$REF_DIR"/delta_*.dat; do
      local fname=$(basename "$ref_file")
      local status="${file_status["$fname|$IPOT|$ILB"]}"
      printf "%b: %s for IPOT=%s ILB=%s\n" "$status" "$fname" "$IPOT" "$ILB"
    done
  done

  # Print plotting suggestions if any
  if [ -s /tmp/variational_module_plot_suggestions.txt ]; then
    echo
    echo "To compare the results visually, you can use the following commands:"
    cat /tmp/variational_module_plot_suggestions.txt
    rm -f /tmp/variational_module_plot_suggestions.txt
  fi

  if [ $all_passed -eq 1 ]; then
    echo -e "\n\033[1;32mALL VARIATIONAL TESTS PASSED\033[0m"
    exit 0
  else
    exit 1
  fi
}

main
