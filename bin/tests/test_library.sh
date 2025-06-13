#!/bin/bash
# Script to test the static library and module files
# 1. Generate the zip library using the Makefile
# 2. Copy the zip to the script directory
# 3. Extract it
# 4. Compile test_library.f90 using the library and mod files
# 5. Run it
# 6. Compare output files with reference files using a numerical diff

echo "========================================"
echo "LIBRARY TEST"
echo "========================================"

DEBUG=0
if [[ "$1" == "--debug" ]]; then
  DEBUG=1
fi

# set -e    # Uncomment to stop on errors

ROOT_DIR="$(dirname $(dirname $(dirname $(realpath "$0"))))"
if [[ DEBUG -eq 1 ]]; then
  BUILD_DIR="$ROOT_DIR/build/debug"
  ZIP_NAME="libvariational_debug.zip"
else
  BUILD_DIR="$ROOT_DIR/build/release"
  ZIP_NAME="libvariational_release.zip"
fi
SCRIPT_DIR="$(dirname $(realpath "$0"))"
LIB_DIR="libvariational"

# 1. Generate the zip library
if [ $DEBUG -eq 1 ]; then
  echo -e "\033[0;34m[DEBUG] Running: make -C $ROOT_DIR $BUILD_DIR/$ZIP_NAME\033[0m"
  make -C "$ROOT_DIR" export DEBUG=1
else
  make -C "$ROOT_DIR" export 
fi

# 2. Copy the zip to the script directory
if [ $DEBUG -eq 1 ]; then
  echo -e "\033[0;34m[DEBUG] Copying $BUILD_DIR/$ZIP_NAME to $SCRIPT_DIR/\033[0m"
  cp -v "$BUILD_DIR/$ZIP_NAME" "$SCRIPT_DIR/"
else
  cp "$BUILD_DIR/$ZIP_NAME" "$SCRIPT_DIR/"
fi

# 3. Extract it
if [ $DEBUG -eq 1 ]; then
  echo -e "\033[0;34m[DEBUG] Extracting $ZIP_NAME in $SCRIPT_DIR and deleting the zip file\033[0m"
fi
cd "$SCRIPT_DIR"
unzip -o "$ZIP_NAME" >/dev/null
rm -f "$ZIP_NAME"

# 4. Compile test_library.f90 using the library and mod files
FC=gfortran
MODDIR="$SCRIPT_DIR/$LIB_DIR"
LIBFILE="$MODDIR/libvariational.a"
MODS="$MODDIR"/*.mod
SRC="$SCRIPT_DIR/test_library.f90"
EXE="$SCRIPT_DIR/test_library.x"
FFLAGS="-O3 -march=native -funroll-loops -ftree-vectorize -fopenmp -Wall -fdefault-real-8 -fdefault-double-8 -ffpe-trap=invalid,zero,overflow -finit-real=snan -I$MODDIR -J$MODDIR"
LDFLAGS="-lgsl -lgslcblas -llapack -lblas"
if [ $DEBUG -eq 1 ]; then
  echo -e "\033[0;34m[DEBUG] Compiling $SRC with:"
  echo -e "\033[0m $FC $FFLAGS -o $EXE $SRC $LIBFILE $LDFLAGS"
  $FC $FFLAGS -o "$EXE" "$SRC" "$LIBFILE" $LDFLAGS
else
  $FC $FFLAGS -o "$EXE" "$SRC" "$LIBFILE" $LDFLAGS >/dev/null 2>&1
fi

# 5. Show debug info and stop before running the executable
if [ $DEBUG -eq 1 ]; then
  echo -e "\033[0;34m[DEBUG] Preparation complete.\033[0m"
  echo -e "\033[0;34m[DEBUG] Files in $SCRIPT_DIR:\033[0m"
  ls -l "$SCRIPT_DIR"
  echo -e "\033[0;34m[DEBUG] Files in $MODDIR:\033[0m"
  ls -l "$MODDIR"
  echo -e "\033[0;34m[DEBUG] Compiled executable: $EXE\033[0m"
fi
if [ -f "$EXE" ]; then
  if [ $DEBUG -eq 1 ]; then
    echo -e "\033[0;34m[DEBUG] Executable exists.\033[0m"
  fi
else
  echo -e "\033[0;31m[ERROR] Executable not found!\033[0m"
  exit 1
fi

# Animation function
run_with_dots() {
  CMD=("$@")
  if [ $DEBUG -eq 1 ]; then
    "${CMD[@]}" &
    pid=$!
    dots=1
    while kill -0 $pid 2>/dev/null; do
      case $dots in
        1) printf "\r." ;;
        2) printf "\r.." ;;
        3) printf "\r..." ;;
      esac
      dots=$((dots % 3 + 1))
      sleep 0.4
    done
    wait $pid
    status=$?
    printf "\r   \r" # clear line
    return $status
  else
    # Show dots animation in ALL modes, but suppress program output unless in debug
    "$@" >.tmp_run_with_dots_out 2>&1 &
    pid=$!
    dots=1
    while kill -0 $pid 2>/dev/null; do
      case $dots in
        1) printf "\r." ;;
        2) printf "\r.." ;;
        3) printf "\r..." ;;
      esac
      dots=$((dots % 3 + 1))
      sleep 0.4
    done
    wait $pid
    status=$?
    printf "\r   \r" # clear line
    if [ $DEBUG -eq 1 ]; then
      cat .tmp_run_with_dots_out
    fi
    rm -f .tmp_run_with_dots_out
    return $status
  fi
}

# Run the executable with IPOT 18
IPOT=18
echo "Preparing outputs with ipot = $IPOT"
run_with_dots "$EXE" $IPOT
if [ $? -ne 0 ]; then
  echo -e "\033[0;31m[ERROR] Execution failed!\033[0m"
  exit 1
fi
echo "Preparation for potential $IPOT done."
if [ $DEBUG -eq 1 ]; then
  echo -e "\033[0;32m[DEBUG] Execution successful.\033[0m"
fi

# Run the executable with IPOT 19
IPOT=19
echo "Preparing outputs with ipot = $IPOT"
run_with_dots "$EXE" $IPOT
if [ $? -ne 0 ]; then
  echo -e "\033[0;31m[ERROR] Execution failed!\033[0m"
  exit 1
fi
echo "Preparation for potential $IPOT done."
if [ $DEBUG -eq 1 ]; then
  echo -e "\033[0;32m[DEBUG] Execution successful.\033[0m"
fi

# 6. Compare output files with reference files
compare_files() {
  local ref_file="$1"
  local out_file="$2"
  local tol=1e-2
  local threshold=1e-5

  if command -v tput >/dev/null 2>&1; then
    RED="$(tput setaf 1)"
    GREEN="$(tput setaf 2)"
    RESET="$(tput sgr0)"
  else
    RED="\033[0;31m"
    GREEN="\033[0;32m"
    RESET="\033[0m"
  fi

  diff-numerics  $ref_file $out_file --tolerance $tol --threshold $threshold -q > tmp.diff
  if [ $? -ne 0 ]; then
    rel_ref_file="${ref_file#*$ROOT_DIR/}"
    rel_out_file="${out_file#*$ROOT_DIR/}"
    echo -e "${RED}Files $rel_ref_file and $rel_out_file differ.${RESET}"
    rm -f tmp.diff
    return 1
  else
    rel_ref_file="${ref_file#*$ROOT_DIR/}"
    rel_out_file="${out_file#*$ROOT_DIR/}"
    echo -e "${GREEN}Files $rel_ref_file and $rel_out_file are identical.${RESET}"
    rm -f tmp.diff
    return 0
  fi
}

# Compare all output files with reference files (AV18 and EFT_pless_15)
for refdir in "$SCRIPT_DIR/test_files/scattering_phase_shifts/AV18" "$SCRIPT_DIR/test_files/scattering_phase_shifts/EFT_pless_15"; do
  if [ -d "$refdir" ]; then
    for ref in "$refdir"/delta_*.dat; do
      fname=$(basename "$ref")
      # Skip files containing 'Stapp' in the name
      if [[ "$fname" == *BB* ]]; then
        continue
      fi
      echo "Comparing $fname in $refdir"
      outdir="$SCRIPT_DIR/tmp_output_library/$(basename $refdir)"
      out="$outdir/$fname"
      if [ -f "$out" ]; then
        if [ $DEBUG -eq 1 ]; then
          echo -e "\033[0;34m[DEBUG] Comparing $out to $ref...\033[0m"
        fi
        compare_files "$ref" "$out"
      else
        if [ $DEBUG -eq 1 ]; then
          echo -e "\033[0;34m[DEBUG] Output file $out not found.\033[0m"
        fi
      fi
    done
  fi
done

rm -rf "$SCRIPT_DIR/tmp_output_library"
rm -rf "$MODDIR"
rm -f "$EXE"

echo "Test complete."
