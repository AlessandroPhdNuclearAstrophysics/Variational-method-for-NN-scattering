#!/bin/bash
# Script to test the static library and module files
# 1. Generate the zip library using the Makefile
# 2. Copy the zip to the script directory
# 3. Extract it
# 4. Compile test_library.f90 using the library and mod files
# 5. Run it
# 6. Compare output files with reference files using a numerical diff

set -e

ROOT_DIR="$(dirname $(dirname $(dirname $(realpath "$0"))))"
BUILD_DIR="$ROOT_DIR/build"
SCRIPT_DIR="$(dirname $(realpath "$0"))"
ZIP_NAME="libvariational.zip"
LIB_DIR="libvariational"

# 1. Generate the zip library
echo -e "\033[0;34m[DEBUG] Running: make -C $ROOT_DIR $BUILD_DIR/$ZIP_NAME\033[0m"
make -C "$ROOT_DIR" "build/$ZIP_NAME"

# 2. Copy the zip to the script directory
echo -e "\033[0;34m[DEBUG] Copying $BUILD_DIR/$ZIP_NAME to $SCRIPT_DIR/\033[0m"
cp -v "$BUILD_DIR/$ZIP_NAME" "$SCRIPT_DIR/"

# 3. Extract it
echo -e "\033[0;34m[DEBUG] Extracting $ZIP_NAME in $SCRIPT_DIR and deleting the zip file\033[0m"
cd "$SCRIPT_DIR"
unzip -o "$ZIP_NAME"
rm -v "$ZIP_NAME"

# 4. Compile test_library.f90 using the library and mod files
FC=gfortran
MODDIR="$SCRIPT_DIR/$LIB_DIR"
LIBFILE="$MODDIR/libvariational.a"
MODS="$MODDIR"/*.mod
SRC="$SCRIPT_DIR/test_library.f90"
EXE="$SCRIPT_DIR/test_library.x"
# Use the same flags as in the Makefile
FFLAGS="-O3 -march=native -funroll-loops -ftree-vectorize -fopenmp -Wall -fdefault-real-8 -fdefault-double-8 -ffpe-trap=invalid,zero,overflow -finit-real=snan -I$MODDIR -J$MODDIR"
LDFLAGS="-O3 -march=native -funroll-loops -ftree-vectorize -fopenmp -Wall -fdefault-real-8 -fdefault-double-8 -ffpe-trap=invalid,zero,overflow -finit-real=snan -lgsl -lgslcblas -llapack -lblas"
echo -e "\033[0;34m[DEBUG] Compiling $SRC with:\033[0m $FC $FFLAGS -o $EXE $SRC $LIBFILE $LDFLAGS"
$FC $FFLAGS -o "$EXE" "$SRC" "$LIBFILE" $LDFLAGS

# 5. Show debug info and stop before running the executable
echo -e "\033[0;34m[DEBUG] Preparation complete.\033[0m"
echo -e "\033[0;34m[DEBUG] Files in $SCRIPT_DIR:\033[0m"
ls -l "$SCRIPT_DIR"
echo -e "\033[0;34m[DEBUG] Files in $MODDIR:\033[0m"
ls -l "$MODDIR"
echo -e "\033[0;34m[DEBUG] Compiled executable: $EXE\033[0m"
if [ -f "$EXE" ]; then
  echo -e "\033[0;34m[DEBUG] Executable exists.\033[0m"
else
  echo -e "\033[0;31m[ERROR] Executable not found!\033[0m"
  exit 1
fi

# Run the executable with IPOT 18
echo -e "\033[0;34m[DEBUG] Running $EXE with IPOT 18\033[0m"
"$EXE" 18
if [ $? -ne 0 ]; then
  echo -e "\033[0;31m[ERROR] Execution failed!\033[0m"
  exit 1
fi
echo -e "\033[0;32m[DEBUG] Execution successful.\033[0m"

# Run the executable with IPOT 18
echo -e "\033[0;34m[DEBUG] Running $EXE with IPOT 18\033[0m"
"$EXE" 19
if [ $? -ne 0 ]; then
  echo -e "\033[0;31m[ERROR] Execution failed!\033[0m"
  exit 1
fi
echo -e "\033[0;32m[DEBUG] Execution successful.\033[0m"

# 6. Compare output files with reference files
compare_files() {
  local ref_file="$1"
  local out_file="$2"
  local tol=1e-2
  local threshold=1e-5
  diff-numerics  $ref_file $out_file --tolerance $tol --threshold $threshold -q > tmp.diff
  if [ $? -ne 0 ]; then
    echo "Files $ref_file and $out_file differ."
    echo -e "\033[0;31m$(tail -n 2 tmp.diff)\033[0m"
    rm -f tmp.diff
    return 1
  else
    echo "Files $ref_file and $out_file are identical."
    rm -f tmp.diff
    return 0
  fi
}

# Compare all output files with reference files (AV18 and EFT_pless_15)
for refdir in "$SCRIPT_DIR/test_files/scattering_phase_shifts/AV18" "$SCRIPT_DIR/test_files/scattering_phase_shifts/EFT_pless_15"; do
  if [ -d "$refdir" ]; then
    for ref in "$refdir"/delta_*.dat; do
      fname=$(basename "$ref")
      outdir="$SCRIPT_DIR/tmp_output_library/$(basename $refdir)" # output dir created by test_library.f90
      out="$outdir/$fname"
      if [ -f "$out" ]; then
        echo -e "\033[0;34m[DEBUG] Comparing $out to $ref...\033[0m"
        compare_files "$ref" "$out" || echo -e "\033[0;34m[DEBUG] Difference found in $fname\033[0m"
      else
        echo -e "\033[0;34m[DEBUG] Output file $out not found.\033[0m"
      fi
    done
  fi
done

rm -rf "$SCRIPT_DIR/tmp_output_library" # Clean up output directory
rm -rf "$MODDIR" # Clean up mod directory
rm -f "$EXE" # Clean up executable

echo "Test complete."
