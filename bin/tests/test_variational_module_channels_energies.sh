#!/bin/bash


echo 
echo "============================================================"
echo "Testing variational module channels energies"
echo "============================================================"

AV18_2="output/test_AV18_2"
EFT_pless_15_2="output/test_EFT_pless_15_2"
EFT_pless_15_dynamic_2="output/test_EFT_pless_15_dynamic_2"
EXE="build/tests/test_NN_SCATTERING_VARIATIONAL_ENERGIES_CHANNELS.x"

need_make=0
for dir in "$AV18_2" "$EFT_pless_15_2" "$EFT_pless_15_dynamic_2"; do
  if [ ! -d "$dir" ] || [ "$(ls -1 "$dir"/delta_* 2>/dev/null | wc -l)" -ne 8 ]; then
    echo "Directory $dir missing or does not contain 8 files. Will run make -j and $EXE..."
    need_make=1
    break
  fi
done

if [ $need_make -eq 1 ]; then
  make -j
  $EXE
fi

for i in $AV18_2/delta_* ; do 
  diff-numerics $i bin/tests/test_files/scattering_phase_shifts/AV18/$(basename $i) -q
  if [ $? -eq 0 ]; then
    echo -e "\e[32m[PASSED]\e[0m No difference between $i and bin/tests/test_files/scattering_phase_shifts/AV18/$(basename $i)"
  else
    echo -e "\e[31m[FAILED]\e[0m Difference found between $i and bin/tests/test_files/scattering_phase_shifts/AV18/$(basename $i)"
  fi 
done

for i in $EFT_pless_15_2/delta_* ; do 
  diff-numerics $i bin/tests/test_files/scattering_phase_shifts/EFT_pless_15/$(basename $i) -q
  if [ $? -eq 0 ]; then
    echo -e "\e[32m[PASSED]\e[0m No difference between $i and bin/tests/test_files/scattering_phase_shifts/EFT_pless_15/$(basename $i)"
  else
    echo -e "\e[31m[FAILED]\e[0m Difference found between $i and bin/tests/test_files/scattering_phase_shifts/EFT_pless_15/$(basename $i)"
  fi 
done

for i in $EFT_pless_15_dynamic_2/delta_* ; do 
  diff-numerics $i bin/tests/test_files/scattering_phase_shifts/EFT_pless_15/$(basename $i) -q
  if [ $? -eq 0 ]; then
    echo -e "\e[32m[PASSED]\e[0m No difference between $i and bin/tests/test_files/scattering_phase_shifts/EFT_pless_15/$(basename $i)"
  else
    echo -e "\e[31m[FAILED]\e[0m Difference found between $i and bin/tests/test_files/scattering_phase_shifts/EFT_pless_15/$(basename $i)"
  fi 
done