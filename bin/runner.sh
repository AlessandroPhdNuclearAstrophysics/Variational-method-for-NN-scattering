#!/bin/bash

make

# Default values for the parameters
E=1.D0
J=1
L=0
S=1
TZ=0
IPOT=18
ILB=1
LEMP=1
VCOUL=false

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -e|--E)
      E="$2"
      shift 2
      ;;
    -j|--J)
      J="$2"
      shift 2
      ;;
    -l|--L)
      L="$2"
      shift 2
      ;;
    -s|--S)
      S="$2"
      shift 2
      ;;
    -tz|--TZ)
      TZ="$2"
      shift 2
      ;;
    -ipot|--IPOT)
      IPOT="$2"
      shift 2
      ;;
    -ilb|--ILB)
      ILB="$2"
      shift 2
      ;;
    -lemp|--LEMP)
      LEMP="$2"
      shift 2
      ;;
    -vcoul|--VCOUL)
      VCOUL="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Create a temporary namelist input file
NAMELIST_FILE=$(mktemp)
cat > "$NAMELIST_FILE" << EOF
&IN
  E = $E,
  J = $J,
  L = $L,
  S = $S,
  TZ = $TZ,
  IPOT = $IPOT,
  ILB = $ILB,
  LEMP = $LEMP,
  VCOUL = $VCOUL
/
EOF

# Run the Fortran program with the namelist file
./build/main_scattering_NN_variazional_method.x "$NAMELIST_FILE"

# Clean up the temporary file
rm -f "$NAMELIST_FILE"