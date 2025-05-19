#!/bin/bash

make

# Default values for the parameters
EMAX=1.D0
NE=200
J=1
L=0
S=1
TZ=0
IPOT=18
ILB=1
LEMP=0
PRINT_COEFF=".FALSE."

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -emax|--Emax)
      EMAX="$2"
      shift 2
      ;;
    -ne|--NE)
      NE="$2"
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
    -print|--PRINT)
      PRINT_COEFF=".TRUE."
      shift
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
  EMAX = $EMAX,
  NE = $NE,
  J = $J,
  L = $L,
  S = $S,
  TZ = $TZ,
  IPOT = $IPOT,
  ILB = $ILB,
  LEMP = $LEMP,
  PRINT_COEFFICIENTS = $PRINT_COEFF,
/
EOF

# Run the Fortran program with the namelist file
./build/main_scattering_NN_variazional_method.x "$NAMELIST_FILE"

# Clean up the temporary file
rm -f "$NAMELIST_FILE"