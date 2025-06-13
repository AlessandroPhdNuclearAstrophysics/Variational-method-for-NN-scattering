#!/bin/bash

if [ "$1" == "--debug" ]; then
  BUILD_DIR="build/debug"
else
  BUILD_DIR="build/release"
fi
# Check if the build directory exists
if [ ! -d "$BUILD_DIR" ]; then
  echo "Build directory $BUILD_DIR does not exist. Please run 'make' first."
  exit 1
fi

EXE="$BUILD_DIR/main_scattering_NN_variazional_method.x"
# Check if the executable exists
if [ ! -f "$EXE" ]; then
  echo "Executable $EXE does not exist. Please run 'make' first."
  exit 1
fi

# Default values for the parameters
EMAX=1.D0
NE=200
TZ=0
IPOT=18
ILB=1
LEMP=0
PRINT_COEFF=".FALSE."
OUTPUT_DIR="output/POT$IPOT/"

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
    -out_dir|--OUT_DIR)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

echo "OUTPUT_DIR: $OUTPUT_DIR"
# Create a temporary namelist input file
NAMELIST_FILE=$(mktemp)
cat > "$NAMELIST_FILE" << EOF
&IN
  EMAX = $EMAX,
  NE = $NE,
  TZ = $TZ,
  IPOT = $IPOT,
  ILB = $ILB,
  LEMP = $LEMP,
  PRINT_COEFFICIENTS = $PRINT_COEFF,
  OUT_DIR = "$OUTPUT_DIR",
/
EOF

# Create the output directory only if it does not exist
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

# Run the Fortran program with the namelist file
"$EXE" "$NAMELIST_FILE"

# Clean up the temporary file
rm -f "$NAMELIST_FILE"