# Fortran compiler
FC = gfortran

# Compiler and linker flags
# Set DEBUG to 1 for debug build, 0 for release (default)
DEBUG ?= 0
GMON ?= 0

ifeq ($(DEBUG),1)
  FFLAGS = -c -J$(BUILD_DIR) -I$(BUILD_DIR) \
           -Wall -g -fcheck=all -finit-real=snan -finit-local-zero \
           -fbacktrace -fdefault-real-8 -fdefault-double-8
  LDFLAGS = -Wall -g -fcheck=all -finit-real=snan -finit-local-zero \
            -fdefault-real-8 -fdefault-double-8
else
  FFLAGS = -c -J$(BUILD_DIR) -I$(BUILD_DIR) \
           -Wall -O3 -march=native -funroll-loops -ftree-vectorize \
           -fdefault-real-8 -fdefault-double-8
  LDFLAGS = -Wall -O3 -march=native -funroll-loops -ftree-vectorize \
            -fdefault-real-8 -fdefault-double-8
endif

OMP ?= 0
ifeq ($(OMP),1)	
	FFLAGS += -fopenmp
	LDFLAGS += -fopenmp
endif

ifeq ($(GMON),1)
	FFLAGS += -pg
	LDFLAGS += -pg
endif

LDFLAGS += -lgsl -lgslcblas -llapack -lblas


# Define directory paths
SRC_DIR := src
BUILD_DIR := build
LOG_DIR := $(BUILD_DIR)/logs
DEP_DIR := $(BUILD_DIR)/dep

# Find all .f90 files (recursively in src/)
ALL_SOURCES := $(shell find $(SRC_DIR) -name "*.f90")

# Select main program sources (e.g., src/main1.f90, src/main2.f90, ...)
MAIN_SOURCES := $(filter $(SRC_DIR)/main%.f90, $(ALL_SOURCES))

# Create lists of corresponding object files and dependency files
OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(ALL_SOURCES))
DEPFILES := $(patsubst $(SRC_DIR)/%.f90,$(DEP_DIR)/%.d,$(ALL_SOURCES))
EXECUTABLES := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.x,$(MAIN_SOURCES))

# Default rule
all: $(LOG_DIR) $(DEP_DIR) $(EXECUTABLES)

# Rule to build executables: one for each main program
$(BUILD_DIR)/%.x: $(BUILD_DIR)/%.o $(filter-out $<, $(OBJECTS))
	@echo "Linking $@"
	# $^ = all dependencies (main.o and module .o files)
	$(FC) $^ $(LDFLAGS) -o $@ > $(LOG_DIR)/$(notdir $@).link.log 2>&1 || (cat $(LOG_DIR)/$(notdir $@).link.log && false)
# $(BUILD_DIR)/%.x: $(BUILD_DIR)/%.o $(BUILD_DIR)/libs/math/coulomb_FG.o $(BUILD_DIR)/libs/math/algebra.o $(BUILD_DIR)/libs/utils.o
# 	@echo "Linking $@"
# 	$(FC) $^ $(LDFLAGS) -o $@ > $(LOG_DIR)/$(notdir $@).link.log 2>&1 || (cat $(LOG_DIR)/$(notdir $@).link.log && false)

# Rule to compile each .f90 file into a .o object
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	$(FC) $(FFLAGS) -o $@ $< > $(LOG_DIR)/$(subst /,_,$(notdir $<)).log 2>&1

# Rule to generate dependency files using the Python script
$(DEP_DIR)/%.d: $(SRC_DIR)/%.f90
	@echo "Generating dependencies for $<"
	@mkdir -p $(dir $@)
	@python3 tools/scan_deps.py $< || echo "$<: dependency scan failed" > $@

# Create required directories
$(LOG_DIR):
	mkdir -p $(LOG_DIR)

$(DEP_DIR):
	mkdir -p $(DEP_DIR)

# Run the script in bin directory
run:
	./bin/runner.sh

# Clean rule to delete all build artifacts
clean:
	rm -rvf $(BUILD_DIR)

# Include all dependency files if they exist
-include $(DEPFILES)

.PHONY: all clean run
