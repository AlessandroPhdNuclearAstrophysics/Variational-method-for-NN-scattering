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

OMP ?= 1
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
TEST_DIR := $(SRC_DIR)/tests
LOG_DIR := $(BUILD_DIR)/logs
DEP_DIR := $(BUILD_DIR)/dep
OUT_DIR := output

# Find all .f90 files (recursively in src/)
ALL_SOURCES := $(shell find $(SRC_DIR) -name "*.f90")

# Select main program sources: main_*.f90 in src/ and all .f90 in $(TEST_DIR)
MAIN_SOURCES := $(filter $(SRC_DIR)/main_%.f90, $(ALL_SOURCES))
TEST_SOURCES := $(shell find $(TEST_DIR) -name "*.f90")
ALL_MAIN_SOURCES := $(MAIN_SOURCES) $(TEST_SOURCES)

# Create lists of corresponding object files and dependency files
OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(ALL_SOURCES))
DEPFILES := $(patsubst $(SRC_DIR)/%.f90,$(DEP_DIR)/%.d,$(ALL_SOURCES))
TEST_OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(TEST_SOURCES))
TEST_DEPFILES := $(patsubst $(SRC_DIR)/%.f90,$(DEP_DIR)/%.d,$(TEST_SOURCES))

# Executables for main_*.f90 and test programs
EXECUTABLES := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.x,$(MAIN_SOURCES))
TEST_EXECUTABLES := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.x,$(TEST_SOURCES))
ALL_EXECUTABLES := $(EXECUTABLES) $(TEST_EXECUTABLES)

# List of all main and test object files
MAIN_OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(MAIN_SOURCES))
ALL_MAIN_OBJECTS := $(MAIN_OBJECTS) $(TEST_OBJECTS)

# List of all object files except main and test objects
NON_MAIN_OBJECTS := $(filter-out $(ALL_MAIN_OBJECTS),$(OBJECTS))

# Default rule
all: $(LOG_DIR) $(DEP_DIR) $(ALL_EXECUTABLES)

# Rule to build executables: link only its own object and NON_MAIN_OBJECTS
$(BUILD_DIR)/%.x: $(BUILD_DIR)/%.o $(NON_MAIN_OBJECTS)
	@echo "Linking $@"
	$(FC) $^ $(LDFLAGS) -o $@ > $(LOG_DIR)/$(notdir $@).link.log 2>&1 || (cat $(LOG_DIR)/$(notdir $@).link.log && false)

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

# Check the graphs
POT ?= AV18
check_graphs:
	@echo "Checking graphs..."
	xmgrace $(OUT_DIR)/$(POT)/delta_1S0.dat ../07b-Potential_revised/scattering/output/$(POT)_Laura/np/0p.dat; \
	xmgrace $(OUT_DIR)/$(POT)/delta_3P0.dat ../07b-Potential_revised/scattering/output/$(POT)_Laura/np/0m.dat; \
	xmgrace -nxy $(OUT_DIR)/$(POT)/delta_3S1-3D1.dat -nxy ../07b-Potential_revised/scattering/output/$(POT)_Laura/np/1p.dat; \
	xmgrace $(OUT_DIR)/$(POT)/delta_1P1.dat ../07b-Potential_revised/scattering/output/$(POT)_Laura/np/1m.dat; \
	xmgrace $(OUT_DIR)/$(POT)/delta_3P1.dat -block ../07b-Potential_revised/scattering/output/$(POT)_Laura/np/1m.dat -bxy 1:3; \
	xmgrace -nxy $(OUT_DIR)/$(POT)/delta_3P2-3F2.dat -nxy ../07b-Potential_revised/scattering/output/$(POT)_Laura/np/2m.dat; \
	xmgrace $(OUT_DIR)/$(POT)/delta_1D2.dat ../07b-Potential_revised/scattering/output/$(POT)_Laura/np/2p.dat; \
	xmgrace $(OUT_DIR)/$(POT)/delta_3D2.dat -block ../07b-Potential_revised/scattering/output/$(POT)_Laura/np/2p.dat -bxy 1:3; \

# Clean rule to delete all build artifacts
clean:
	rm -rvf $(BUILD_DIR)

# Delete output files
delete_out:
	@rm -rvf $(shell find $(OUT_DIR) -type f)

# Doxygen documentation target (Doxyfile generated dynamically)
doc:
	@echo "Generating Doxyfile..."
	@echo "PROJECT_NAME           = Variational scattering" > Doxyfile
	@echo "OUTPUT_DIRECTORY       = doc" >> Doxyfile
	@echo "GENERATE_HTML          = YES" >> Doxyfile
	@echo "GENERATE_LATEX         = NO" >> Doxyfile
	@echo "RECURSIVE              = YES" >> Doxyfile
	@echo "INPUT                  = src" >> Doxyfile
	@echo "FILE_PATTERNS          = *.f90" >> Doxyfile
	@echo "EXTRACT_ALL            = YES" >> Doxyfile
	@echo "OPTIMIZE_FOR_FORTRAN   = YES" >> Doxyfile
	@echo "QUIET                  = YES" >> Doxyfile
	@echo "JAVADOC_AUTOBRIEF      = YES" >> Doxyfile
	@echo "MULTILINE_CPP_IS_BRIEF = YES" >> Doxyfile
	@echo "DOT_IMAGE_FORMAT       = svg" >> Doxyfile
	@echo "INTERACTIVE_SVG        = YES" >> Doxyfile
	@echo "HAVE_DOT           		= YES" >> Doxyfile
	@echo "CALL_GRAPH         		= YES" >> Doxyfile
	@echo "CALLER_GRAPH       		= YES" >> Doxyfile
	@echo "DOT_GRAPH_MAX_NODES		= 100" >> Doxyfile
	@echo "DOT_CLEANUP        		= YES" >> Doxyfile
	@echo "DOT_GROUP_GRAPHS   		= YES" >> Doxyfile
	@doxygen Doxyfile
	@rm -vf Doxyfile
	@echo "Doxygen documentation generated in doc/"
	@echo "Open doc/html/index.html in a web browser to view the documentation."
	@echo "Doxygen documentation generation complete."


delete_doc:
	@rm -rvf doc

# Test target: build and run all test executables in $(TEST_DIR)
test: $(TEST_EXECUTABLES)
	@for t in $(TEST_EXECUTABLES); do \
		echo "Running $$t"; \
		$$t || exit 1; \
	done

# Include all dependency files if they exist
-include $(DEPFILES) $(TEST_DEPFILES)

.PHONY: all clean run delete_out check_graphs doc delete_doc test