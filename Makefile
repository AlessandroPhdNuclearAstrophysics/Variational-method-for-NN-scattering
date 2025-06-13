.DEFAULT_GOAL := all

# Fortran compiler
FC = gfortran

# Compiler and linker flags
# Set DEBUG to 1 for debug build, 0 for release (default)
DEBUG ?= 0
GMON ?= 0

# Define directory paths
SRC_DIR := src
BUILD_DIR ?= build/$(if $(filter 1,$(DEBUG)),debug,release)
TEST_DIR := $(SRC_DIR)/tests
LOG_DIR := $(BUILD_DIR)/logs
DEP_DIR := $(BUILD_DIR)/dep
OUT_DIR := output

# Common flags
COMMON_FLAGS = -fimplicit-none -ffree-line-length-none -Wall -Wextra
COMMON_FLAGS+= -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-unused-label -Wno-unused-dummy-argument
COMMON_FLAGS+= -Wno-implicit-interface -Wno-argument-mismatch
COMMON_FLAGS+= -Wno-compare-reals -Wno-maybe-uninitialized
COMMON_FLAGS+= -fdefault-real-8 -fdefault-double-8 -ffpe-trap=invalid,zero,overflow -finit-real=snan

# Set the flags for the Fortran compiler
FFLAGS = -c -J$(BUILD_DIR) -I$(BUILD_DIR) $(COMMON_FLAGS)
LDFLAGS = $(COMMON_FLAGS)

ifeq ($(DEBUG),1)
	FFLAGS  += -g -fcheck=all -finit-local-zero -fbacktrace
	LDFLAGS += -g -fcheck=all
else
	FFLAGS  += -O3 -march=native -funroll-loops -ftree-vectorize
	LDFLAGS += -O3 -march=native -funroll-loops -ftree-vectorize
endif

# Enable OpenMP
FFLAGS += -fopenmp
LDFLAGS += -fopenmp

ifeq ($(GMON),1)
	FFLAGS += -pg
	LDFLAGS += -pg
endif

LDFLAGS += -lgsl -lgslcblas -llapack -lblas

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
all: $(LOG_DIR) $(DEP_DIR) $(OUT_DIR) $(ALL_EXECUTABLES)

# Rule to build executables: link only its own object and NON_MAIN_OBJECTS
$(BUILD_DIR)/%.x: $(BUILD_DIR)/%.o $(NON_MAIN_OBJECTS)
	@echo "\033[0;32mLinking $@\033[0m"
	$(FC) $^ $(LDFLAGS) -o $@ > $(LOG_DIR)/$(notdir $@).link.log 2>&1 || (cat $(LOG_DIR)/$(notdir $@).link.log && false)

# Rule to compile each .f90 file into a .o object
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(dir $@)
	@echo "\033[0;32mCompiling $<\033[0m"
	$(FC) $(FFLAGS) -o $@ $< > $(LOG_DIR)/$(subst /,_,$(notdir $<)).log 2>&1 || (cat $(LOG_DIR)/$(subst /,_,$(notdir $<)).log && false)

# Rule to generate dependency files using the Python script
$(DEP_DIR)/%.d: $(SRC_DIR)/%.f90
	@echo "Generating dependencies for $<"
	@mkdir -p $(dir $@)
	@python3 tools/scan_deps.py $<  $(BUILD_DIR) || echo "$<: dependency scan failed" > $@

# Create required directories
$(LOG_DIR):
	mkdir -p $(LOG_DIR)

$(DEP_DIR):
	mkdir -p $(DEP_DIR)

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

# Generate dependency graphs using Graphviz
generate_dependency_graphs: $(DEP_FILES)
	@echo "Generating dependency graphs..."
	@mkdir -p dependency_graphs
	@python3.12 tools/dependency_tree.py $(DEP_DIR) -d dependency_graphs --auto-adjust -l circular -f -g -p 3 -s 50 --remove-transitive

SRC_XMGRACE := $(shell find . -name "*.agr")
OUT_XMGRACE := $(patsubst %.agr, %.pdf, $(SRC_XMGRACE))

graphs: $(OUT_XMGRACE)

$(OUT_XMGRACE): %.pdf: %.agr
	@echo "Converting $< to PDF..."
	xmgrace -hardcopy $< -hdevice EPS -printfile $(subst .agr,.eps,$<)
	epstopdf $(subst .agr,.eps,$<) -o $@
	@rm -vf $(subst .agr,.eps,$<) 

# Clean rule to delete all build artifacts
clean:
	rm -rvf $(BUILD_DIR) $(shell find $(OUT_DIR) -mindepth 2 -name "*.eps") dependency_graphs/* libvariational.a

# Delete output files
delete_out:
	@rm -rvf $(shell find $(OUT_DIR) -type d)

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

GPROF_DIR := $(LOG_DIR)/gprof
analyze_gprof_dynamic_LECS:
	@echo "Analyzing gprof data for build/tests/test_variational_module_dynamic_LECS.x..."
	@make clean
	@echo "Rebuilding executable with profiling enabled..."
	@$(MAKE) GMON=1 DEBUG=$(DEBUG) -j
	@echo
	@echo "\033[0;32mRunning gprof analysis on build/tests/test_variational_module_dynamic_LECS.x...\033[0m"
	@mkdir -p $(GPROF_DIR)
	@if [ -f build/tests/test_variational_module_dynamic_LECS.x ]; then \
		echo "Analyzing build/tests/test_variational_module_dynamic_LECS.x"; \
		time ./build/tests/test_variational_module_dynamic_LECS.x > /dev/null ; \
		gprof build/tests/test_variational_module_dynamic_LECS.x gmon.out > $(GPROF_DIR)/test_variational_module_dynamic_LECS.gprof.log; \
		rm -vf gmon.out; \
	else \
		echo "Executable build/tests/test_variational_module_dynamic_LECS.x not found, skipping."; \
	fi
	@echo "Gprof analysis complete. Output saved to $(GPROF_DIR)/test_variational_module_dynamic_LECS.gprof.log"


delete_doc:
	@rm -rvf doc

# Test target: build and run all test executables in $(TEST_DIR)
test: $(TEST_EXECUTABLES)
	@for t in $(TEST_EXECUTABLES); do \
		echo ; \
		echo "\033[1;34mRunning $$t\033[0m"; \
		$$t || (echo "\033[1;31mERROR:Something went wrong\033[0m"); \
	done
	@echo "\033[1;34mRunning test_variational_module.x \033[0m"
	@./bin/tests/test_variational_module.sh || (echo "\033[1;31mERROR:Something went wrong\033[0m")
	@echo "\033[1;34mRunning test_library.x \033[0m"
	@./bin/tests/test_library.sh || (echo "\033[1;31mERROR:Something went wrong\033[0m")
	@echo "\033[1;34mRunning test_variational_module_channels_energies.x \033[0m"
	@./bin/tests/test_variational_module_channels_energies.sh || (echo "\033[1;31mERROR:Something went wrong\033[0m")


# Make static library from all the modules
LIBRARY_FILE := libvariational.a
$(BUILD_DIR)/$(LIBRARY_FILE): $(NON_MAIN_OBJECTS)
	@echo "\033[0;32mCreating static library $(LIBRARY_FILE)\033[0m"
	ar rcs $(BUILD_DIR)/$(LIBRARY_FILE) $^

# Export library with mode-specific zip name
export: $(BUILD_DIR)/$(LIBRARY_FILE)
	@echo "\033[0;32mExporting static library and module files...\033[0m"
	@mkdir -p $(BUILD_DIR)/libvariational/public_interfaces
	@for i in $$(find $(SRC_DIR)/libs -name "*.f90"); do \
		echo "Generating module interface for $$i..."; \
		fortran-module-interface-generator $$i -o $(BUILD_DIR)/libvariational/public_interfaces/$$(basename $$i) || true; \
	done
	@cp $(BUILD_DIR)/$(LIBRARY_FILE) $(BUILD_DIR)/libvariational/
	@cp $(BUILD_DIR)/*.mod $(BUILD_DIR)/libvariational/ 2>/dev/null || true
	@find . -name "lecs_eft.dat" -exec cp {} $(BUILD_DIR)/libvariational/ \; 2>/dev/null || true
	@MODE=$$( [ "$(DEBUG)" = "1" ] && echo "debug" || echo "release" ); \
	cd $(BUILD_DIR) && zip -r libvariational_$$MODE.zip libvariational
	@rm -rvf $(BUILD_DIR)/libvariational
	@echo "\033[0;32mLibrary exported as $(BUILD_DIR)/libvariational_`[ "$(DEBUG)" = "1" ] && echo debug || echo release`.zip\033[0m"


# Include all dependency files if they exist
-include $(DEPFILES) $(TEST_DEPFILES)

.PHONY: all clean delete_out check_graphs doc delete_doc test analyze_gprof_dynamic_LECS