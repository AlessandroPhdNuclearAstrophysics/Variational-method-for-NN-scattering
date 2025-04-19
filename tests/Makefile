.PHONY: all copy_structure clean

GF := gfortran
F90FLAGS := -g -Wall -fdefault-real-8 -fdefault-double-8 -frange-check \
		-ffree-form -fcheck=all -fbacktrace -fcheck=bounds \
		-fcheck=recursion 

SRC_TESTS_DIR := ../src/tests
BUILD_DIR := ../build


TEST_DIRS := $(shell find $(SRC_TESTS_DIR) -mindepth 1 -type d | sed 's|$(SRC_TESTS_DIR)/||')
SRC_F90 := $(shell find $(SRC_TESTS_DIR) -name '*.f90')
TAR_F90 := $(patsubst $(SRC_TESTS_DIR)/%.f90, %.x, $(SRC_F90))

OBJ_F90 := $(BUILD_DIR)/libs/potentials/av18.o \
					 $(BUILD_DIR)/libs/potentials/av18_f.o 

all: copy_structure make_build $(TAR_F90)


make_build:
	@$(MAKE) -j -C $(BUILD_DIR) all

# Target to copy the subfolder structure from src/tests to tests
copy_structure:
	@$(foreach dir, $(TEST_DIRS), $(shell mkdir -p ./$(dir)))

$(TAR_F90): %.x: $(SRC_TESTS_DIR)/%.f90 $(OBJ_F90) 
	$(GF) $(F90FLAGS) $^ -o $@

clean:
	@echo "Cleaning up the tests directory..."
	@find $(TESTS_DIR) -mindepth 1 -type d -exec rm -rf {} +
	@echo "Tests directory cleaned."