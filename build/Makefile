.PHONY: all clean print_dirs

RM := rm -vf

BIN_DIR := ../bin
SRC_DIR := ../src
OUT_DIR := ../output

DIRS := $(shell find $(SRC_DIR) -mindepth 1 -type d | sed 's|$(SRC_DIR)/||')
SRC_F90 := $(shell find $(SRC_DIR) -path $(SRC_DIR)/tests -prune -o -name '*.f90' -print)
SRC_F   := $(shell find $(SRC_DIR) -path $(SRC_DIR)/tests -prune -o -name '*.f' -print)
OBJ_F90 := $(patsubst $(SRC_DIR)/%.f90, %.o, $(SRC_F90))
OBJ_F   := $(patsubst $(SRC_DIR)/%.f, %_f.o, $(SRC_F))
TARGET  := ../bin/scattering_NN_variazional_method.x

CF := gfortran
F90FLAGS := -g -Wall -fdefault-real-8 -fdefault-double-8 -frange-check \
		-ffree-form -fcheck=all -fbacktrace -fcheck=bounds \
		-fcheck=recursion -fcheck=pointer
FFLAGS := -g -Wall -fdefault-real-8 -fdefault-double-8 -frange-check \
		-fcheck=all -fbacktrace -fcheck=bounds \
		-fcheck=recursion -fcheck=pointer

all: $(TARGET)

prepare:
	$(foreach dir, $(DIRS), $(shell mkdir -p $(dir)))

$(OBJ_F90): %.o: $(SRC_DIR)/%.f90
	$(CF) $(F90FLAGS) -c $< -o $@

$(OBJ_F): %_f.o: $(SRC_DIR)/%.f
	$(CF) $(FFLAGS) -c $< -o $@

$(TARGET): ../bin/%.x: %.o $(OBJ_F90)
	$(CF) $(F90FLAGS) $^ -o $@
	@echo "Build complete: $@"

print_dirs:
	@echo "Dirs for the object files: $(DIRS)"

print:
	@echo $(EXES)

clean:
	@$(RM) $(OBJ_F90) $(OBJ_F) $(TARGET)
