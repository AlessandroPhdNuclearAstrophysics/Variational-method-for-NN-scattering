.PHONY: prepare all clean

# This Makefile is used to build the project using a separate build directory.
# It assumes that the build directory already contains a Makefile.
RM := rm -vf
MKDIR := mkdir -p
MAKE := make


BUILD_DIR := build	# Must exist before running this Makefile and contain a Makefile
BIN_DIR := bin
SRC_DIR := src
OUT_DIR := output
TEST_DIR:= tests

$(shell $(MKDIR) $(BIN_DIR) $(SRC_DIR) $(OUT_DIR) $(TEST_DIR))

all: 
	@$(MAKE) -j -C build all

test:
	@$(MAKE) -j -C tests

run: all
	./$(BIN_DIR)/scattering_NN_variazionale.x

clean:
	@$(MAKE) -C build clean
	@$(MAKE) -C tests clean