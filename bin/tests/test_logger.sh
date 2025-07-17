#!/bin/bash
# Test suite for logger.f90 module
# This script runs Fortran test programs for the logger module and checks basic output

set -e

TEST_DIR="$(dirname "$0")"
LOG_TEST_EXE="$TEST_DIR/../../build/release/tests/test_logger.x"
LOG_OUTPUT_FILE="$TEST_DIR/test_logger_output.txt"

# Use make to compile all
make

# Run the test program
if ! "$LOG_TEST_EXE" > "$LOG_OUTPUT_FILE"; then
    echo "Test program failed!"
    exit 1
fi

# Check for expected output
if  grep -a -q "This is an error message" "$LOG_OUTPUT_FILE" && \
    grep -a -q "This is a warning message" "$LOG_OUTPUT_FILE" && \
    grep -a -q "This is an info message" "$LOG_OUTPUT_FILE" && \
    grep -a -q "This is a debug message" "$LOG_OUTPUT_FILE"; then
    echo -e "\e[32mPASSED\e[0m: Logger test passed."
else
    echo -e "\e[31mFAILED\e[0m: Logger test failed: expected log messages not found."
    exit 1
fi

rm -vf "$LOG_OUTPUT_FILE"
