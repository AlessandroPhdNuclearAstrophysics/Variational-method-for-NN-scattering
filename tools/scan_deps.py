#!/usr/bin/env python3
import sys                # for command-line arguments
import re                 # for regular expressions
from pathlib import Path  # for path operations (safer and cleaner than os.path)
import os                 # for environment variables

# Get the Fortran source file passed as argument
src_file = Path(sys.argv[1])

# Define root folders for source and build
# Try to read from environment variables if set, otherwise use defaults
src_root = Path(os.environ.get("SRC_ROOT", "src"))
build_root = Path(os.environ.get("BUILD_ROOT", "build"))

# Debug info if needed
debug_mode = os.environ.get("DEBUG_DEPS", "0") == "1"

if debug_mode:
    print(f"Processing file: {src_file}")
    print(f"Source root: {src_root}")
    print(f"Build root: {build_root}")

# Dictionary to map module name → path where it is defined
module_to_path = {}

# Dictionary to track which files contain a main program
program_files = {}

# List of known system modules to ignore
system_modules = ["omp_lib", "iso_c_binding", "iso_fortran_env"]

# Search every Fortran file in src/ to find module definitions and main programs
for path in src_root.rglob("*.f90"):
    try:
        with open(path, "r", encoding="utf-8") as f:
            content = f.read().lower()  # Convert to lowercase for case-insensitive matching
            
            # Check for module definitions
            for match in re.finditer(r"^\s*module\s+([a-z0-9_]+)(?:\s|$)", content, re.MULTILINE):
                modname = match.group(1)
                module_to_path[modname] = path
                if debug_mode:
                    print(f"Found module '{modname}' in {path}")
            
            # Check for program declarations
            if re.search(r"^\s*program\s+([a-z0-9_]+)(?:\s|$)", content, re.MULTILINE):
                program_files[path] = True
                if debug_mode:
                    print(f"Found program in {path}")
    except Exception as e:
        print(f"Warning: Error processing {path}: {e}")

if debug_mode:
    print(f"Found {len(module_to_path)} modules and {len(program_files)} programs")

# Find modules used in the current file
used_modules = []
is_program_file = False

try:
    with open(src_file, "r", encoding="utf-8") as f:
        content = f.read().lower()  # Convert to lowercase for case-insensitive matching
        
        # Check if this is a main program file
        if re.search(r"^\s*program\s+([a-z0-9_]+)(?:\s|$)", content, re.MULTILINE):
            is_program_file = True
            if debug_mode:
                print(f"This file contains a main program")
        
        # Find all used modules
        for match in re.finditer(r"^\s*use\s+([a-z0-9_]+)", content, re.MULTILINE):
            modname = match.group(1)
            if modname not in used_modules:  # Avoid duplicates
                used_modules.append(modname)
                if debug_mode:
                    print(f"This file uses module '{modname}'")
except Exception as e:
    print(f"Error reading {src_file}: {e}")
    sys.exit(1)

# Build path to dependency file
rel_path = src_file.relative_to(src_root).with_suffix(".d")  # e.g. src/libs/utils.f90 → libs/utils.d
dep_file = build_root / "dep" / rel_path
dep_file.parent.mkdir(parents=True, exist_ok=True)

# Build path to object file for the current source
obj_file = build_root / src_file.relative_to(src_root).with_suffix(".o")

# List of missing modules
missing_modules = []

# List of found module objects
found_module_objs = []

# Process dependencies
for modname in used_modules:
    if modname in system_modules:
        # Skip system modules - don't add them as dependencies
        if debug_mode:
            print(f"Skipping system module '{modname}'")
        continue
    elif modname in module_to_path:
        # Where is that module defined?
        mod_path = module_to_path[modname]
        # Where is its object file?
        mod_obj = build_root / mod_path.relative_to(src_root).with_suffix(".o")
        found_module_objs.append(str(mod_obj))
        if debug_mode:
            print(f"Module '{modname}' found in {mod_path} → {mod_obj}")
    else:
        # If the module isn't found, add to missing modules list
        missing_modules.append(modname)
        if debug_mode:
            print(f"Missing module: {modname}")

# Write the .d file (Make-style dependency format)
try:
    with open(dep_file, "w") as f:
        f.write(f"{obj_file}:")
        
        # Add found module objects as dependencies
        for mod_obj in found_module_objs:
            f.write(f" {mod_obj}")
        
        # Add missing modules as comments
        if missing_modules:
            f.write("\n# Missing modules: " + " ".join(missing_modules))
        
        # Mark if this is a main program file
        if is_program_file:
            f.write("\n# MAIN_PROGRAM: True")
        
        f.write("\n")
        
    if debug_mode:
        print(f"Dependency file written to {dep_file}")
except Exception as e:
    print(f"Error writing dependency file {dep_file}: {e}")

if debug_mode:
    print("Dependency processing complete")