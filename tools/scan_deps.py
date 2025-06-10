#!/usr/bin/env python3
import sys                # for command-line arguments
import re                 # for regular expressions
from pathlib import Path  # for path operations (safer and cleaner than os.path)

# Get the Fortran source file passed as argument
src_file = Path(sys.argv[1])

# Define root folders for source and build
src_root = Path("src")
build_root = Path("build")

# Dictionary to map module name → path where it is defined
module_to_path = {}

# List of known system modules to ignore
system_modules = ["omp_lib", "iso_c_binding", "iso_fortran_env"]

# Search every Fortran file in src/ to find module definitions
for path in src_root.rglob("*.f90"):
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            # Match lines like: module utils
            match = re.match(r"\s*module\s+([a-zA-Z0-9_]+)", line, re.IGNORECASE)
            if match:
                modname = match.group(1).lower()
                module_to_path[modname] = path

# Find modules used in the current file
used_modules = []
with open(src_file, "r", encoding="utf-8") as f:
    for line in f:
        # Match lines like: use utils
        match = re.match(r"\s*use\s+([a-zA-Z0-9_]+)", line, re.IGNORECASE)
        if match:
            used_modules.append(match.group(1).lower())

# Build path to dependency file
rel_path = src_file.relative_to(src_root).with_suffix(".d")  # e.g. src/libs/utils.f90 → libs/utils.d
dep_file = build_root / "dep" / rel_path
dep_file.parent.mkdir(parents=True, exist_ok=True)

# Build path to object file for the current source
obj_file = build_root / src_file.relative_to(src_root).with_suffix(".o")

# Write the .d file (Make-style dependency format)
with open(dep_file, "w") as f:
    f.write(f"{obj_file}:")
    for mod in used_modules:
        if mod in system_modules:
            # Skip system modules - don't add them as dependencies
            continue
        elif mod in module_to_path:
            # Where is that module defined?
            mod_path = module_to_path[mod]
            # Where is its object file?
            mod_obj = build_root / mod_path.relative_to(src_root).with_suffix(".o")
            f.write(f" {mod_obj}")
        else:
            # If the module isn't found, add as a separate comment line, not inline
            f.write(f"\n# Missing module: {mod}")
    f.write("\n")