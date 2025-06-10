# Two-Nucleon Scattering 

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Language](https://img.shields.io/badge/language-Fortran-orange.svg)

A comprehensive Fortran framework for computing nucleon-nucleon scattering phase shifts using variational methods with different potential models.

## Overview

This project implements a variational approach to nuclear scattering calculations, providing tools to:

- Calculate phase shifts and mixing angles for nucleon-nucleon scattering
- Transform phase shifts to k²*cot(δ) format for effective range analysis
- Evaluate low-energy scattering observables and extract scattering lengths and effective ranges
- Support multiple potential models (AV14, AV18, EFT_pless)

The code handles all physical nucleon-nucleon channels, including both uncoupled (e.g., ¹S₀, ³P₀) and coupled channels (e.g., ³S₁-³D₁, ³P₂-³F₂).

## Prerequisites

### Compiler Requirements
- Fortran compiler (gfortran ≥ 8.0 recommended)
- OpenMP support

### Required Libraries
- GSL (GNU Scientific Library) for special functions
- LAPACK/BLAS for linear algebra operations

### Development Tools
- Python 3 (for dependency scanning)
- make (for build system)

### Optional Tools
- xmgrace (for plotting results)
- Doxygen (for generating documentation)

### Developer Tools
- `tools/scan_deps.py`: Scans Fortran source files to map module dependencies and main programs. Used in the build system for automatic dependency tracking.
- `tools/dependency_tree.py`: Visualizes module dependencies as tree diagrams (text and graphical output). Consumes outputs from `scan_deps.py` and can generate PNG diagrams in the `dependency_graphs/` directory.

#### Example Usage (from the Makefile)
To generate dependency files for a Fortran source file:
```bash
python3 tools/scan_deps.py src/main_scattering_NN_variazional_method.f90
```

To generate dependency graphs (as used in the Makefile's `generate_dependency_graphs` target):
```bash
python3 tools/dependency_tree.py build/dep -d dependency_graphs --auto-adjust -l circular -f -g -p 3 -s 50 --remove-transitive
```
This will analyze all dependency files in `build/dep/` and generate diagrams in `dependency_graphs/`.

### Installation on Debian/Ubuntu Systems
```bash
sudo apt-get install gfortran libgsl-dev liblapack-dev libblas-dev python3 grace doxygen
```

### Directory Structure
```
.
├── bin/                # Shell scripts for running and processing the code
│   ├── runner.sh
│   └── fitter.sh
│
├── src/                # Source code
│   ├── libs/           # Library modules
│   │   ├── math/           # Mathematical utilities (integration, special functions, etc.)
│   │   ├── physics/        # Physical quantities and constants
│   │   ├── potentials/     # Nuclear potential implementations
│   │   └── utils/          # Utility and helper functions
│   │
│   ├── main_*.f90      # Main program files
│   └── tests/          # Test programs
│
└── Makefile            # Build configuration
```

## Building the Project

The project uses a Makefile with several options and automation features:

### Standard Build (Optimized)
```bash
make
```

### Debug Build
```bash
make DEBUG=1
```

### With Profiling Information
```bash
make GMON=1
```

### Building Documentation
```bash
make doc
```

### Running Tests
```bash
make test
```

### Cleaning Build Artifacts
```bash
make clean
```
Removes all build files, logs, and generated graphs.

### Deleting Output Files
```bash
make delete_out
```
Removes all files in the output directory.

### Generating Dependency Graphs
```bash
make generate_dependency_graphs
```
Creates graphical representations of module dependencies in the `dependency_graphs/` directory.

### Library Management and Export

The Makefile can build a static library (`libvariational.a`) from all modules:
```bash
make build/libvariational.a
```

To export the static library **together with all module files** (required for use in other Fortran projects), you can create a zip archive:
```bash
make build/libvariational.zip
```
This will package `libvariational.a` and all `.mod` files from the build directory into `build/libvariational.zip` for easy distribution and reuse.

### Additional Features
- Automatic dependency tracking and correct compilation order using generated `.d` files.
- Support for static library creation and export.
- Automated test execution and logging.
- Automated documentation generation with Doxygen.
- Profiling support with gprof.
- **OpenMP parallelization is always enabled by default.**

## Usage

### Basic Execution
```bash
./bin/runner.sh
```

### With Specific Parameters
```bash
./bin/runner.sh -emax 50.0 -ne 500 -ipot 18 -out_dir output/AV18_50MeV/
```

#### Parameters:
- `-emax`: Maximum energy in MeV
- `-ne`: Number of energy points
- `-tz`: Isospin (0 for np, 1 for pp/nn)
- `-ipot`: Potential model (14 for AV14, 18 for AV18, etc.)
- `-ilb`: Potential submodel
- `-lemp`: Electromagnetic potential flag (0 for pure Coulomb, 1 for full EM)
- `-print`: Print info during calculations
- `-out_dir`: Output directory

### Processing Phase Shifts
```bash
./bin/fitter.sh output/AV18_50MeV/
```

### Transforming Data Format
```bash
./build/main_transform_to_k2_kcotd.x output/AV18_50MeV/
```

## Output Files

### Phase Shift Files
Located in your specified output directory (default: `output/AV18/`):
- `delta_1S0.dat`: Phase shifts for ¹S₀ channel
- `delta_3P0.dat`: Phase shifts for ³P₀ channel
- `delta_3S1-3D1.dat`: Phase shifts and mixing angle for ³S₁-³D₁ coupled channel
- etc.

Format: `energy  delta1  [delta2  epsilon]` (columns depend on channel type)

### k²*cot(δ) Files
Generated by the transform program:
- `k2_kcotd_1S0.dat`
- `k2_kcotd_3S1.dat`
- etc.

Format: `k²  k^(2L+1)*cot(delta)`

### Analysis Files
Generated by the fitter program in the `fits/` directory:
- Effective range parameters
- Polynomial expansion coefficients

## Mathematical Background

The code implements the variational approach to scattering, which:

1. Expands the wave function in a suitable basis
2. Minimizes the functional to find phase shifts
3. Uses Laguerre polynomials as the expansion basis
4. Handles different nuclear potentials through a modular interface

### Effective Range Expansion

For uncoupled channels, the effective range expansion is implemented as:
```
k^(2L+1)cot(δ) = Σ c_i k^(2i)
```
From this expansion, we extract:
- Scattering length: a = -1/c₀
- Effective range: r_e = 2c₁

For coupled channels, mixing angles are represented by:
```
ε_J = Σ e_i k^(2(ΔL+i))
```

### Numerical Integration

The code uses several numerical integration methods:
- Block-adaptive integration for radial functions
- Gauss-Laguerre quadrature for exponentially decaying functions
- Exponentially growing grids for high resolution near the origin

## To-Do List
- [x] Solve the problems with `OMP` parallelization
- [x] Implement possibility to fit new potential models using calculating the potential radial functions once
- [x] Optimize performance for large-scale calculations
- [ ] Implement "nn" and "pp" cases in the same way as "np" (currently only np is implemented)
- [ ] Implement the `lemp=1` option to include the electromagnetic potential in the calculations and the correct asymptotic behavior of the wave function
- [ ] Implement additional potential models (e.g., EFT_pionfull)
- [ ] Add more test cases for robustness
- [ ] Improve documentation and examples

## License

You are welcome to use this project and contribute improvements, provided that you collaborate by sharing your enhancements with the community. If you use this project or its results in any publication, please acknowledge and cite the original author.

## Contact

Alessandro Grassi

alessandro.grassi@df.unipi.it