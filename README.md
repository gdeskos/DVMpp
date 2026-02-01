# DVMpp

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/gdeskos/DVMpp)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/96229133.svg)](https://zenodo.org/badge/latestdoi/96229133)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)

**2D Discrete Vortex Method Code written in modern C++17**

A high-performance implementation of the classic Discrete Vortex Method (DVM) with random-walk diffusion for simulating 2D incompressible viscous flows around bluff bodies.

[YouTube Demo: Flow past a cylinder](https://www.youtube.com/watch?v=xckO4AxQIf8&feature=youtu.be)

---

## Features

- **Random-walk diffusion model** based on Chorin (1973)
- **Vortex sheet boundary conditions** following Morgenthal (2002)
- **Multiple time-stepping schemes**: Euler explicit and RK3
- **Surface crossing algorithms**: Delete, Absorb, or Reflect
- **Regularized Biot-Savart kernel** for smooth velocity fields
- **OpenMP parallelization** for performance
- **Optional FMM acceleration** for O(N) velocity computation
- **Comprehensive test suite** with Google Test
- **Modern C++17 codebase** with proper exception handling

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/gdeskos/DVMpp.git
cd DVMpp

# Build
mkdir build && cd build
cmake ..
make -j4

# Run cylinder example
./bin/DVMpp ../examples/cylinder/cylinder.xml
```

---

## Requirements

### Dependencies

| Dependency | Minimum Version | Purpose |
|------------|-----------------|---------|
| CMake | 3.14 | Build system |
| C++ Compiler | C++17 support | GCC 7+, Clang 5+, MSVC 2017+ |
| Armadillo | 7.6 | Linear algebra |
| Boost | 1.56.0 | Program options |
| OpenMP | - | Parallelization |
| BLAS | - | Linear algebra backend |

### Optional Dependencies

| Dependency | Purpose |
|------------|---------|
| Google Test | Testing (fetched automatically) |
| exafmm-t | Fast Multipole Method acceleration |
| Doxygen | Documentation generation |

### Installation (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install -y \
    cmake \
    g++ \
    libarmadillo-dev \
    libboost-program-options-dev \
    libomp-dev \
    libblas-dev \
    liblapack-dev
```

### Installation (macOS with Homebrew)

```bash
brew install cmake armadillo boost libomp
```

---

## Building

### Standard Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `BUILD_TESTING` | ON | Build test suite |
| `USE_FMM` | OFF | Enable Fast Multipole Method |
| `BUILD_DOC` | ON | Build Doxygen documentation |

### Build with FMM Support

```bash
cmake .. -DUSE_FMM=ON
make -j$(nproc)
```

### Build with Tests

```bash
cmake .. -DBUILD_TESTING=ON
make -j$(nproc)
ctest --output-on-failure
```

### Release Build

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

---

## Usage

### Command Line

```bash
# Run a simulation
./bin/DVMpp input.xml

# Generate example XML configuration
./bin/DVMpp --example-xml example.xml

# Print help
./bin/DVMpp --help
```

### XML Configuration

```xml
<?xml version="1.0"?>
<DVMpp>
    <algorithms>
        <surface_crossing string="REFLECT" />
        <!-- Options: DELETE, ABSORB, REFLECT -->
    </algorithms>
    <constants>
        <density val="1.0" />
        <nu val="0.001" />              <!-- Kinematic viscosity -->
        <max_NumPanelVort val="3" />
        <cutoff_exp val="0.675" />
        <seed val="-1" />               <!-- -1 for random, >=0 for fixed -->
    </constants>
    <io>
        <input_dir string="./input/" />
        <output_dir string="./output/" />
        <domain_file string="body.txt" />
    </io>
    <flow>
        <ux val="3.0" />
        <uz val="0.0" />
    </flow>
    <probe>
        <x val="20" />
        <z val="1" />
    </probe>
    <time>
        <scheme string="euler" />       <!-- Options: euler, RK3 -->
        <dt val="0.005" />
        <steps val="1000" />
    </time>
</DVMpp>
```

### Body Geometry File

The body geometry is specified as a list of (x, z) coordinates defining the closed contour:

```
0.5 0.0
0.475 0.154
...
0.5 0.0
```

---

## Output Files

| File | Description |
|------|-------------|
| `*_vortex.dat` | Vortex blob positions and circulations |
| `*_vortex_num.dat` | Number of vortices per timestep |
| `*_gamma.dat` | Surface vorticity distribution |
| `*_loads.dat` | Force coefficients (Cd, Cl) |
| `*_probe.dat` | Velocity at probe points |
| `*_xml_in.xml` | Copy of input configuration |

---

## Testing

```bash
# Build and run all tests
cd build
cmake .. -DBUILD_TESTING=ON
make -j$(nproc)
ctest --output-on-failure

# Run specific test categories
ctest -R unit          # Unit tests only
ctest -R integration   # Integration tests only
ctest -R regression    # Regression tests only

# Run with verbose output
ctest -V
```

---

## Project Structure

```
DVMpp/
├── CMakeLists.txt          # Root CMake configuration
├── README.md               # This file
├── src/                    # Source code
│   ├── CMakeLists.txt
│   ├── BaseTypes.hpp       # Type definitions
│   ├── DVMBase.hpp/cpp     # Main DVM solver
│   ├── VortexBlobs.hpp/cpp # Point vortices
│   ├── VortexSheet.hpp/cpp # Boundary elements
│   ├── VelocityKernel.hpp/cpp # Biot-Savart kernels
│   ├── Exceptions.hpp      # Exception hierarchy
│   ├── Probe.hpp/cpp       # Velocity probes
│   ├── Random.hpp/cpp      # RNG wrapper
│   ├── XmlHandler.hpp/cpp  # XML configuration
│   └── pugi*.hpp/cpp       # PugiXML library
├── tests/                  # Test suite
│   ├── unit/               # Unit tests
│   ├── integration/        # Integration tests
│   ├── regression/         # Regression tests
│   └── fixtures/           # Test data
├── examples/               # Example configurations
│   └── cylinder/
├── doc/                    # Documentation
│   ├── CMakeLists.txt
│   └── Doxyfile
└── bin/                    # Compiled executables
```

---

## Citation

If you use DVMpp in your research, please cite:

```bibtex
@software{dvmpp,
  author = {Deskos, Georgios},
  title = {DVMpp: 2D Discrete Vortex Method in C++},
  year = {2017},
  publisher = {GitHub},
  url = {https://github.com/gdeskos/DVMpp},
  doi = {10.5281/zenodo.XXXXXX}
}
```

---

## References

1. Chorin, A. J. (1973). Numerical study of slightly viscous flow. *Journal of Fluid Mechanics*, 57, 785-796.

2. Morgenthal, G. (2002). *Aerodynamic Analysis of Structures Using High-Resolution Vortex Particle Methods*. PhD thesis, University of Cambridge.

3. Perlman, M. (1985). On the accuracy of vortex methods. *Journal of Computational Physics*, 59(2), 200-223.

4. Smith, P. A., & Stansby, P. K. (1989). An efficient surface algorithm for random-particle simulation of vorticity and heat transport. *Journal of Computational Physics*, 81(2), 349-371.

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

## Acknowledgments

- Developed at Imperial College London (2015-2017)
- Based on classical DVM theory by Chorin and extensions by Morgenthal
