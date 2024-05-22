# Concentric Interpolation

## Content of this file

1. [What's this](#whats-this)
2. [How to use](#how-to-use)
3. [Code style and paradigm](#code-style-and-paradigm)
4. [Bug reporting](#bug-reporting)
5. [License information](#license-information)

## What's this

The Concentric Interpolation method efficiently interpolates multi-dimensional data on multi-dimensional spaces. It is based on a concentric placement of the function samples, i.e. centered around the origin and placed along equidistributed directions at certain radii. This structure is exploited by means of a two-staged interpolation: first, a radial interpolation is performed, i.e. a one-dimensional interpolation along each direction by means of piecewise defined polynomial. Second, a tangential interpolation is performed using classical spherical basis functions. For more information, see the associated research article:

|||
|---------|-------------|
| Authors | Oliver Kunc and Felix Fritzen |
| Title   | Generation of energy-minimizing point sets on spheres and their application in mesh-free interpolation and differentiation |
| Journal | Advances in Computational Mathematics |
| Year    | 2019 |
| Volume  | 45 |
| Issue   | 5-6 |
| Pages   | 3021-3056 |
| DOI     | [10.1007/s10444-019-09726-5](https://doi.org/10.1007/s10444-019-09726-5) |

## How to use

### Requirements

A C++ compiler (e.g. g++), the CMake utility, the GNU GSL library as well as an implementation of LAPACKE are required. In Debian or Ubuntu systems, the command

```bash
sudo apt install g++ cmake liblapacke liblapacke-dev libgsl-dev libgsl23 libgslcblas0
```

should suffice.

### Documentation

Build command for the documentation:

```bash
doxygen ConcentricInterpolation.doxy
```

This requires doxygen software which can be installed by running

```bash
sudo apt install doxygen
```

on Debian or Ubuntu systems. Then, see doxygen output in `./doxygen/html/index.html` or in `./doxygen/latex/*.pdf`. See inline comments for further details.

### Compilation and Example

For the library, in the project root run (1. configure cmake and generate build files, 2. perform actual build)

```bash
cmake --preset ci-linux [-DFLAGS...=...]
cmake --build --preset ci-linux
```

This builds the project. The .a binary will be located in `build/release`. By default a static library is built. To make the library shared, use the flag `-DCONINTER_BUILD_SHARED=ON`.

Optional: To make the library available systemwide you can install it by running  (add `--prefix <install_dir>` for a test install):

```bash
sudo cmake --install build/release
```

Alternatively you can use the CMake feature `FetchContent` and refer to this git repository to pull the code and compile it as part of your project (See CMake documentation for details).

Independent of the way you chose to incorporate the library into your project you can use `find_package(ConInter)` and `target_link_libraries(tgt ConInter)` to find and link to the library.

For the examples, use the flag `-DCONINTER_BUILD_EXAMPLES=ON`. The executables will be in `build/release/examples`.

### How to Cite

When referring to this software, please cite the article mentioned above in Section 1 [What's this](#whats-this).

## Code style and paradigm

- Indentation by 4 spaces, no tabs.
- Camel case function names, variable names less strict.
- No trailing whitespaces.
- Function arguments begin with "a_" in order to strictly separate from local identifiers. Local identifiers do not begin with "a_".
- Error handling not rigorous yet. Unclear where to go from here, depends on yet unknown user requirements.
- Working on low-level pointer objects: the user is required to assure that memory bounds are respected. Setup functions and destructors take care of the basic management, but access to components is unrestricted. Read or write access is controlled via the respective Get... or Set... functions.
- Const correctness: top-level pointers are not const, all lower levels should be const. This is mostly but not 100% consistent atm.
- For essential conditional checks in functions, use UTILITY::assert_msg()
- If a scope (e.g. body of a loop) consists of only one line, don't use braces "{}" but only indentation.

## Bug reporting

Please report bugs and suggestions to the administrators of the corresponding github repository: [https://github.com/EMMA-Group/ConcentricInterpolation](https://github.com/EMMA-Group/ConcentricInterpolation) or directly e-mail [fritzen@mechbau.uni-stuttgart.de](mailto:fritzen@mechbau.uni-stuttgart.de) or [kunc@mechbau.uni-stuttgart.de](mailto:kunc@mechbau.uni-stuttgart.de)

## License information

See file [LICENSE](LICENSE)
