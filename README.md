# Parallel machine scheduling with decision diagrams

In this repository, you will find the implementation used to produce the results
in the following paper:

We present a new flow-based formulation for identical parallel machine
scheduling with a regular objective function and no idle time. We use a decision
diagram containing all the possible sequences of jobs that follow some ordering
rules to construct the new formulation. These rules do not exclude the optimal
solution of a given instance but constrain the optimal solution to some
canonical form. To define these rules, we need to partition the planning horizon
into non-uniform periods. The novel formulation will have numerous variables and
constraints. Hence, we apply a Dantzig-Wolfe decomposition to compute the linear
programming relaxation of the new flow-based formulation in a reasonable amount
of time.

Moreover, we will see that the resulting lower bound will be stronger than the
lower bound provided by the linear relaxation of the classical time-indexed
formulation. We use a Branch-and-Price framework to solve the new formulation.
Several instances from the literature will be solved for the first time.​

## Instructions

In order to reproduce the results, you need to perform the following steps after
downloading the package.  These instructions will get you a copy of the project
up and running on your local machine for development and testing purposes.

### Prerequisites
* **CMake v3.21+** - found at [https://cmake.org/](https://cmake.org/)

* **C++ Compiler** - needs to support at least the **C++20** standard, i.e.
  *MSVC* at least version 142, *GCC* at least version 10 for linux based systems

* **Vcpkg** - C/C++ dependency manager from Microsoft. For more information on
  how to install and to use vcpkg, we refer to [https://vcpkg.io](https://vcpkg.io).

* **Gurobi** - We use the Gurobi optimizer to compute the linear programming
  (LP) relaxations. You can download Gurobi from the company's
  [https://www.gurobi.com/](website). Follow the installation instructions that
  are provide by gurobi. Per Operating system you will need to adjust
  environment variables such that cmake can find gurobi on your computer. On
  windows for example the environment variables are automatically set.

* **Osi** - Also [https://github.com/coin-or/Osi](Osi) is used. Osi (Open
  Solver Interface) provides an abstract class to generic linear LP, together
  with derived classes for specific solvers. In theory, we can use an arbitrary
  LP solver (Gurobi, CPLEX, XPress, Soplex, ...) to solve the LP relaxation, but
  this is not implemented yet. For now, we can only use gurobi. To install Osi,
  we use coinbrew which can be found
  [https://coin-or.github.io/coinbrew/](here).

### Installing


### Building the project

To build the project, all you need to do, ***after correctly
[installing the project](README.md#Installing)***, is run a similar **CMake** routine
to the the one below:

```bash
mkdir build/ && cd build/
cmake .. -DCMAKE_INSTALL_PREFIX=/absolute/path/to/custom/install/directory
cmake --build . --target install
```

> ***Note:*** *The custom ``CMAKE_INSTALL_PREFIX`` can be omitted if you wish to
install in [the default install location](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html).*

More options that you can set for the project can be found in the
[`cmake/StandardSettings.cmake` file](cmake/StandardSettings.cmake). For certain
options additional configuration may be needed in their respective `*.cmake` files
(i.e. Conan needs the `CONAN_REQUIRES` and might need the `CONAN_OPTIONS` to be setup
for it work correctly; the two are set in the [`cmake/Conan.cmake` file](cmake/Conan.cmake)).

## Generating the documentation

In order to generate documentation for the project, you need to configure the build
to use Doxygen. This is easily done, by modifying the workflow shown above as follows:

```bash
mkdir build/ && cd build/
cmake .. -D<project_name>_ENABLE_DOXYGEN=1 -DCMAKE_INSTALL_PREFIX=/absolute/path/to/custom/install/directory
cmake --build . --target doxygen-docs
```

> ***Note:*** *This will generate a `docs/` directory in the **project's root directory**.*

## Running the tests

By default, the template uses [Google Test](https://github.com/google/googletest/)
for unit testing. Unit testing can be disabled in the options, by setting the
`ENABLE_UNIT_TESTING` (from
[cmake/StandardSettings.cmake](cmake/StandardSettings.cmake)) to be false. To run
the tests, simply use CTest, from the build directory, passing the desire
configuration for which to run tests for. An example of this procedure is:

```bash
cd build          # if not in the build directory already
ctest -C Release  # or `ctest -C Debug` or any other configuration you wish to test

# you can also run tests with the `-VV` flag for a more verbose output (i.e.
#GoogleTest output as well)
```

### End to end tests

If applicable, should be presented here.

### Coding style tests

If applicable, should be presented here.

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our how you can
become a contributor and the process for submitting pull requests to us.

## Authors

* **Daniel Kowalczyk** - [@danielkowalczyk](https://gitlab.kuleuven.be/u0056096)

## License

This project is licensed under the [Unlicense](https://unlicense.org/) - see the
[LICENSE](LICENSE) file for details
