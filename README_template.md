[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/bsamseth/cpp-project.svg?branch=master)](https://travis-ci.org/bsamseth/cpp-project)
[![Build status](https://ci.appveyor.com/api/projects/status/g9bh9kjl6ocvsvse/branch/master?svg=true)](https://ci.appveyor.com/project/bsamseth/cpp-project/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/bsamseth/cpp-project/badge.svg?branch=master)](https://coveralls.io/github/bsamseth/cpp-project?branch=master)
[![codecov](https://codecov.io/gh/bsamseth/cpp-project/branch/master/graph/badge.svg)](https://codecov.io/gh/bsamseth/cpp-project)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/eb004322b0d146239a57eb242078e179)](https://www.codacy.com/app/bsamseth/cpp-project?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bsamseth/cpp-project&amp;utm_campaign=Badge_Grade)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/bsamseth/cpp-project.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/bsamseth/cpp-project/context:cpp)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/bsamseth/cpp-project.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/bsamseth/cpp-project/alerts/)
[![license](https://img.shields.io/badge/license-Unlicense-blue.svg)](https://github.com/bsamseth/cpp-project/blob/master/LICENSE)
[![Lines of Code](https://tokei.rs/b1/github/bsamseth/cpp-project)](https://github.com/Aaronepower/tokei)
[![Average time to resolve an issue](http://isitmaintained.com/badge/resolution/bsamseth/cpp-project.svg)](http://isitmaintained.com/project/bsamseth/cpp-project "Average time to resolve an issue")
[![Percentage of issues still open](http://isitmaintained.com/badge/open/bsamseth/cpp-project.svg)](http://isitmaintained.com/project/bsamseth/cpp-project "Percentage of issues still open")

<br />
<p align="center">
  <a href="https://github.com/tractatus/microLRM">
    <img src="https://github.com/tractatus/microLRM/blob/main/logo.png" alt="Logo" width="300">
  </a>
  <p align="center">
     Turn your microscope into a sequencing machine.
    <br />
    <a href="https://github.com/github_username/repo_name"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name">View Demo</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</p>


# What are the aims of µLRM_

This is a boiler plate for C++ projects. What you get:

-   Sources, headers and mains separated in distinct folders
-   Use of modern [CMake](https://cmake.org/) for much easier compiling
-   Setup for tests using [doctest](https://github.com/onqtam/doctest)
-   Continuous testing with [Travis-CI](https://travis-ci.org/), [Appveyor](https://www.appveyor.com) and [GitHub Actions](https://github.com/features/actions), with support for C++17.
-   Code coverage reports, including automatic upload to [Coveralls.io](https://coveralls.io/) and/or [Codecov.io](https://codecov.io)
-   Code documentation with [Doxygen](http://www.stack.nl/~dimitri/doxygen/)

## Structure
``` text
.
├── CMakeLists.txt
├── app
│   └── main.cpp
├── include
│   ├── example.h
│   └── exampleConfig.h.in
├── src
│   └── example.cpp
└── tests
    ├── dummy.cpp
    └── main.cpp
```

Sources go in [src/](src/), header files in [include/](include/), main programs in [app/](app), and
tests go in [tests/](tests/) (compiled to `unit_tests` by default). 

If you add a new executable, say `app/hello.cpp`, you only need to add the following two lines to [CMakeLists.txt](CMakeLists.txt): 

``` cmake
add_executable(main app/main.cpp)   # Name of exec. and location of file.
target_link_libraries(main PRIVATE ${LIBRARY_NAME})  # Link the executable to lib built from src/*.cpp (if it uses it).
```

You can find the example source code that builds the `main` executable in [app/main.cpp](app/main.cpp) under the `Build` section in [CMakeLists.txt](CMakeLists.txt). 
If the executable you made does not use the library in [src/](src), then only the first line is needed.



## Building

Build by making a build directory (i.e. `build/`), run `cmake` in that dir, and then use `make` to build the desired target.

Example:

``` bash
> mkdir build && cd build
> cmake .. -DCMAKE_BUILD_TYPE=[Debug | Coverage | Release]
> make
> ./main
> make test      # Makes and runs the tests.
> make coverage  # Generate a coverage report.
> make doc       # Generate html documentation.
```

## .gitignore

The [.gitignore](.gitignore) file is a copy of the [Github C++.gitignore file](https://github.com/github/gitignore/blob/master/C%2B%2B.gitignore),
with the addition of ignoring the build directory (`build/`).

## Services

If the repository is activated with Travis-CI, then unit tests will be built and executed on each commit.
The same is true if the repository is activated with Appveyor.

If the repository is activated with Coveralls/Codecov, then deployment to Travis will also calculate code coverage and
upload this to Coveralls.io and/or Codecov.io

