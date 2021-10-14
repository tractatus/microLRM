[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/bsamseth/cpp-project.svg?branch=master)](https://travis-ci.org/bsamseth/cpp-project)
[![Build status](https://ci.appveyor.com/api/projects/status/g9bh9kjl6ocvsvse/branch/master?svg=true)](https://ci.appveyor.com/project/bsamseth/cpp-project/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/bsamseth/cpp-project/badge.svg?branch=master)](https://coveralls.io/github/bsamseth/cpp-project?branch=master)
[![codecov](https://codecov.io/gh/bsamseth/cpp-project/branch/master/graph/badge.svg)](https://codecov.io/gh/bsamseth/cpp-project)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/eb004322b0d146239a57eb242078e179)](https://www.codacy.com/app/bsamseth/cpp-project?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bsamseth/cpp-project&amp;utm_campaign=Badge_Grade)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/tractatus/microLRM.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/tractatus/microLRM/context:cpp)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/tractatus/microLRM.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/tractatus/microLRM/alerts?mode=list)
[![license](https://img.shields.io/badge/license-Unlicense-blue.svg)](https://github.com/bsamseth/cpp-project/blob/master/LICENSE)
[![Lines of Code](https://tokei.rs/b1/github/tractatus/microlrm)](https://github.com/Aaronepower/tokei)
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
    <a href="https://github.com/tractatus/microLRM"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/tractatus/microLRM">View Demo</a>
    ·
    <a href="https://github.com/tractatus/microLRM/issues">Report Bug</a>
    ·
    <a href="https://github.com/tractatus/microLRM/issues">Request Feature</a>
  </p>
</p>

<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="https://github.com/tractatus/microLRM/blob/main/CODE_OF_CONDUCT.md">Code of conduct</a>
    </li>
    <li>
      <a href="#about-the-project">What are the aims of µLRM?</a>
    </li>
    <li>
      <a href="#dependencies">Dependencies</a>
    </li>
    <li><a href="#usage">How is µLRM input and output organized?</a>  
      <ol>
      <li><a href="#output">Output</a>
      <li><a href="#bcl">BCL file format</a></li>
      <li><a href="#bcl">LOCS file format</a></li>
      </ol>

    </li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>


# About the project

What are the aims of µLRM?

µLRM or microLRM is an open source alternative to Illumina's [Local Run Manager](https://www.illumina.com/products/by-type/informatics-products/local-run-manager.html). µLRM is a compiled software for real-time microscope and fluidics control that turns any microscope supported by micro-manager's [mmCoreAndDevices](https://github.com/micro-manager/mmCoreAndDevices) into a sequencing machine for _in situ_ sequencing.

With micro Local Run Manager you will be able to: create, monitor, and analyze microscope sequencing runs.  

## Dependencies

* [mmCoreAndDevices](https://github.com/micro-manager/mmCoreAndDevices)
* [OpenCV](https://github.com/opencv/opencv)
* [imgui](https://github.com/ocornut/imgui)
* [Eigen](https://eigen.tuxfamily.org/)
* [TensorRT](https://github.com/NVIDIA/TensorRT)
* [PolyScope](https://polyscope.run/)


## How is µLRM input and output organized?

### Output

Output format closely follows the file folder structure of most contemporary sequencing machines.

Traditionally a "Lane" in Next Generation Sequencing context is the physical lane on a flow cell. 

Lanes might be used to separate different experiments through physical separation of the sequencing reaction. 

Splitting up the analysis based on lanes has traditionally been important for QC purposes in figuring out when fluidics for a given flow cell lane is failing.

However, for _in situ_ sequencing context the lane, or rather the [reaction well](https://gracebio.com/products/hybridization-and-incubation/hybriwell-sealing-system-hybridization-and-incubation/) or [fluidic wall](https://www.pnas.org/content/115/26/E5926.short),  is not the only important level of analysis since realtime sequencing performance is also dependent on the tissue environment of the single sample being sequenced. 

Rather than structuring up into more subfolders a "lane" here is considered to be an individual piece of tissue covered by a set of field of views. Note that with this definition several "lanes" can belong to a single reaction well. To further designate which tissue pieces belong to the same reaction a `📝 LaneGroups.json` file provides information on the groupings of "lanes".

``` text
.
├── 📂 Data
│   ├── 📂 Intensities
│   │   ├── 📂 BaseCalls
│   │   │   └── 📂 L001
│   │   │       └── 📂 C1.1
│   │   │           └── 📝 s_1_1.bcl
│   │   └── 📂 L001
│   │       └── 📝 s_1_1.locs
│   └── 📂 RTALogs
├── 📂 InterOp
│   └── 📂 cache
├── 📂 Logs
├── 📂 Recipe
├── 📝 RunInfo.json
├── 📝 RunParameters.json
├── 📝 LaneGroups.json
├── 📝 SampleSheet.csv
└── 📂 Thumbnail_images
    └── 📂 L001
        └── 📝 s_1_1.html
```
At the lowest level base calls for a given cycle, e.g. 📂 C1.1, lane, e.g.  📂 L001, and field of view (FOV) is represented by a single `.bcl` file.

Locations of individual amplicons and their point set registration results across cycles are summarized for a given lane, e.g.  📂 L001, and FOV by a single `.locs` file.

#### BCL file format
The binary base call (BCL) sequence file format is a binary format that can easily be converted to a human readable FASTQ file. 
µLRM software writes the base and the confidence in the call as a quality score to base call (.bcl) files. This is done in real time, i.e. for every cycle of the sequencing run a call for every location identified on the flow cell is added. BCL files are stored in binary format and represent the raw data output of a sequencing run.

BCL files are compressed with the [gzip (*.gz)](https://www.gnu.org/software/gzip/) or [blocked GNU zip (`*.bgzf`)](https://github.com/lh3/samtools/blob/master/bgzf.h) format.
Blocked gzip files are larger in size but [improves](https://blastedbio.blogspot.com/2011/11/bgzf-blocked-bigger-better-gzip.html) random access.

*Table 1. Byte specification of the BCL format*
| Bytes                         | Description                                                                                                                                                              | Type   |
|-------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------|
| 0–3                     | Number of N cluster                                                                                                                                                      | uint32 |
| 4-(N+3) <br>N-Cluster index | Bits 0–1 are the bases, [A, C, G, T] for [0, 1, 2, 3]. Bits 2–7 are shifted by 2 bits and contain the quality score. All bits with 0 in a byte are reserved for no call. | uint8  |

#### LOCS file format
The `locs` file format stores position data exclusively. `locs` files store position data for successive clusters in 4 byte float pairs, described as follows:

*Table 2. Byte specification of the LOCS format*
| Bytes 	| Description                   	| Type  	|
|-------	|-------------------------------	|-------	|
| 0-3   	| Version number                	| uint8   |
| 4-7   	| 2.0                           	| double 	|
| 8-11  	| Number of clusters            	| uint32  |
| 12-15 	| X coordinate of first cluster 	| double 	|
| 16-19 	| Y coordinate of first cluster 	| double 	|
| 20-23 	| Z coordinate of first cluster 	| double 	|

The remaining bytes of the file store the X, Y and Z coordinates of the remaining clusters.

### Input

µLRM, unlike traditional fluorescent NGS platforms, is specifically designed to handle 3D data acquired through z-stacks a tiles with confocal quality. That is not to say that 2D data cant be used as input. However the default trained neural networks differ from 3D and 2D scenario. See segmentation and base calling for more details.

## Structure of build
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

Additional sources go in [src/](src/), header files in [include/](include/), main programs in [app/](app), and
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
<!--
## Services

If the repository is activated with Travis-CI, then unit tests will be built and executed on each commit.
The same is true if the repository is activated with Appveyor.

If the repository is activated with Coveralls/Codecov, then deployment to Travis will also calculate code coverage and
upload this to Coveralls.io and/or Codecov.io
-->

