# [SiPaR-SiskinCC-RobinCC]()

### Documentation
**[docs/readme.md](docs/readme.md)**

### Cloning 
`git clone https://github.com/MohsenKoohi/SiPaR-SiskinCC-RobinCC.git --recursive`

### Updating through `pull`
- `git pull --recurse-submodules` or
- You may set recursive submodule update globally using `git config --global submodule.recurse true` and then `git pull` fetches all updates.

### Requirements
1. Libraries: `mpi`, `openmp`, and `libnuma`
2. Compiler: `gcc` with a version greater than 9 are required.
3. For using ParaGrapher: `JDK` with a version greater than 15 and `libfuse`.
4. Bash tools: `unzip`, `bc`,  and `wget`.

### Compiling and Executing Code
 - Make sure the requried libraries are accessible through `$LD_LIBRARY_PATH`.
 - Allocating machines in SLURM: 
 - Run `make alg...` (e.g. `make alg1_sapco_sort`). This builds the executible file and runs it for the test graph. 
 - For identifying input graph and other options, please refer to [Loading Graphs](docs/0.2-loading.md).
 
### Supported Graph Types & Loading Graphs
 - [CompBin format](https://doi.org/10.48550/arXiv.2507.00716), 
 using [WG2CompBin library](https://github.com/MohsenKoohi/WG2CompBin)
 - [WebGraph format](https://webgraph.di.unimi.it/), using
[ParaGrapher](https://github.com/MohsenKoohi/ParaGrapher) library.

Please refer to [Graph Loading Documentation](docs/0.2-loading.md).

### Evaluating a Number of Graph Datasets

Please refer to [Launcher Script Documentaion](docs/0.3-launcher.md).
 
### Bugs & Support

If you receive wrong results or you are suspicious about parts of the code, 
please [contact us](https://orcid.org/0000-0002-7465-8003).

### License

Licensed under the GNU v3 General Public License, as published by the Free Software Foundation. 
You must not use this Software except in compliance with the terms of the License. 
Unless required by applicable law or agreed upon in writing, this Software is distributed 
on an "as is" basis, without any warranty; without even the implied warranty of 
merchantability or fitness for a particular purpose, neither express nor implied. 
For details see terms of the License (see attached file: LICENSE). 

#### Copyright 2025
