# MDSScope

Research code associated with the following two publications.
This code is meant to explore decycling sets, minimum decycling set, the effect of using canonical k-mers and sketching methods for computational biology.

A Minimum Decycling Set (MDS) is a set of k-mers (word of size k) that is unavoidable and of minimum size.
Part of this code is associated with the manuscript ["Sketching methods with small window guarantee using minimum decycling sets"](https://arxiv.org/abs/2311.03592).
Please cite:

``` bibtex
@misc{sketchingmds,
      title={Sketching methods with small window guarantee using minimum decycling sets}, 
      author={Guillaume Marçais and Dan DeBlasio and Carl Kingsford},
      year={2023},
      eprint={2311.03592},
      archivePrefix={arXiv},
      primaryClass={cs.DS}
}
```

The canonical representation $m^c$ of k-mer $m$ is either a $m$ itself or its reverse complement $\overline{m}^r$, whichever is least lexicographically.
Part of this code is associated with the research on the effect of using canonical k-mers with sketching methods in the manuscript [$k$-nonical space: sketching with reverse complements](https://www.biorxiv.org/content/10.1101/2024.01.25.577301).
Please cite:

``` bibtex
@article{knonicalspace,
  author = {Guillaume Mar\c{c}ais and Shane Elder and Carl Kingsford},
  title = {k-nonical space: sketching with reverse complements},
  elocation-id = {2024.01.25.577301},
  year = {2024},
  doi = {10.1101/2024.01.25.577301},
  publisher = {Cold Spring Harbor Laboratory},
  journal = {bioRxiv}
}
```

# Compiling

## An apptainer/singularity

Create an apptainer image with

``` shell
apptainer build mdsscope.sif mdsscope.def
```

The alphabet size and k-mer size are compile time constants, not runtime parameter.
By default the code is compiled with the following parameters: alpha=2, k=6 and alpha=4, k=3.

To change which parameters are used, pass the `params` build argument.
For example, to compile for alpha=2, k=6,7,8 and alpha=4, k=4,5

``` shell
apptainer build --build-arg params="2:6 2:7 2:8 4:4 4:5" mdsscope.sif mdsscope.def
```

To get, for example, the Mykkelteit set for alpha=2 and k=6, run as follows:

``` shell
apptainer run mdsscope.sif /A2K6/mykkeltveit_set
```

See below for a list of programs available.

## Locally with build script

### Requirements

The XXHash library and `g++` >= 12 or `clang++` >= 16 is required.

``` shell
sudo apt install build-essential libxxhash-dev pkg-config
```

The alphabet size and k-mer size are compile time constant (not runtime parameters).
So one needs to compile new executable for every combination of (alphabet, k-mer) they want to use.
For example, to build for the DNA alphabet and k=6:

``` shell
./build.sh 4 6
```

This will create a directory `A4K6` with all the executables in it.

The environmental variables `CXX`, `CXXFLAGS`, `LDFLAGS` and `LDLIBS` are supported to change the compiler frmo the default `g++` or pass extra compiler flags.
For example:

``` shell
CXX=g++-12 CXXFLAGS=-g ./build.sh 4 6
```


## Locally with tup

Building with `tup` is the recommended method to make modification to the code and to run the experiments as in the publications.

### Requirements

The `tup` build system, the XXHash library and `g++` >= 12 or `clang++` >= 16 is required.
On Ubuntu 24.04, install with:

``` shell
sudo apt install tup fuse3 build-essential libxxhash-dev pkg-config
```

### Alphabet and K-mer size

The alphabet size and k-mer size are compile time constants, not runtime parameter.
Hence the programs must be compiled for any pair of alphabet size / k-mer size combination desired.
For example, to generate the programs for the binary alphabet (A=2) and k-mer size 4 (K=4), run:

``` shell
./configure.sh 2 4
tup
```

Then use the programs in the directory `build-A2K4`

To enable debugging (e.g., with GDB), use the `-d` switch:

``` shell
./configure.sh -d 2 4
tup
```

Then use the programs in `build-A2K4-debug`.

Extra compilation flags can be passed to `configure.sh` via the usual environment variables `CXXFLAGS`, `LDFLAGS` and `LDLIBS`.

# Programs

* `mykkeltveit_set`: generate the Mykkeltveit MDS.
* `champarnaud_set`: generate the Champarnaud MDS.
* `syncmer_set`: generate a syncmer set.
* `frac_set`: generate a fractional k-mer set.*
* `comp2rankdot`: given an MDS M, generate a dot plot of all the MDSs reachable from M using F-moves.
* `traverse_components`: given an MDS M, traverse components of the MDS graph using I-moves.
* `fms2mds`: convert a list of F-moves (e.g., as in the output of `traverse_components`) into a corresponding MDS.
* `optimize_rem_path_len`: simulated annealing algorithm to find MDS with minimum or maximum remaining path length.
* `sketch_components`: find the strongly connected components in the de Bruijn graph minus the method's set.
* `sketch_histo`: create a histogram of the distances between selected k-mers for a given set.
* `opt_canon`: greedy optimization procedure to remove k-mers from a set as long as it doesn't create SCCs.
* `create_seed`: create a valid seed file to be used as input to other programs.

Every program has a `--help`.

# Canonical space

The following is not supported by the apptainer, the code must be compiled/run locally.

## Experiments

The `-f` flags triggers running experiments on the sequence given to `f` (the *absolute* path to a file containing the sequence).
The build directory will end in `-exp`, and the results are in the `experiment` subdirectory: see the aggregated files `experiment/{histos,sccs}`.

The sequence file must contain 1 DNA sequence (alphabet = `{A, C, G, T}`) and nothing else (except maybe white space / new lines).
No `N` or other ambiguous bases, no `>` or `@` (not fasta or fastq), no headers, etc.
Only the sequence itself.

Configuration parameters in `config/A?K?-exp.conf`:
* `-f` can be given multiple times to analyse more multiple sequences.
  Alternatively: `CONFIG_EXP_FILES` is a space separated list of file paths.
* `CONFIG_EXP_REPEAT` give the number of random iteration for schemes that are randomized (e.g., syncmers random s-mer order, frac random k-mer order).
* `CONFIG_EXP_HISTO_THRESH` is the threshold for "unmappable" sequences: any sequence of length at least this threshold are considered unmappable.
* `CONFIG_EXP_ILP_MAX_LEN` is the ILP maximum remaining path length parameter.
  See below.
* `CONFIG_EX_TRANSCRIPTS` one (multi-)fasta file.
  Creates an output fasta file with histogram of gaps, one entry for each input entry in the fasta file.
* `CONFIG_EXP_STREAM` if not empty, don't create the sets representing the sketching methods explicitely.
  Create the histograms of gaps in a streaming fashion by selecting the k-mers on the fly.
  Also do not compute the SCCs in the de Bruijn graph.
  Usueful for large values of k.
* `CONFIG_EXP_S` overwrite the value of the s parameter for syncmers from its default of k/2-1

## MDS ILP

This is not well tested and likely not working.

If `MDS-ILP.py -h` works (requires python's modules gurobi, numpy, scipy, igraph, networkx and tqdm), then it will be used to generate a set.
The configuration scripts looks for a `gurobienv` virtual environment in the local directory and sets the `PATH` to use python from that environment.

The maximum default path length with be set to `2*k^3` by default, but can be adjusted with `CONFIG_EXP_ILP_MAX_LEN`.
Recommended setting:

``` shell
pytnon3 -m venv gurobienv
./gurobienv/bin/pip install numpy igraph networkx gurobipy tqdm scipy
```

# Disclaimer

This is research quality code.
We strive for it to be correct for the narrow purpose of our research and publication.
We make no claims of code robustness or usefulness beyond this scope.
Use at your own risk.

<!--  LocalWords:  unmappable
 -->
