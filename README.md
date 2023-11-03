# MDSScope

A Minimum Decycling Set (MDS) is a set of k-mers (word of size k) that is unavoidable and of minimum size.
This code is research code associated with the manuscript "Sketching methods with small window guarantee using minimum decycling sets".

# Compiling

## Requirements

The `tup` build system and the `yaggo` argument parsing library are required.
Both are available from Ubuntu repositories:

``` shell
sudo apt install tup yaggo
```

For the compilation commands to work, initialize `tup` in the directory (needed to be run only once): `tup init`.

## Alphabet and K-mer size

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
./configure.sh -h 2 4
tup
```

Then use the programs in `build-A2K4-debug`.

Extra compilation flags can be passed to `configure.sh` via the usual environment variables `CXXFLAGS`, `LDFLAGS` and `LDLIBS`.

# Programs

* `mykkeltveit_set`: generate the Mykkeltveit MDS.
* `champarnaud_set`: generate the Champarnaud MDS.
* `comp2rankdot`: given an MDS M, generate a dot plot of all the MDSs reachable from M using F-moves.
* `traverse_components`: given an MDS M, traverse components of the MDS graph using I-moves.
* `fms2mds`: convert a list of F-moves (e.g., as in the output of `traverse_components`) into a corresponding MDS.
* `optimize_rem_path_len`: simulated annealing algorithm to find MDS with minimum or maximum remaining path length.

# Disclaimer

This is research quality code.
We strive for it to be correct for the narrow purpose of our research and publication.
We make no claims of code robustness or usefulness beyond this scope.
Use at your own risk.
