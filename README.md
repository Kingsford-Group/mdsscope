# MDSScope

A Minimum Decycling Set (MDS) is a set of k-mers (word of size k) that is unavoidable and of minimum size.
This code is research code associated with the manuscript ["Sketching methods with small window guarantee using minimum decycling sets"](https://arxiv.org/abs/2311.03592).
Please cite:

``` bibtex
@misc{marçais2023sketching,
      title={Sketching methods with small window guarantee using minimum decycling sets}, 
      author={Guillaume Marçais and Dan DeBlasio and Carl Kingsford},
      year={2023},
      eprint={2311.03592},
      archivePrefix={arXiv},
      primaryClass={cs.DS}
}
```

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
./configure.sh -d 2 4
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

# Canonical space

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

## MDS ILP

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
