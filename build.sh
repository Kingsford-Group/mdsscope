#! /usr/bin/env bash

# This file is based on the file automatically generated with: tup generate
# --config configs/A2K2-debug.config --builddir BUILDDIR build.sh

set -e


usage() {
    echo "Usage: $0 ALPHA K"
}

help() {
    cat <<EOF

Build software for MDS and canonical k-mer exploration. The parameters ALPHA
(the alphabet size) and K (the k-mer length) are required. Extra compilation
flags or option can be passed in the following environment variables:

  CXX             Path to g++ version at least 12
  CXXFLAGS        Compilation flags
  LDFLAGS         Linker flags
  LDLIBS          Extra libraries flags
  PKG_CONFIG_PATH Used by pkg-config

The following programs are built:

traverse_comp
mdss2dot
comp2rankdot
fms2mds
optimize_rem_path_len
mykkeltveit_set
find_longest_path
champarnaud_set
sketch_components
syncmer_set
syncmer_sketch
frac_set
create_seed
sketch_histo
old_champarnaud_set
opt_canon

Options:
  -h          Help, this message
EOF
}

while getopts "h" o; do
    case "$o" in
        (h) usage; help; exit 0 ;;
        (*) usage >&2; exit 1 ;;
    esac
done
shift $((OPTIND-1))

(( $#  != 2 )) && { usage >&2; exit 1; }
ALPHA=$1
K=$2
BUILDDIR="A${ALPHA}K${K}"

set -x

mkdir -p "${BUILDDIR}"/
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/backtrace.o backtrace.cc)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/common.o common.cc)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/sequence.o sequence.cc)
(ar sruv    "${BUILDDIR}"/common.ar "${BUILDDIR}"/backtrace.o "${BUILDDIR}"/common.o "${BUILDDIR}"/sequence.o)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/traverse_comp.o traverse_comp.cc)
(g++ -pthread   "${BUILDDIR}"/traverse_comp.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/traverse_comp)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/mdss2dot.o mdss2dot.cc)
(g++ -pthread   "${BUILDDIR}"/mdss2dot.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/mdss2dot)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/comp2rankdot.o comp2rankdot.cc)
(g++ -pthread   "${BUILDDIR}"/comp2rankdot.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/comp2rankdot)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/fms2mds.o fms2mds.cc)
(g++ -pthread   "${BUILDDIR}"/fms2mds.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/fms2mds)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/optimize_rem_path_len.o optimize_rem_path_len.cc)
(g++ -pthread   "${BUILDDIR}"/optimize_rem_path_len.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/optimize_rem_path_len)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/mykkeltveit_set.o mykkeltveit_set.cc)
(g++ -pthread   "${BUILDDIR}"/mykkeltveit_set.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/mykkeltveit_set)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/find_longest_path.o find_longest_path.cc)
(g++ -pthread   "${BUILDDIR}"/find_longest_path.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/find_longest_path)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/champarnaud_set.o champarnaud_set.cc)
(g++ -pthread   "${BUILDDIR}"/champarnaud_set.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/champarnaud_set)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/sketch_components.o sketch_components.cc)
(g++ -pthread   "${BUILDDIR}"/sketch_components.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/sketch_components)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/syncmer_set.o syncmer_set.cc)
(g++ -pthread   "${BUILDDIR}"/syncmer_set.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/syncmer_set)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/syncmer_sketch.o syncmer_sketch.cc)
(g++ -pthread   "${BUILDDIR}"/syncmer_sketch.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/syncmer_sketch)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/frac_set.o frac_set.cc)
(g++ -pthread   "${BUILDDIR}"/frac_set.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/frac_set)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/create_seed.o create_seed.cc)
(g++ -pthread   "${BUILDDIR}"/create_seed.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/create_seed)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/sketch_histo.o sketch_histo.cc)
(g++ -pthread   "${BUILDDIR}"/sketch_histo.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/sketch_histo)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/old_champarnaud_set.o old_champarnaud_set.cc)
(g++ -pthread   "${BUILDDIR}"/old_champarnaud_set.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/old_champarnaud_set)
(g++    -Wall -DHAVE_EXECINFO_H -I"${BUILDDIR}" -pthread -std=gnu++20 -DHAVE_INT128 -DK="$K" -DALPHA="$ALPHA"   -O0 -g -c -o "${BUILDDIR}"/opt_canon.o opt_canon.cc)
(g++ -pthread   "${BUILDDIR}"/opt_canon.o "${BUILDDIR}"/common.ar   -lxxhash -o "${BUILDDIR}"/opt_canon)
