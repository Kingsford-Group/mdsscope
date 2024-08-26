#!/bin/bash

set -e

usage() {
    echo "Usage: $0 [options] ALPHA K" >&2
}

help() {
    usage
    cat <<EOF

Generate a tup configuration in configs/name.config, and create the tup variant.
ALPHA (the alphabet size) and K (the k-mer length) are required. Extra
compilation flags or option can be passed in the following environment
variables:

  CXX             Path to g++ version at least 12
  CXXFLAGS        Compilation flags
  LDFLAGS         Linker flags
  LDLIBS          Extra libraries flags
  YAGGO           Path to yaggo
  PKG_CONFIG_PATH Used by pkg-config

Enable testing the -r and -f flags. -r gives the number of repeats for an
experiment that is randomized (e.g., syncmer sketches). -f gives the sequence
files to work on.

Options:
  -d          Debugging exec
  -j          Force ^j^ flag to create compiledb
  -r          Number of repeats
  -f          Testing sequence file
  -n          Dry run
  -h          Help
EOF
}

OPTFLAGS="-O3 -DNDEBUG"
SUFFIX=
COMPILEDB=
REPEAT=1
FILES=()
DRYRUN=
ILP=

while getopts "djr:f:inh" o; do
    case "${o}" in
        (d) OPTFLAGS="-O0 -g"; SUFFIX="-debug" ;;
        (h) help; exit 0 ;;
        (j) COMPILEDB=yes ;;
        (r) REPEAT=$OPTARG ;;
        (f) FILES+=("$(realpath $OPTARG)") ; SUFFIX="-exp" ;;
        (i) ILP=1 ;;
        (n) DRYRUN=1 ;;
        (*) usage; exit 1 ;;
    esac
done
shift $((OPTIND-1))

ALPHA=$1
K=$2

if ! grep -qP '^\d+$' <<< "$ALPHA" || ! grep -qP '^\d+$' <<< "$K" || ! grep -qP '^\d+$' <<< "$REPEAT"; then
    usage
    exit 1
fi

NAME="A${ALPHA}K${K}${SUFFIX}"

[ -z "$YAGGO" ] && YAGGO=$(which yaggo || true)
[ -z "$TUP" ] && TUP=$(which tup || true)

[[ -z "$TUP" || -z "$YAGGO" ]] && { echo >&2 "Missing required dependencies: tup and/or yaggo"; false; }

detect_compiledb() {
  tup compiledb >& /dev/null && echo yes || true
}

[ -z "$COMPILEDB" ] && COMPILEDB="$(detect_compiledb)"

# Detect python for ilp_set
ILPPYTHON=
if [ -n "$ILP" ]; then
for p in "$(realpath -s gurobienv/bin/python)" "$(which python3)"; do
  if [ -x "$p" ] && "$p" MDS-ILP.py -h &> /dev/null; then
    ILPPYTHON=$p
    break
  fi
done
  [ -z "$ILPPYTHON" ] && echo >&2 "No ILP: didn't find a satisfying Python interpreter and packages"
fi

# Find a valid version for g++
GCXX=
for gcc in $CXX g++ g++-12; do
    $gcc -o check_gcc_version -O0 check_gcc_version.cc
    ./check_gcc_version 12 0 && GCXX=$gcc && break
done
rm -f check_gcc_version
[ -z "$GCXX" ] && { echo >&2 "Didn't find g++ version at least 12.0"; false; }

# Finc xxhash via pkg-config
XXHASH_CFLAGS=$(pkg-config --cflags libxxhash)
XXHASH_LDFLAGS=$(pkg-config --libs-only-L libxxhash) $(pkg-config --libs-only-L libxxhash | sed 's/-L/-Wl,-rpath/g')
XXHASH_LDLIBS=$(pkg-config --libs libxxhash)

mkdir -p configs
confFile=configs/${NAME}.config
tmpFile=${confFile}.tmp
cat > "$tmpFile" <<EOF
CONFIG_ALPHA=$ALPHA
CONFIG_K=$K
CONFIG_CXX=$GCXX
CONFIG_CXXFLAGS=$XXHASH_CFLAGS $OPTFLAGS $CXXFLAGS
CONFIG_LDFLAGS=$XXHASH_LDFLAGS $LDFLAGS
CONFIG_LDLIBS=$XXHASH_LDLIBS $LDLIBS
CONFIG_YAGGO=$YAGGO
CONFIG_COMPILEDB=$COMPILEDB
EOF

[ -n "$ILPPYTHON" ] && echo "CONFIG_ILP_PYTHON=${ILPPYTHON}" >> "$tmpFile"

# Testing
if [[ ${#FILES[@]} -gt 0 ]]; then
cat >> "$tmpFile" <<EOF
CONFIG_EXP_REPEAT=$REPEAT
CONFIG_EXP_FILES=${FILES[@]}
CONFIG_EXP_HISTO_THRESH=50
EOF
fi

if [ -n "$DRYRUN" ]; then
echo "Dryrun. Config: ${confFile} ${NAME}"
cat "$tmpFile"
rm "$tmpFile"
exit 0
fi

mv -f "$tmpFile" "$confFile"
[ -d ".tup" ] || { echo "Initialize tup"; tup init; }
[ -d "build-${NAME}" ] || "$TUP" variant "$confFile"
