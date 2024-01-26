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

[ -z "$YAGGO" ] && YAGGO=$(which yaggo 2>/dev/null)
[ -z "$TUP" ] && TUP=$(which tup 2>/dev/null)

[[ -z "$TUP" || -z "$YAGGO" ]] && { echo "Missing required dependencies: tup and/or yaggo"; false; }

detect_compiledb() {
  v=$("$TUP" --version | sed -e 's/^tup v\?//' -e 's/-.*$//')
  a=( ${v//./ } )
  if [ "${a[1]}" -gt "7" ] || [[ "${a[1]}" -eq "7" && "${a[2]}" -ge "11" ]]; then
  echo yes
  fi
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

mkdir -p configs
confFile=configs/${NAME}.config
tmpFile=${confFile}.tmp
cat > "$tmpFile" <<EOF
CONFIG_ALPHA=$ALPHA
CONFIG_K=$K
CONFIG_CXXFLAGS=$OPTFLAGS $CXXFLAGS
CONFIG_LDFLAGS=$LDFLAGS
CONFIG_LDLIBS=$LDLIBS
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
[ -d ".tup" ] || { echo "Initialize tup"; tup init }
[ -d "build-${NAME}" ] || "$TUP" variant "$confFile"
