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
  -h          Help
EOF
}

OPTFLAGS="-O3 -DNDEBUG"
SUFFIX=
COMPILEDB=
REPEAT=1
FILES=()

while getopts "djr:f:h" o; do
    case "${o}" in
        (d) OPTFLAGS="-O0 -g"; SUFFIX="-debug" ;;
        (h) help; exit 0 ;;
        (j) COMPILEDB=yes ;;
        (r) REPEAT=$OPTARG ;;
        (f) FILES+=($OPTARG) ; SUFFIX="-exp";;
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

mkdir -p configs
confFile=configs/${NAME}.config
cat > "$confFile" <<EOF
CONFIG_ALPHA=$ALPHA
CONFIG_K=$K
CONFIG_CXXFLAGS=$OPTFLAGS $CXXFLAGS
CONFIG_LDFLAGS=$LDFLAGS
CONFIG_LDLIBS=$LDLIBS
CONFIG_YAGGO=$YAGGO
CONFIG_COMPILEDB=$COMPILEDB
EOF

# Testing
if [[ ${#FILES[@]} -gt 0 ]]; then
cat >> "$confFile" <<EOF
CONFIG_EXP_REPEAT=$REPEAT
CONFIG_EXP_FILES=${FILES[@]}
EOF
fi


"$TUP" variant "configs/${NAME}.config"
