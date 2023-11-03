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

Options:
  -d          Debugging exec
  -h          Help
EOF
}

OPTFLAGS="-O3 -DNDEBUG"
SUFFIX=

while getopts "dh" o; do
    case "${o}" in
        (d) OPTFLAGS="-O0 -g"; SUFFIX="-debug" ;;
        (h) help; exit 0 ;;
        (*) usage; exit 1 ;;
    esac
done
shift $((OPTIND-1))

ALPHA=$1
K=$2

if [ -z "$ALPHA" ] || [ -z "$K" ]; then
    usage
    exit 1
fi

NAME="A${ALPHA}K${K}${SUFFIX}"

[ -z "$YAGGO" ] && YAGGO=$(which yaggo 2>/dev/null)

mkdir -p configs
cat > "configs/${NAME}.config" <<EOF
CONFIG_ALPHA=$ALPHA
CONFIG_K=$K
CONFIG_CXXFLAGS=$OPTFLAGS $CXXFLAGS
CONFIG_LDFLAGS=$LDFLAGS
CONFIG_LDLIBS=$LDLIBS
CONFIG_YAGGO=$YAGGO
EOF

tup variant "configs/${NAME}.config"
