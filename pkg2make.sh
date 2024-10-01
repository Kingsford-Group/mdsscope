#!/usr/bin/env bash

# Bad one to parse a .pc file and create Makefile variables

set -e

FILE=$1
echo >&2 "Processing $FILE"

source <(grep -iv '^[a-z]*:' "$FILE")
cflags=$(grep '^Cflags:' "$FILE" | sed 's/Cflags: //')
libs=$(grep '^Libs:' "$FILE" | sed 's/^Libs: //')
libsl=$(eval "echo $libs" | sed 's/-L[^ ]*//g')
libsL=$(eval "echo $libs" | sed 's/-l[^ ]*//g')
echo "CXXFLAGS += " $(eval "echo $cflags")
echo "LDFLAGS += " "$libsL" $(echo $libsL | sed 's/-L/-Wl,-rpath,/')
echo "LDLIBS += " "$libsl"
