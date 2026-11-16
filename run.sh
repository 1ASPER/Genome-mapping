#!/bin/bash
# Usage: ./run.sh A [-d]

[ $# -ge 1 ] || { echo "Usage: ./run.sh <problem> [-d]"; exit 1; }
[ -f "$1.cpp" ] || { echo "Error: $1.cpp not found"; exit 1; }

for v in 18 17 16 15 14; do command -v g++-$v >/dev/null && { cxx=g++-$v; break; }; done
[ -n "$cxx" ] || { echo "Error: GCC not found (brew install gcc)"; exit 1; }

flags="-std=c++17 -DLOCAL -Wall -Wextra -Wshadow"
[ "$2" = "-d" ] && flags+=" -O0 -g -D_GLIBCXX_DEBUG" || flags+=" -O2"

mkdir -p outs
$cxx $flags "$1.cpp" -o "outs/$1.out" || exit 1
ulimit -s 65532 2>/dev/null

if [ -f "$1.in" ]; then "./outs/$1.out" < "$1.in"; else "./outs/$1.out"; fi
code=$?; [ $code -eq 0 ] || echo "Exited with code $code"
