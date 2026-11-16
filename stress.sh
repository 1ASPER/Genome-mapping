#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./stress.sh <problem_letter>"
    exit 1
fi

PROB=$1
MAIN="${PROB}.cpp"
BRUTE="${PROB}_BruteForce.cpp"
GEN="${PROB}_TestGen.cpp"

if [ ! -f "$BRUTE" ] || [ ! -f "$GEN" ]; then
    if [ ! -f "$BRUTE" ]; then
        cp "$MAIN" "$BRUTE"
        echo "Created $BRUTE"
    fi

    if [ ! -f "$GEN" ]; then
        cat > "$GEN" << 'GEN_EOF'
#include <bits/stdc++.h>

using namespace std;

const int MAX_T = 10;
const int MAX_N = 10;
const int MAX_A = 100;

mt19937_64 rng;

int rnd(int L, int R) {
    return uniform_int_distribution<int>(L, R)(rng);
}

void solve() {
    int n = rnd(1, MAX_N);
    cout << n << "\n";
    for (int i = 0; i < n; i++) {
        cout << rnd(1, MAX_A) << (i == n - 1? "": " ");
    }
    cout << "\n";
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    rng.seed(argc > 1? atoi(argv[1]): 1337);

    int t = rnd(1, MAX_T);
    cout << t << "\n";
    while (t--) {
        solve();
    }
    return 0;
}
GEN_EOF
        echo "Created $GEN"
    fi

    echo "Edit $BRUTE and $GEN, then run ./stress.sh $PROB again to start testing."
    exit 0
fi

for v in 18 17 16 15 14; do command -v g++-$v >/dev/null && { cxx=g++-$v; break; }; done
[ -n "$cxx" ] || { echo "Error: GCC not found (brew install gcc)"; exit 1; }

flags="-std=c++17 -O2 -Wall -Wextra -Wshadow -D_GLIBCXX_DEBUG"

echo "Compiling..."
mkdir -p outs
$cxx $flags "$MAIN" -o outs/main_out || { echo "Compilation failed: $MAIN"; exit 1; }
$cxx $flags "$BRUTE" -o outs/brute_out || { echo "Compilation failed: $BRUTE"; exit 1; }
$cxx $flags "$GEN" -o outs/gen_out || { echo "Compilation failed: $GEN"; exit 1; }

echo "Starting stress test for $PROB..."
for ((i=1; ; ++i)); do
    outs/gen_out $i > outs/test.in
    outs/main_out < outs/test.in > outs/main.out
    outs/brute_out < outs/test.in > outs/brute.out

    if diff -w outs/main.out outs/brute.out > /dev/null; then
        echo -ne "Test $i: OK\r"
    else
        echo -e "\nTest $i: FAILED! (seed $i)"
        echo "================ INPUT ================"
        cat outs/test.in
        echo "============= MAIN OUTPUT ============="
        cat outs/main.out
        echo "============= BRUTE OUTPUT ============"
        cat outs/brute.out
        break
    fi
done
