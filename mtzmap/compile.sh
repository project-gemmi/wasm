#!/bin/bash -eu

em++ --bind -O3 -Wall -Wextra -std=c++14 -I../../gemmi/include \
    -s WASM=1 \
    -s DISABLE_EXCEPTION_CATCHING=0 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s 'EXTRA_EXPORTED_RUNTIME_METHODS=["writeArrayToMemory"]' \
    mtzmap.cpp -o mtzmap.js \
    #--emrun
