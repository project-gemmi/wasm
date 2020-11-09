#!/bin/bash -eu

em++ -Os -Wall -Wextra -std=c++20 -I../../gemmi/include \
    -s WASM=1 \
    -s STRICT=1 \
    -s DISABLE_EXCEPTION_CATCHING=0 \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s 'EXPORTED_FUNCTIONS=["_pdb2cif", "_mtz2cif", "_mtzpair2cif", "_clear_string", "_get_str2", "_get_version"]' \
    -s 'EXTRA_EXPORTED_RUNTIME_METHODS=["ccall", "writeArrayToMemory"]' \
    pdb2cif.cpp -o pdb2cif.js \
    #--emrun
