#!/bin/bash

set -euo pipefail

rm -rf build install
mkdir -p build install

cd build
# cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../ -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../ -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX

#cmake .. -DCMAKE_INSTALL_PREFIX=../instal
