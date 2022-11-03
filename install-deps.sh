#!/bin/bash

set -euo pipefail

rm -rf CombBLAS
git submodule init
git submodule update
git submodule foreach --recursive git reset --hard

cd seqan
git fetch
git checkout develop

cd ..

cd CombBLAS
mkdir build
mkdir install
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../ -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX
make -j$(nproc)
make install 