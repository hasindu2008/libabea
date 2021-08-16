#!/bin/bash


make clean
GCC_ASAN_PRELOAD=$(gcc -print-file-name=libasan.so)
CFLAGS="-fsanitize=address -fno-omit-frame-pointer" python3 setup.py build
cp build/lib.*/*.so  ./
echo $GCC_ASAN_PRELOAD
LD_PRELOAD=$GCC_ASAN_PRELOAD  python3 example.py > /dev/null
