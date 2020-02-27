#!/bin/bash

python3 setup.py clean
make clean
make
python3 setup.py build && cp build/lib.linux-x86_64-3.5/test/abea.cpython-35m-x86_64-linux-gnu.so ./ && python3 test.py



