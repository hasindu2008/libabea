#!/bin/bash

make clean
make
make pylib

echo "C test"
./abea_example > out.txt || echo "execution fail"
diff test/stdout.exp out.txt || echo "C test failed"

echo "Python test"
python3 example.py | tr -d '[],' | tr ' ' '\t' >  outpy.txt || echo "execution fail"
diff test/stdout.exp outpy.txt || echo "Python test failed"
