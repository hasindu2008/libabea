# Adaptive Banded Event Alignment (abea)

Python library bindings for Adaptive Banded Event Alignment (ABEA) for Nanopore reads. Used for aligning a base-called read to the raw nanopore signal. 

## Building

```sh
pip3 install cython
python3 setup.py build 
cp build/lib.linux-x86_64-3.5/test/abea.cpython-35m-x86_64-linux-gnu.so ./
```

# Launch a test run

```sh           
./test.sh
```

# Usage

See [example.py](example.py)

## Acknowledgement 
This reuses code and methods from [Nanopolish](https://github.com/jts/nanopolish).
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).

