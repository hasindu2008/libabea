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

# How to use:

Experimental use on RAGE-seq single-cell data

```sh
/usr/bin/time -v python3 -m cProfile -s cumtime motif2sig.py -f demuxed.fastq -s sequencing_summary.txt -p fast5s/ -q ../test_squig.tsv -v 1 -t 0 > test.out.tsv
```

# Usage

See [example.py](example.py)

## Acknowledgement
This reuses code and methods from [Nanopolish](https://github.com/jts/nanopolish).
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
