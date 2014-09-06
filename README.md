hts-python
==========

pythonic wrapper for [htslib](https://github.com/samtools/htslib.git) C-API using python [cffi](https://cffi.readthedocs.org/).

There is enough functionality for this to be useful, but it still needs a lot of work.

[![Build Status](https://travis-ci.org/brentp/hts-python.svg?branch=v0.4.1)](https://travis-ci.org/brentp/hts-python)


A taste
-------

```Python
>>> import os.path as op

>>> from hts import Bam
>>> bam = Bam("hts/test/small.bam") #bam stolen from pybedtools [thanks]
>>> list(bam.header.seqs)
['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']

# region query creates index if needed:
>>> a = next(bam('chr2L:9000-11000'))
>>> a
Alignment('HWUSI-NAME:2:69:512:1017#0')
>>> a.target, a.pos, a.strand
('chr2L', 9329, '-')
>>> a.qlen, a.rlen
(36, 36)
>>> a.strand
'-'
>>> a.seq
'TACAAATCTTACGTAAACACTCCAAGCATGAATTCG'
>>> a.base_q[:10]
[56, 63, 53, 62, 64, 62, 51, 44, 58, 59]

>>> a.flag, a.flag_str
(16, 'REVERSE')

>>> a.cigar
Cigar('36M')

>>> str(a)[:40]
'HWUSI-NAME:2:69:512:1017#0\t16\tchr2L\t9330'

```

There are also wrappers for:
+ Fai for fasta querying fasta files.
+ Tbx for tabix files (indexed bed/gff/sam, etc.).
+ fisher for fisher's exact test.

Installation
------------

1. Install [htslib](https://github.com/samtools/htslib.git htslib) using `make install`
2. pip/easy_install python [cffi](https://cffi.readthedocs.org/).
3. run `python setup.py install` (--user) from this directory.

Development
-----------

This is a work in progress that relies on the hts library. All of the wrapped functions are included in `hts/hts_concat.h` and then available from python as, e.g. `htslib.sam_read1`

When C-functions not provided by the api are needed, they are added to `hts_extra.c/.h`.

One can run the tests with: `python -c "import hts; hts.doctests()"`

Things to work on:

1. Make properties settable in hts.bam. 
   Currently, they are read-only properties. At very least, it will be useful
   to have setters for seq, base_q, qname, tname, pos, strand, flag.

2. Wrap B/VCF stuff?


Why
---

Why use this when `pysam` exists? It's an experiment with python cffi and to provide a pythonic access to htslib.
