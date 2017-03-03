.. hts-python documentation master file, created by
   sphinx-quickstart on Wed Mar  1 09:40:03 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

hts-python
==========

hts-python is a pythonic wrapper for `htslib <https://github.com/samtools/htslib/>`_ using python `cffi <https://cffi.readthedocs.io/en/latest/>`_.

A taste of hts-python...

.. code-block:: python

    >>> import os.path as op

    >>> from hts import Bam
    >>> bam = Bam("hts/test/small.bam") 
    >>> list(bam.header.seqs)
    ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']

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
    >>> a.qual[:10]
    [56, 63, 53, 62, 64, 62, 51, 44, 58, 59]

    >>> a.flag, a.flag_str
    (16, 'REVERSE')

    >>> a.cigar
    Cigar('36M')

    >>> str(a)[:40]
    'HWUSI-NAME:2:69:512:1017#0\t16\tchr2L\t9330'

See Also
========

.. toctree::
   :maxdepth: 2

