from __future__ import print_function
import sys
from .htsffi import libhts, ffi, _raise_if_null


class VCF(object):
    r"""
    >>> vcf = VCF('/usr/local/src/gocode/src/github.com/brentp/vcfgo/examples/test.query.vcf')
    >>> vcf #doctest: +ELLIPSIS
    VCF('...')

    >>> v = next(vcf)
    >>> v
    Variant('chr1:30547')

    >>> v.genotypes[:10]
    [2, 2, 2, 2, 2, 2, 0, 2, 2, 2]

    """

    def __init__(self, fname, mode="r"):
        htf = self._htf = libhts.hts_open(fname, mode)
        hdr = self._hdr = libhts.bcf_hdr_read(htf)
        assert libhts.bcf_hdr_set_samples(hdr, "-", 0) == 0, ("error setting samples")
        self.fname = fname

    def __repr__(self):
        return "%s('%s')" % (self.__class__.__name__, self.fname)

    def __iter__(self):
        return self

    def __next__(self):
        bcf_t = libhts.bcf_init()
        if libhts.bcf_read(self._htf, self._hdr, bcf_t) >= 0:
            return Variant(bcf_t, self)
        raise StopIteration

    next = __next__

    def seq(self, tid):
        s = libhts.bcf_hdr_id2name(self._hdr, tid)
        return ffi.string(s)

    @property
    def n_samples(self):
        return self._hdr.n[libhts.BCF_DT_SAMPLE]

class Variant(object):
    """Variants are created from VCF"""

    def __init__(self, bcf_t, vcf):
        self._bcf = bcf_t
        libhts.bcf_unpack(bcf_t, 15)  # unpack everything.

        self.vcf = vcf
        self.chrom = vcf.seq(bcf_t.rid)
        self.pos = bcf_t.pos

    def __repr__(self):
        return "%s('%s:%d')" % (self.__class__.__name__, self.chrom, self.pos)

    @property
    def formats(self):
        return Format(self._bcf.d.fmt)

    @property
    def genotypes(self):

        dst = ffi.new("int **")
        ndst = ffi.new("int *")
        n_gts = libhts.bcf_get_genotypes(self.vcf._hdr, self._bcf, dst, ndst)

        if n_gts < 0:
            raise Exception("cant pull genotypes")

        n_samples = self.vcf.n_samples
        gts_per_sample = n_gts / n_samples
        if gts_per_sample != 2:
            return [None] * n_samples

        n = libhts.as_gts(dst[0], n_samples)
        assert n == n_samples, ("incorrect number of genotypes, got", n, "expected", n_samples)
        try:
            return [dst[0][i] for i in range(n_samples)]
        finally:
            libhts.free(dst[0])

    gt_types = genotypes

class Format(object):
    def __init__(self, fmt_t):
        self._fmt_t = fmt_t



