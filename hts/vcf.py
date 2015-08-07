from .htsffi import libhts, ffi, _raise_if_null


class VCF(object):
    r"""
    >>> vcf = VCF('/usr/local/src/gocode/src/github.com/brentp/vcfgo/examples/test.query.vcf')
    >>> vcf #doctest: +ELLIPSIS
    VCF('...')

    >>> v = next(vcf)
    >>> v
    Variant('chr1:30547')

    >>> v.formats

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
        assert libhts.bcf_read(self._htf, self._hdr, bcf_t) == 0, ("error reading")
        return Variant(bcf_t, self)

    next = __next__

    def seq(self, tid):
        s = libhts.bcf_hdr_id2name(self._hdr, tid)
        return ffi.string(s)

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

class Format(object):
    def __init__(self, fmt_t):
        self._fmt_t = fmt_t



