import os.path as op
from .htsffi import libhts, ffi, _raise_if_null
import atexit

class Tbx(object):
    """
    >>> tbx = Tbx('%s/test/example.gtf.gz' % op.dirname(__file__))

    >>> tbx.sequences
    [u'chr1', u'chr2']

    >>> for s in tbx('chr1:1-1800'):
    ...     print(s[:5])
    ['chr1', 'ENSEMBL', 'UTR', 1737, 2090]
    ['chr1', 'ENSEMBL', 'exon', 1737, 2090]
    ['chr1', 'ENSEMBL', 'transcript', 1737, 4275]
    ['chr1', 'HAVANA', 'gene', 1737, 4275]
    """
    def __init__(self, fname):
        assert op.exists(fname), ("Tbx: no file", fname)
        if not op.exists("%s.tbi" % fname):
            if fname.endswith('.bed.gz'):
                Tbx.build(fname)
            elif fname.endswith(('.gff.gz', '.gtf.gz')):
                Tbx.build(fname, 1, 4, 5, '#', 0)
            elif fname.endswith(('.vcf.gz')):
                Tbx.build(fname, 1, 2, 0, '#', 0)
            else:
                raise Exception('%s.tbi not found and filetype not known' % fname)

        tbx = self._tbx = libhts.tbx_index_load(fname);
        _raise_if_null(tbx, "Tbx:unable to find %s.tbi" % fname)

        htf = self._htf = libhts.hts_open(fname, "r");
        _raise_if_null(tbx, "Tbx:unable to find %s" % fname)

        atexit.register(libhts.hts_close, htf)
        atexit.register(libhts.tbx_destroy, tbx)

    @classmethod
    def build(cls, fname, seq_col=1, start_col=2, end_col=3, comment="#",
                 line_skip=0):
        """
        1-based columns
        """
        conf = ffi.new('tbx_conf_t *')
        conf.sc, conf.bc, conf.ec = seq_col, start_col, end_col
        conf.meta_char = ffi.cast('char', comment)
        conf.preset = 0
        conf.line_skip = line_skip
        assert libhts.tbx_index_build(fname, -1, conf) != -1
        return Tbx(fname)


    @property
    def sequences(self):
        n = ffi.new("int *")
        cnames = libhts.tbx_seqnames(self._tbx, n)
        names = [ffi.string(cnames[i]).decode() for i in range(n[0])]
        try:
            return names
        finally:
            libhts.free(cnames)

    def __call__(self, region):
        itr = libhts.tbx_itr_querys(self._tbx, region)
        s = ffi.new("kstring_t *")
        conf = self._tbx.conf

        # get 0-based cols so we can convert to ints
        start_col, end_col = conf.bc - 1, conf.ec - 1

        try:
            slen = libhts.tbx_itr_next(self._htf, self._tbx, itr, s)
            while slen > 0:
                toks = ffi.string(s.s, slen).split("\t")

                # convert the start and end columns to int
                toks[start_col] = int(toks[start_col])
                toks[end_col] = int(toks[end_col])
                yield toks

                slen = libhts.tbx_itr_next(self._htf, self._tbx, itr, s)
        finally:
            libhts.free(s.s)
            libhts.hts_itr_destroy(itr)

    query = __call__

