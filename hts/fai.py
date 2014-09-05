from .htsffi import libhts, ffi
import os.path as op
import atexit

class Fai(object):

    """
    Random access to fasta files.

    uses 1-based starts for regions

    Arguments
    ---------

    fn: str
        fasta file path.

    Examples
    --------

    >>> f = Fai('%s/test/t.fa' % op.dirname(__file__))
    >>> f('chr1:1-10')
    'AGAAAACCCCCCCAC'
    >>> f.nseqs
    1
    >>> 'chr1' in f
    True
    >>> 'chr2' in f
    False

    >>> for chrom, clen in f:
    ...     print chrom, clen
    chr1 46

    >>> dict(f)
    {'chr1': 46}
    """

    def __init__(self, fn):
        """Create a new Fai object

        Examples
        --------

        >>> f = Fai('%s/test/t.fa' % op.dirname(__file__))
        """
        if fn.endswith(".fai"):
            fn = fn[:-4]
        assert op.exists(fn), (fn, "does not exist")
        if not op.exists("%s.fai"):
            assert libhts.fai_build(fn) == 0
        self._fai = libhts.fai_load(fn)
        atexit.register(libhts.fai_destroy, self._fai)
        self.fn = fn

    def __call__(self, region):
        """Extract a region.

        Arguments
        ---------

        region: str
        """
        rlen = ffi.new("int *")
        seq = libhts.fai_fetch(self._fai, 'chr1:1-15', rlen)
        if rlen == -2:
            raise Exception("sequence not present")
        elif rlen == -1:
            raise Exception("error fetching sequence")
        return ffi.string(seq)

    @property
    def nseqs(self):
        """Count number of sequences in the fastq."""
        return libhts.faidx_fetch_nseq(self._fai)

    def __contains__(self, seq):
        """Test if a sequence name is in the fasta."""
        return bool(libhts.faidx_has_seq(self._fai, seq))

    def __iter__(self):
        """Generate names and lengths of the sequences."""
        with open(self.fn + ".fai") as fai:
            for toks in (l.rstrip("\r\n").split("\t") for l in fai):
                yield toks[0], int(toks[1])
