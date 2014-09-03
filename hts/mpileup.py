from .htsffi import libhts, ffi
from .bam import Bam, Alignment

@ffi.callback("int(void *, bam1_t*)")
def _callback(data, b):
    aln = Alignment(b, ffi.cast("bam_hdr_t", data))
    print aln
    return




class Pile(object):
    def __init__(self, bams, merge_overlaps=True, max_count=None):
        n = len(bams)
        self._data = ffi.new("void **data")
        self._iter = libhts.bam_mplp_init(n, _callback, self._data)
        if merge_overlaps:
            libhts.bam_mplp_init_overlaps(self._iter)
        if max_counts:
            libhts.bam_mplp_set_maxcnt(self._iter, max_count)

    def __iter__(self):
        libhts.bam_plp_reset(self._iter)
        return self

    def next(self):
        libhts.bam_plp_reset(self._iter)
        pass



    __next__ = next

