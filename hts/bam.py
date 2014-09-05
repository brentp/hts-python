from __future__ import print_function
from .htsffi import libhts, ffi, _raise_if_null
import atexit
import os.path as op
import sys

class BamHeader(object):
    def __init__(self, h):
        self._h = h

    @property
    def seqs(self):
        for i in range(self._h.n_targets):
            yield ffi.string(self._h.target_name[i])

class Cigar(object):

    """
    A cigar object usually created from `Alignment`.

    Takes the uin32_t * from bam_get_cigar()

    """

    ops = "MIDNSHP=XB"

    # cigar is the uint32_t * from bam_get_cigar()
    def __init__(self, cigar, n_cigar):
        self._c = cigar
        self.n_cigar = n_cigar
    def __str__(self):
        """Str."""
        cig, s = self._c, []
        for k in range(self.n_cigar):
            s.extend((str(libhts.bam_cigar_oplen(cig[k])),
                      libhts.bam_cigar_opchr(cig[k])))
        return "".join(s)

    def __repr__(self):
        """repr."""
        return "Cigar('%s')" % str(self)

class Alignment(object):

    """Alignment object; usually created by iterating over a `Bam` object."""

    def __init__(self, bam1_t, bam_hdr_t):
        self._b = bam1_t
        self._h = bam_hdr_t

    def copy(self):
        """Create a copy of the alignment for modification or storage."""
        b = libhts.bam_dup1(self._b)
        atexit.register(libhts.bam_destroy1, b)
        return Alignment(b, self._h)

    @property
    def aux(self):
        """The auxillary tags from the alignment."""
        auxs = []
        l = libhts.bam_get_l_aux(self._b)
        if not l : return None

        s_aux_ptr = libhts.bam_get_aux(self._b)
        i = 0
        #print(str(self), file=sys.stderr)
        while i + 4 < l:
            key = "".join(chr(s) for s in s_aux_ptr[i:i + 2])
            ftype = chr(s_aux_ptr[i + 2])
            if ftype == 'Z':
                j = 0
                while chr(s_aux_ptr[i + 4 + j]) != '\0':
                    j += 1
                val = s_aux_ptr[i + 3: i + 4 + j]
                val = "".join("%c" % val[i] for i in range(1 + j))
                i += j + 1
            elif ftype in "iI":
                t = "uint32_t" if ftype == "I" else "int32_t"
                val = ffi.cast("%s *" % t, s_aux_ptr + i + 3)[0]
            else:
                val = s_aux_ptr[i + 3: i + 4][0]

            #print(l, key, ftype, val, file=sys.stderr)
            auxs.append((key, ftype, val))
            i += 4
        # TODO: other types from htslib/sam.c
        return auxs

    def adjust_overlap_quality(self, other):
        """
        Adjust base-quality of overlapping reads.

        For paired-end reads, we may not want to double-count both ends.
        This function will set (in-place), the base-quality of `other` to
        0 and self to other + self for any overlapping bases.
        If bases are not equal then quality of self is set to 80% of its
        original value.
        """
        libhts.tweak_overlap_quality(self._b, other._b)

    @property
    def qname(self):
        """Read (query) name."""
        return ffi.string(libhts.bam_get_qname(self._b))

    @qname.setter
    def qname(self, name):
        """Set (query) name."""
        raise NotImplemented
        # need to check if proposed name is longer/shorter
        # than current and calloc self._b.data based on that!!!
        current_name = libhts.bam_get_qname(self._b)

    @classmethod
    def from_sam_str(cls, sam_str, bam_hdr):
        s = ffi.new('kstring_t *', {'m': 0, 'l': 0, 's': ffi.NULL})
        libhts.kputsn(sam_str, len(sam_str) + 1, s)

        b = libhts.bam_init1()
        res = libhts.sam_parse1(s, bam_hdr, b)
        assert res <= 0, ("SAM parse error", res)
        h2 = libhts.bam_hdr_dup(bam_hdr)
        return Alignment(b, h2)

    @property
    def tname(self):
        """Chromosome or sequence aligned to."""
        tid = self._b.core.tid
        return ffi.string(self._h.target_name[tid])
    target = tname

    @property
    def strand(self):
        """Strand of alignment."""
        return "-" if libhts.bam_is_rev(self._b) else "+"

    @property
    def base_qualities(self):
        """Base-qualities of query sequence."""
        q_ptr = libhts.bam_get_qual(self._b)
        has_qual = q_ptr[0] != 0xff
        if not has_qual: return None
        return [q_ptr[i] for i in range(self._b.core.l_qseq)]

    base_q = base_qualities

    @property
    def pos(self):
        """Left-most position of alignment."""
        return self._b.core.pos

    @pos.setter
    def pos(self, value):
        """Set the position."""
        self._b.core.pos = value

    @property
    def isize(self):
        """Insert size of alignment if applicable."""
        return self._b.core.isize

    @property
    def mapping_quality(self):
        """Mapping quality of alignment."""
        return self._b.core.qual

    @mapping_quality.setter
    def mapping_quality(self, value):
        """Set the mapping quality to a new value."""
        self._b.core.qual = value

    map_q = mapping_quality

    @property
    def cigar(self):
        """Cigar object of alignment."""
        return Cigar(libhts.bam_get_cigar(self._b), self._b.core.n_cigar)

    @property
    def qlen(self):
        """Length of alignment by query."""
        return libhts.bam_cigar2qlen(self._b.core.n_cigar,
                                  libhts.bam_get_cigar(self._b))

    @property
    def rlen(self):
        """Length of alignment on reference."""
        return libhts.bam_cigar2rlen(self._b.core.n_cigar,
                                  libhts.bam_get_cigar(self._b))

    @property
    def seq(self):
        """Nucleotide sequence of read."""
        kstr = ffi.new('kstring_t *', {'m': 0, 'l': 0, 's': ffi.NULL})
        r = libhts.bam_get_read_seq(self._b, kstr)
        assert r == kstr.l, (r, kstr.l)
        return ffi.string(kstr.s)

    @property
    def flag_str(self):
        """Alignment flag as a string."""
        return ffi.string(libhts.bam_flag2str(self.flag))

    @property
    def flag(self):
        """Alignment flag."""
        return self._b.core.flag

    def __eq__(self, other):
        """test equality."""
        if self.qname != other.qname: return False
        if (self.pos, self.seq, self.map_q, self.qname) == \
                (other.pos, other.seq, other.map_q, other.qname):
            return True
        return False

    def __repr__(self):
        """repr of aln."""
        return "%s('%s')" % (self.__class__.__name__, self.qname)

    def __str__(self):
        """str of aln."""
        s = ffi.new("kstring_t *")
        l = libhts.sam_format1(self._h, self._b, s)
        return ffi.string(s.s, l)

class Bam(object):

    r"""
    A BAM class with properties that call the C-API.

    Parameters
    ----------

    fname : string
        file name of bam

    mode : string
        mode for opening file r/w

    create_index: bool, string, optional
        if True, then force create index. If "auto", then only create index
        if it doesn't exist

    header: bam_hdr_t * or str
        ffi bam_hdr_t object for use when mode == "w"
        or a SAM header string.

    Examples
    --------

    >>> import os.path as op
    >>> bam = Bam("%s/test/small.bam" % op.dirname(__file__))
    >>> next(bam)
    Alignment('HWUSI-NAME:2:69:512:1017#0')
    >>> region_iter = bam('chr2L:9000-11000')
    >>> a = next(region_iter)

    # save copy for later test.
    >>> asav = a.copy()

    >>> list(bam.header.seqs)
    ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']


    Query a specific Region
    >>> a
    Alignment('HWUSI-NAME:2:69:512:1017#0')
    >>> a.qname
    'HWUSI-NAME:2:69:512:1017#0'
    >>> a.tname # or a.target
    'chr2L'

    >>> a.cigar
    Cigar('36M')
    >>> a.qlen, a.rlen
    (36, 36)

    >>> a.strand
    '-'

    >>> a.seq
    'TACAAATCTTACGTAAACACTCCAAGCATGAATTCG'
    >>> b = next(region_iter)

    >>> b.seq
    'TGTAGAATGCAAAAATTACATTTGTGAGTATCATCA'

    # NOTE that a.seq has changed. this is because we updated the underlying
    # pointer and htslib assumes streaming stuff.
    >>> a.seq
    'TGTAGAATGCAAAAATTACATTTGTGAGTATCATCA'

    # We can use the copy to avoid these problems:
    >>> asav.seq
    'TACAAATCTTACGTAAACACTCCAAGCATGAATTCG'
    >>> a = asav

    # in the future, we may default to copy semantics or allow that to be specified
    # when creating a BAM.

    >>> a.flag, a.flag_str
    (16, 'REVERSE')
    >>> a.base_qualities[:10]
    [56, 63, 53, 62, 64, 62, 51, 44, 58, 59]

    >>> a.mapping_quality # or a.map_q
    3
    >>> a.pos
    9329

    >>> a.isize
    0

    >>> str(a)
    'HWUSI-NAME:2:69:512:1017#0\t16\tchr2L\t9330\t3\t36M\t*\t0\t0\tTACAAATCTTACGTAAACACTCCAAGCATGAATTCG\tY`V_a_TM[\\_V`abb`^^Q]QZaaaaa_aaaaaaa\tNM:i:0\tNH:i:2\tCC:Z:chrX\tCP:i:19096815'

    >>> a.aux
    [('NM', 'C', 0), ('NH', 'C', 2), ('CC', 'Z', 'chrX'), ('CP', 'I', 19096815)]

    >>> a.base_q[:10]
    [56, 63, 53, 62, 64, 62, 51, 44, 58, 59]

    # for paired end reads, we may want to avoid double-counting
    # the overlapping regions. We can use the samtools algorithm
    # to give a, the base-quality of a + b when the sequence is
    # the same or 0.8 * a when the seq is different. In either case,
    # the base-quality of b is set to 0.
    >>> b = a.copy()

    >>> a.adjust_overlap_quality(b)
    >>> a.base_q[:10]
    [112, 126, 106, 124, 128, 124, 102, 88, 116, 118]

    >>> b.base_q[:10]
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    >>> all(bq == 0 for bq in b.base_q)
    True

Writing.

    >>> out = Bam('out.bam', 'wb', header=a._h)
    >>> c = next(region_iter)
    >>> out.write(c, a, b); out.close()

    >>> cnew = next(Bam('out.bam'))
    >>> cnew.seq == c.seq
    True

    >>> (cnew == c, b.qname == a.qname)
    (True, True)

    we can also create a new writable bam with a fasta to create the bam header:
    >>> outf = Bam('outf.bam', mode='wb', fasta="hts/test/t2.fa")

    >>> outf.write(a); outf.close()
    >>> b = Bam('outf.bam')
    >>> b
    Bam('outf.bam', 'r')

    >>> next(b) == a
    True


    >>> a.pos
    9329
    >>> a.pos += 1
    >>> a.pos
    9330

    >>> assert a.map_q == 3
    >>> a.map_q += 55
    >>> assert a.mapping_quality == 58

    # string to object back to string.
    >>> s = str(a)
    >>> c = Alignment.from_sam_str(s, a._h)
    >>> str(c) == s
    True

    #>>> a.qname == "HELLO"
    #>>> assert "HELLO" in str(a)
    """

    def __init__(self, fname, mode="r", create_index="auto", header=None, fasta=None):
        self.fn, self.mode = fname, mode

        htf = self._htf = libhts.hts_open(fname, mode)
        _raise_if_null(htf, "Bam: bad file %s" % fname)

        if mode[0] == "r":

            idx = self._idx = libhts.sam_index_load(self._htf, fname)
            if (idx == ffi.NULL and create_index == "auto") or create_index is True:
                libhts.bam_index_build(fname, -1)
                idx = self._idx = libhts.sam_index_load(self._htf, fname)
                if idx == ffi.NULL:
                    # no querying.
                    print("NO querying", file=sys.stderr)
                    pass

            self.header = BamHeader(libhts.sam_hdr_read(htf))
            b = self._b = libhts.bam_init1()
            atexit.register(libhts.bam_destroy1, b)

        else:
            assert mode[0] == "w" and (header or fasta)
            if hasattr(header, "_h"):
                header = header._h
            elif fasta:
                assert op.exists(fasta) # it's a fasta file.
                header_str = Bam.header_from_fasta(fasta)
                header = libhts.sam_hdr_parse(len(header_str), header_str)
                header.l_text = len(header_str)
                header.text = ffi.new("char[]", header_str)

            self.header = header
            libhts.sam_hdr_write(self._htf, self.header);

    def __repr__(self):
        return "Bam('%s', '%s')" % (self.fn, self.mode)

    def write(self, *alns):
        """
        Write alignments to the current file.

        Parameters
        ----------
        *alns: bam1_t *

        """
        for aln in alns:
            assert libhts.sam_write1(self._htf, self.header, aln._b) > 0, "error"

    @classmethod
    def header_from_fasta(self, fasta, sort_order='unknown'):
        """Create a sam header from a fasta file.
        
        >>> import os.path as op
        >>> fa = '%s/test/t.fa' % op.dirname(__file__)
        >>> print(Bam.header_from_fasta(fa)) #doctest: +NORMALIZE_WHITESPACE
        @HD VN:1.0 SO:unknown
        @SQ SN:chr1 LN:46
        """
        from .fai import Fai
        header = ["@HD\tVN:1.0\tSO:%s" % sort_order]
        for chrom, length in Fai(fasta):
            header.append("@SQ\tSN:%s\tLN:%i" % (chrom, length))
        return "\n".join(header) + "\n"

    def close(self):
        """Close the current file."""
        libhts.hts_close(self._htf)

    def flush(self):
        """Flush the current file."""
        libhts.bgzf_flush(self._htf.fp.bgzf)

    def __iter__(self):
        """Return self as an iterator of alignments."""
        return self

    def next(self):
        """Iterate over all alignments returning Alignment object."""
        libhts.sam_read1(self._htf, self.header._h, self._b)
        return Alignment(self._b, self.header._h)

    def __call__(self, region):
        """Perform a region query.
        
        Parameters
        ----------
        region : str
            region to query, e.g. chr1:1234-5678
        """
        qiter = libhts.sam_itr_querys(self._idx, self.header._h, region);
        try:
            slen = libhts.sam_itr_next(self._htf, qiter, self._b)

            while slen > 0:
                yield Alignment(self._b, self.header._h)
                slen = libhts.sam_itr_next(self._htf, qiter, self._b)

        finally:
            libhts.hts_itr_destroy(qiter)
