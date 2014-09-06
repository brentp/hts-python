from __future__ import print_function, division
from cffi import FFI
import os.path as op
import atexit

def _raise_if_null(v, msg):
    if v == ffi.NULL:
        raise Exception(msg)

ffi = FFI()

HERE = op.relpath(op.dirname(__file__))

with open(op.join(HERE, "hts_concat.h")) as headers:
    ffi.cdef(headers.read())
with open(op.join(HERE, "hts_extra.h")) as headers:
    ffi.cdef(headers.read())


libhts = ffi.verify('''
#include "stdlib.h"
#include <zlib.h>
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kfunc.h"
#include "hts_extra.h"
''',
    libraries=['c', 'z', 'hts'],
    ext_package='htsffii',
    depends=["%s/hts_extra.h" % HERE],
    sources=["%s/hts_extra.c" % HERE],
    include_dirs=[HERE, "/usr/include/", "/usr/local/include/"],
)
