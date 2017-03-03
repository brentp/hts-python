import os
from hts import Tbx
from nose.tools import assert_raises

HERE = os.path.dirname(__file__)
GTF = os.path.join(HERE, "example.gtf.gz")
BAM = os.path.join(HERE, "e.sam")

def test_build():
    if os.path.exists(GTF + ".tbi"):
        os.unlink(GTF + ".tbi")
    t = Tbx(GTF)
    assert t._tbx 
    assert os.path.exists(GTF + ".tbi")

def test_error_on_bad_file():
    assert_raises(Exception, Tbx, BAM)

