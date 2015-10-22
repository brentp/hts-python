import os
from hts import Bam

HERE = os.path.dirname(__file__)

SAM = os.path.join(HERE, "e.sam")

def test_aux():

    b = Bam(SAM)
    r = next(b)
    assert r.tags == [('NH', 'C', 1)], r.aux

def test_tname():
    b = Bam(SAM)
    r = next(b)
    assert r.target is None

