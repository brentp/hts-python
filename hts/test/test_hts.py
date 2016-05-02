import os
from hts import Bam

HERE = os.path.dirname(__file__)

SAM = os.path.join(HERE, "e.sam")
AUX_SAM = os.path.join(HERE, "t.sam")

def test_aux():

    b = Bam(SAM)
    r = next(b)
    assert r.tags == [('NH', 'C', 1)], r.aux

def test_tname():
    b = Bam(SAM)
    r = next(b)
    assert r.target is None



def test_aux_bam():

    bam = Bam(AUX_SAM)

    aln = next(bam)
    assert aln.tags == [('MC', 'Z', '101M'), ('MD', 'Z', '16T2C80'), ('PG', 'Z', 'bwa-meth'), ('RG', 'Z', '44_Mm08_WEAd_Db2_WGBS_E_1_L001__trimmed'), ('NM', 'C', 28), ('MQ', 'C', 60), ('UQ', 'S', 1064), ('AS', 'C', 87), ('XS', 'C', 101)], aln.tags
    aln = next(bam)
    assert aln.tags == [('MC', 'Z', '7M1D94M'), ('MD', 'Z', '9T2C87'), ('PG', 'Z', 'bwa-meth'), ('RG', 'Z', '44_Mm08_WEAd_Db2_WGBS_E_1_L001__trimmed'), ('NM', 'C', 27), ('MQ', 'C', 25), ('UQ', 'S', 981), ('AS', 'C', 87), ('XS', 'C', 101)]
