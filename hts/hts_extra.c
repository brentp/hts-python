#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kstring.h"


// the htslib api only provides bam_seqi to get one base at a time.
// to avoid having to call that, e.g. 100 times to get the sequence for
// a 100 base read, we implement this.
int bam_get_read_seq(bam1_t *b, kstring_t * str){
    int i = 0;
    str->l = 0;
    uint8_t *s = bam_get_seq(b);
    for (; i < b->core.l_qseq; ++i) kputc("=ACMGRSVTWYHKDBN"[bam_seqi(s, i)], str);
    return i;
}

int as_gts(int *gts, int num_samples) {
	int j = 0, i;
	for (i = 0; i < 2 * num_samples; i += 2){
	    if (bcf_gt_is_missing(gts[i]) && bcf_gt_is_missing(gts[i+1])){
			gts[j++] = 2; // unknown
			continue;
		}
		int a = bcf_gt_allele(gts[i]);
		int b = bcf_gt_allele(gts[i+1]);

	    if((a == 0) && (b == 0)) {
			gts[j] = 0; //  HOM_REF
		}
		else if((a == 1) && (b == 1)) {
			gts[j] = 3; //  HOM_ALT
		}
		else if((a  != b)) {
			gts[j] = 1; //  HET
		} else {
			gts[j] = 2; // unknown
		}
		j++;
	}
	// free the latter part of the array that we don't need.
	//free((void *)(&(gts[j])));
	return j;
}


static inline int cigar_iref2iseq_set(uint32_t **cigar, uint32_t *cigar_max, int *icig, int *iseq, int *iref)
{
    int pos = *iref;
    if ( pos < 0 ) return -1;
    *icig = 0;
    *iseq = 0;
    *iref = 0;
    while ( *cigar<cigar_max )
    {
        int cig  = (**cigar) & BAM_CIGAR_MASK;
        int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

        if ( cig==BAM_CSOFT_CLIP ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) { (*cigar)++; *icig = 0; continue; }
        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            pos -= ncig;
            if ( pos < 0 ) { *icig = ncig + pos; *iseq += *icig; *iref += *icig; return BAM_CMATCH; }
            (*cigar)++; *iseq += ncig; *icig = 0; *iref += ncig;
            continue;
        }
        if ( cig==BAM_CINS ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CDEL || cig==BAM_CREF_SKIP )
        {
            pos -= ncig;
            if ( pos<0 ) pos = 0;
            (*cigar)++; *icig = 0; *iref += ncig;
            continue;
        }
    }
    *iseq = -1;
    return -1;
}
static inline int cigar_iref2iseq_next(uint32_t **cigar, uint32_t *cigar_max, int *icig, int *iseq, int *iref)
{
    while ( *cigar < cigar_max )
    {
        int cig  = (**cigar) & BAM_CIGAR_MASK;
        int ncig = (**cigar) >> BAM_CIGAR_SHIFT;

        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            if ( *icig >= ncig - 1 ) { *icig = 0;  (*cigar)++; continue; }
            (*iseq)++; (*icig)++; (*iref)++;
            return BAM_CMATCH;
        }
        if ( cig==BAM_CDEL || cig==BAM_CREF_SKIP ) { (*cigar)++; (*iref) += ncig; *icig = 0; continue; }
        if ( cig==BAM_CINS ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CSOFT_CLIP ) { (*cigar)++; *iseq += ncig; *icig = 0; continue; }
        if ( cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) { (*cigar)++; *icig = 0; continue; }
        fprintf(stderr,"todo: cigar %d\n", cig);
    }
    *iseq = -1;
    *iref = -1;
    return -1;
}

void tweak_overlap_quality(bam1_t *a, bam1_t *b) {
    uint32_t *a_cigar = bam_get_cigar(a), *a_cigar_max = a_cigar + a->core.n_cigar;
    uint32_t *b_cigar = bam_get_cigar(b), *b_cigar_max = b_cigar + b->core.n_cigar;
    int a_icig = 0, a_iseq = 0;
    int b_icig = 0, b_iseq = 0;
    uint8_t *a_qual = bam_get_qual(a), *b_qual = bam_get_qual(b);
    uint8_t *a_seq  = bam_get_seq(a), *b_seq = bam_get_seq(b);

    int iref   = b->core.pos;
    int a_iref = iref - a->core.pos;
    int b_iref = iref - b->core.pos;
    int a_ret = cigar_iref2iseq_set(&a_cigar, a_cigar_max, &a_icig, &a_iseq, &a_iref);
    if ( a_ret<0 ) return;  // no overlap
    int b_ret = cigar_iref2iseq_set(&b_cigar, b_cigar_max, &b_icig, &b_iseq, &b_iref);
    if ( b_ret<0 ) return;  // no overlap


    while ( 1 )
    {
        // Increment reference position
        while ( a_iref>=0 && a_iref < iref - a->core.pos )
            a_ret = cigar_iref2iseq_next(&a_cigar, a_cigar_max, &a_icig, &a_iseq, &a_iref);
        if ( a_ret<0 ) break;   // done
        if ( iref < a_iref + a->core.pos ) iref = a_iref + a->core.pos;

        while ( b_iref>=0 && b_iref < iref - b->core.pos )
            b_ret = cigar_iref2iseq_next(&b_cigar, b_cigar_max, &b_icig, &b_iseq, &b_iref);
        if ( b_ret<0 ) break;   // done
        if ( iref < b_iref + b->core.pos ) iref = b_iref + b->core.pos;

        iref++;
        if ( a_iref+a->core.pos != b_iref+b->core.pos ) continue;   // only CMATCH positions, don't know what to do with indels

        if ( bam_seqi(a_seq,a_iseq) == bam_seqi(b_seq,b_iseq) )
        {
            // we are very confident about this base
            int qual = a_qual[a_iseq] + b_qual[b_iseq];
            a_qual[a_iseq] = qual>200 ? 200 : qual;
            b_qual[b_iseq] = 0;
        }
        else
        {
            if ( a_qual[a_iseq] >= b_qual[b_iseq] )
            {
                a_qual[a_iseq] = 0.8 * a_qual[a_iseq];  // not so confident about a_qual anymore given the mismatch
                b_qual[b_iseq] = 0;
            }
            else
            {
                b_qual[b_iseq] = 0.8 * b_qual[b_iseq];
                a_qual[a_iseq] = 0;
            }
        }
    }
}


int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}


int skip_aux(uint8_t *s)
{
    int size = aux_type2size(*s); ++s; // skip type
    uint32_t n;
	int k = 0;
    switch (size) {
    case 'Z':
    case 'H':
        while (*s) {++s;++k;}
        return k + 1;
    case 'B':
        size = aux_type2size(*s);
		++s;
        memcpy(&n, s, 4); s += 4;
        return n * size + 4;
    case 0:
        abort();
        break;
    default:
        return size;
    }
}

