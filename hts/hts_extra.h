int bam_get_read_seq(bam1_t *b, kstring_t * str);

void tweak_overlap_quality(bam1_t *, bam1_t *);

int as_gts(int *gts, int num_samples);


int aux_type2size(uint8_t type);

int skip_aux(uint8_t *s);

