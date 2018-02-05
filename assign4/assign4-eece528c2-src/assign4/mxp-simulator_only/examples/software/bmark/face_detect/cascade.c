#include <stdio.h>
#include "vbx.h"
#include "haar.h"
#include "lbp.h"
#ifdef __NIOS2__
#include "io.h"
#include "system.h"
#else
#include "xparameters.h"
#endif
#include "vbx_port.h"

#define VLBPLUT VCUSTOM0
#define VLBPLUTWRITE VCUSTOM0
#define VLBPLUTZERO VCUSTOM1
#define VLBPLUTINC  VCUSTOM2

void cascade_print(lbp_stage_t *cascade, int max_stage)
{
    lbp_feat_t feat;
    int s, f, n, g, i, group, idx;
    int fcount = 0;

    for(s = 0; s < max_stage; s++) {
		for (f = 0; f < cascade[s].count; f++) {
            feat = cascade[s].feats[f];

            printf("feat %3d (%2d:%2d)\t", fcount, s, f);
            printf("@ %3d,%3d,%3d\t", feat.pos.src.x, feat.pos.src.y, feat.pos.size.x);
            printf("p/f %3d,%3d\t", feat.pass, feat.fail);
            printf("\n");

            printf("lut ");
            for (n = 0; n < 8; n++) {
                printf("%d ", feat.lut[n]);
            }
            printf("\n");

            printf("hex ");
            for (n = 0; n < 8; n++) {
                printf("%.8X ", feat.lut[n]);
            }
            printf("\n");

            printf("bin ");
#if 0
            for (n = 0; n < 8; n++) {
                for (i = 0; i < 32; i++) {
                    printf("%d", (feat.lut[n] >> (31-i)) & 0x1);
                }
                printf(" ");
            }
#else
            for (n = 0; n < 256; n++) {
                group = n >> 5;
                idx = n & 0x1f;
                if ((n % 32) == 0 && n > 0) {
                    printf(" ");
                }
                printf("%d", (feat.lut[group] >> (31-idx)) & 0x1);
            }
#endif
            printf("\n");

            fcount++;
        }
    }
}

int cascade_check_feature_2w(lbp_stage_t *cascade, int max_stage, int feat_num, unsigned *index_a, unsigned *index_b, int *result, int length, int print)
{
    lbp_feat_t feat_a, feat_b;
    int s, f, n, g, i;
    int group, idx, value_a, value_b;
    int fcount = 0;
    int errors = 0;

    for(s = 0; s < max_stage; s++) {
        for (f = 0; f < (cascade[s].count / 2); f++) {
            if(fcount == feat_num) {
                feat_a = cascade[s].feats[f*2];
                feat_b = cascade[s].feats[f*2+1];
                for (n = 0; n < length; n++) {
                    group = index_a[n] >> 5;
                    idx = index_a[n] & 0x1f;

                    if ((feat_a.lut[group] >> idx) & 0x1) {
                        value_a = feat_a.fail;
                    } else {
                        value_a = feat_a.pass;
                    }
                    group = index_b[n] >> 5;
                    idx = index_b[n] & 0x1f;

                    if ((feat_b.lut[group] >> idx) & 0x1) {
                        value_b = feat_b.fail;
                    } else {
                        value_b = feat_b.pass;
                    }

                    if ((value_a + value_b) != result[n]) {
                        if (print) {
                            printf("f %3d\t", fcount);
                            printf("pat_a[%3d]->%3d\t", n, index_a[n]);
                            printf("pat_b[%3d]->%3d\t", n, index_b[n]);
                            /* if (result[n] != feat.pass && result[n] != feat.fail) { */
                            if (1 || result[n] != (value_a + value_b)) {
                                printf("ans[%3d]->%3d != %3d (%3d + %3d)\ta[p,f = %3d,%3d] b[p,f = %3d,%3d]", 
                                        n, result[n], value_a + value_b, value_a, value_b, feat_a.pass, feat_a.fail, feat_b.pass, feat_b.fail);
                            }
                            printf("\n");
                        }
                        errors++;
                    }
                }
            }
            fcount++;
        }
    }
    return errors;
}

int cascade_check_feature_2b(lbp_stage_t *cascade, int max_stage, int feat_num, unsigned char *index_a, unsigned char *index_b, char *result, int length, int print)
{
    lbp_feat_t feat_a, feat_b;
    int s, f, n, g, i;
    int group, idx, value_a, value_b;
    int fcount = 0;
    int errors = 0;

    for(s = 0; s < max_stage; s++) {
        for (f = 0; f < (cascade[s].count / 2); f++) {
            if(fcount == feat_num) {
                feat_a = cascade[s].feats[f*2];
                feat_b = cascade[s].feats[f*2+1];
                for (n = 0; n < length; n++) {
                    group = index_a[n] >> 5;
                    idx = index_a[n] & 0x1f;

                    if ((feat_a.lut[group] >> idx) & 0x1) {
                        value_a = feat_a.fail;
                    } else {
                        value_a = feat_a.pass;
                    }
                    group = index_b[n] >> 5;
                    idx = index_b[n] & 0x1f;

                    if ((feat_b.lut[group] >> idx) & 0x1) {
                        value_b = feat_b.fail;
                    } else {
                        value_b = feat_b.pass;
                    }

                    if ((value_a + value_b) != result[n]) {
                        if (print) {
                            printf("f %3d\t", fcount);
                            printf("pat_a[%3d]->%3d\t", n, index_a[n]);
                            printf("pat_b[%3d]->%3d\t", n, index_b[n]);
                            /* if (result[n] != feat.pass && result[n] != feat.fail) { */
                            if (1 || result[n] != (value_a + value_b)) {
                                printf("ans[%3d]->%3d != %3d (%3d + %3d)\ta[p,f = %3d,%3d] b[p,f = %3d,%3d]", 
                                        n, result[n], value_a + value_b, value_a, value_b, feat_a.pass, feat_a.fail, feat_b.pass, feat_b.fail);
                            }
                            printf("\n");
                        }
                        errors++;
                    }
                }
            }
            fcount++;
        }
    }
    return errors;
}

int cascade_check_feature_2h(lbp_stage_t *cascade, int max_stage, int feat_num, unsigned short *index_a, unsigned short *index_b, short *result, int length, int print)
{
    lbp_feat_t feat_a, feat_b;
    int s, f, n, g, i;
    int group, idx, value_a, value_b;
    int fcount = 0;
    int errors = 0;

    for(s = 0; s < max_stage; s++) {
        for (f = 0; f < (cascade[s].count / 2); f++) {
            if(fcount == feat_num) {
                feat_a = cascade[s].feats[f*2];
                feat_b = cascade[s].feats[f*2+1];
                for (n = 0; n < length; n++) {
                    group = index_a[n] >> 5;
                    idx = index_a[n] & 0x1f;

                    if ((feat_a.lut[group] >> idx) & 0x1) {
                        value_a = feat_a.fail;
                    } else {
                        value_a = feat_a.pass;
                    }

                    group = index_b[n] >> 5;
                    idx = index_b[n] & 0x1f;

                    if ((feat_b.lut[group] >> idx) & 0x1) {
                        value_b = feat_b.fail;
                    } else {
                        value_b = feat_b.pass;
                    }

                    if ((value_a + value_b) != result[n]) {
                        if (print) {
                            printf("f %3d\t", fcount);
                            printf("pat_a[%3d]->%3d\t", n, index_a[n]);
                            printf("pat_b[%3d]->%3d\t", n, index_b[n]);
                            /* if (result[n] != feat.pass && result[n] != feat.fail) { */
                            if (result[n] != (value_a + value_b)) {
                                printf("ans[%3d]->%3d != %3d (%3d + %3d)\ta[p,f = %3d,%3d] b[p,f = %3d,%3d]", 
                                        n, result[n], value_a + value_b, value_a, value_b, feat_a.pass, feat_a.fail, feat_b.pass, feat_b.fail);
                            }
                            printf("\n");
                        }
                        errors++;
                    }
                }
            }
            fcount++;
        }
    }
    return errors;
}

int cascade_check_feature(lbp_stage_t *cascade, int max_stage, int feat_num, unsigned char *index, char *result, int length, int print)
{
    lbp_feat_t feat;
    int s, f, n, g, i;
    int group, idx, is_fail, value;
    int fcount = 0;
    int errors = 0;

    for(s = 0; s < max_stage; s++) {
		for (f = 0; f < cascade[s].count; f++) {
            if(fcount == feat_num) {
                feat = cascade[s].feats[f];
                for (n = 0; n < length; n++) {
                    group = index[n] >> 5;
                    idx = index[n] & 0x1f;
                    is_fail = (feat.lut[group] >> idx) & 0x1;
                    if (is_fail) {
                        value = feat.fail;
                    } else {
                        value = feat.pass;
                    }
                    if (value != result[n]) {
                        if (print) {
                            printf("f %3d\t", fcount);
                            printf("pat[%3d]->%3d\t", n, index[n]);
                            if (result[n] != feat.pass && result[n] != feat.fail) {
                                printf("ans[%3d]->%3d != %3d\t[p,f = %3d,%3d]", n, result[n], value, feat.pass, feat.fail);
                            }
                            printf("\n");
                        }
                        errors++;
                    }
                }
            }
            fcount++;
        }
    }
    return errors;
}

int cascade_check(lbp_stage_t *cascade, int max_stage, int length)
{
    lbp_feat_t feat;
    int s, f, n, g, i, group, idx;
    int fcount = 0;
    int errors = 0;

    unsigned char *check_pattern = (unsigned char *)vbx_shared_malloc(length * sizeof(char));
    char *check_lut     = (char *)vbx_shared_malloc(length * sizeof(char));
    vbx_ubyte_t* v_out  = (vbx_ubyte_t*)vbx_sp_malloc(length*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_in   = (vbx_ubyte_t*)vbx_sp_malloc(length*sizeof(vbx_ubyte_t));

    for (i = 0; i < length; i++) {
        check_pattern[i] = i & 0xff;
    }

    vbx_dma_to_vector(v_in, check_pattern, length);
    vbx_sync();

    vbx_set_vl(length);
    for(s = 0; s < max_stage; s++) {
        for (f = 0; f < cascade[s].count; f++) {
            if (fcount < 256) {
                vbx(SVBU, VLBPLUT, v_out, fcount & 0xff, v_in);
            } else {
                vbx(SVBS, VLBPLUT, v_out, fcount & 0xff, v_in);
            }
            vbx_dma_to_host(check_lut, v_out, length);
            vbx_sync();
            errors += cascade_check_feature(cascade, max_stage, fcount, check_pattern, check_lut, length, 0);
            fcount++;
        }
    }
#if 0
    for(s = 0; s < max_stage; s++) {
		for (f = 0; f < cascade[s].count; f++) {
            if (use_cvi) {
                if(fcount < 256) {
                    vbx(SVBU, VLBPLUT, (vptr_ubyte)v_lut, fcount & 0xff, v_pattern);
                } else {
                    vbx(SVBS, VLBPLUT, (vptr_ubyte)v_lut, fcount & 0xff, v_pattern);
                }
                vbx(VVB, VADD, v_stage, v_stage, v_lut);
            } else {
                vbx(SVB, VMOV, v_add, feat.fail, NULL);

                /* check if pattern is in lut */
                vbx(SVBU, VSHR, (vptr_ubyte)v_group, 5, v_pattern);
                int n;
                for (n = 0; n < 8; n++) {
                    vbx(SVB, VADD, (vptr_byte)v_sel, -n, (vptr_byte)v_group);
                    vbx(SVBW, VCMV_Z, (vptr_word)v_lut, feat.lut[n], (vptr_byte)v_sel);
                }

                vbx(SVBWU, VAND, (vptr_uword)v_idx, 0x1f, v_pattern);
                vbx(VVWB, VSHR, (vptr_byte)v_lut, (vptr_word)v_idx, (vptr_word)v_lut);
                vbx(SVB, VAND, (vptr_byte)v_lut, 1, (vptr_byte)v_lut);

                vbx(SVB, VCMV_LEZ, v_add, feat.pass, v_lut);
                vbx(VVB, VADD, v_stage, v_stage, v_add);
            }
            fcount++;
        }
    }
#endif
    vbx_sp_free();
    vbx_shared_free(check_lut);
    vbx_shared_free(check_pattern);

    return errors;
}

int cascade_check_2w(lbp_stage_t *cascade, int max_stage, int length)
{
    lbp_feat_t feat;
    int s, f, n, g, i, group, idx;
    int fcount = 0, errors = 0;

    unsigned int *check_pattern_a = (unsigned int*)vbx_shared_malloc(length * sizeof(int));
    unsigned int *check_pattern_b = (unsigned int*)vbx_shared_malloc(length * sizeof(int));
    int *check_lut     = (int*)vbx_shared_malloc(length * sizeof(int));
    vbx_word_t* v_out  = (vbx_word_t*)vbx_sp_malloc(length*sizeof(vbx_word_t));
    vbx_uword_t* v_in_a   = (vbx_uword_t*)vbx_sp_malloc(length*sizeof(vbx_uword_t));
    vbx_uword_t* v_in_b   = (vbx_uword_t*)vbx_sp_malloc(length*sizeof(vbx_uword_t));

    for (i = 0; i < length; i++) {
        check_pattern_b[i] = i & 0xff;
        check_pattern_a[i] = (0xff - i) & 0xff;
    }

    vbx_dma_to_vector(v_in_a, check_pattern_a, length*4);
    vbx_dma_to_vector(v_in_b, check_pattern_b, length*4);

    vbx_set_vl(length);
    for(s = 0; s < max_stage; s++) {
        for (f = 0; f < cascade[s].count /2; f++) {
            /* printf("stage %d, stage feat %d, abs feat %d\n", s, f, fcount); */
            if (fcount == 0) {
                vbx(VVWU, VLBPLUTZERO, v_out, v_in_a, v_in_b);
            } else {
                vbx(VVWU, VLBPLUTINC, v_out, v_in_a, v_in_b);
            }
            vbx_dma_to_host(check_lut, v_out, length*4);
            vbx_sync();
            errors += cascade_check_feature_2w(cascade, max_stage, fcount, check_pattern_a, check_pattern_b, check_lut, length, 1);
            fcount++;
        }
    }

    vbx_sp_free();
    vbx_shared_free(check_lut);
    vbx_shared_free(check_pattern_a);
    vbx_shared_free(check_pattern_b);

    return errors;
}
int cascade_check_2b(lbp_stage_t *cascade, int max_stage, int length)
{
    lbp_feat_t feat;
    int s, f, n, g, i, group, idx;
    int fcount = 0, errors = 0;

    unsigned char *check_pattern_a = (unsigned char *)vbx_shared_malloc(length * sizeof(char));
    unsigned char *check_pattern_b = (unsigned char *)vbx_shared_malloc(length * sizeof(char));
    char *check_lut     = (char *)vbx_shared_malloc(length * sizeof(char));
    vbx_ubyte_t* v_out  = (vbx_ubyte_t*)vbx_sp_malloc(length*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_in_a   = (vbx_ubyte_t*)vbx_sp_malloc(length*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_in_b   = (vbx_ubyte_t*)vbx_sp_malloc(length*sizeof(vbx_ubyte_t));

    for (i = 0; i < length; i++) {
        check_pattern_a[i] = i & 0xff;
        check_pattern_b[i] = (0xff - i) & 0xff;
    }

    vbx_dma_to_vector(v_in_a, check_pattern_a, length);
    vbx_dma_to_vector(v_in_b, check_pattern_b, length);

    vbx_set_vl(length);
    for(s = 0; s < max_stage; s++) {
        for (f = 0; f < cascade[s].count /2; f++) {
            if (fcount == 0) {
                vbx(VVBU, VLBPLUTZERO, v_out, v_in_a, v_in_b);
            } else {
                vbx(VVBU, VLBPLUTINC, v_out, v_in_a, v_in_b);
            }
            vbx_dma_to_host(check_lut, v_out, length);
            vbx_sync();
            errors += cascade_check_feature_2b(cascade, max_stage, fcount, check_pattern_a, check_pattern_b, check_lut, length, 1);
            fcount++;
        }
    }

    vbx_sp_free();
    vbx_shared_free(check_lut);
    vbx_shared_free(check_pattern_a);
    vbx_shared_free(check_pattern_b);

    return errors;
}

int cascade_check_2h(lbp_stage_t *cascade, int max_stage, int length)
{
    lbp_feat_t feat;
    int s, f, n, g, i, group, idx;
    int fcount = 0, errors = 0;

    unsigned short *check_pattern_a = (unsigned short *)vbx_shared_malloc(length * sizeof(short));
    unsigned short *check_pattern_b = (unsigned short *)vbx_shared_malloc(length * sizeof(short));
    short *check_lut     = (short *)vbx_shared_malloc(length * sizeof(short));
    vbx_half_t* v_out  = (vbx_half_t*)vbx_sp_malloc(length*sizeof(vbx_half_t));
    vbx_uhalf_t* v_in_a   = (vbx_uhalf_t*)vbx_sp_malloc(length*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_in_b   = (vbx_uhalf_t*)vbx_sp_malloc(length*sizeof(vbx_uhalf_t));

    for (i = 0; i < length; i++) {
        check_pattern_b[i] = i & 0xff;
        check_pattern_a[i] = (0xff - i) & 0xff;
    }

    vbx_dma_to_vector(v_in_a, check_pattern_a, length*2);
    vbx_dma_to_vector(v_in_b, check_pattern_b, length*2);

    vbx_set_vl(length);
    for(s = 0; s < max_stage; s++) {
        for (f = 0; f < cascade[s].count /2; f++) {
            /* printf("stage %d, stage feat %d, abs feat %d\n", s, f, fcount); */
            if (fcount == 0) {
                vbx(VVHU, VLBPLUTZERO, v_out, v_in_a, v_in_b);
            } else {
                vbx(VVHU, VLBPLUTINC, v_out, v_in_a, v_in_b);
            }
            vbx_dma_to_host(check_lut, v_out, length*2);
            vbx_sync();
            errors += cascade_check_feature_2h(cascade, max_stage, fcount, check_pattern_a, check_pattern_b, check_lut, length, 1);
            fcount++;
        }
    }

    vbx_sp_free();
    vbx_shared_free(check_lut);
    vbx_shared_free(check_pattern_a);
    vbx_shared_free(check_pattern_b);

    return errors;
}

int cascade_max_feature()
{
    int s, fcount = 0;
    for (s = 0; s < face_lbp_max_stage; s++) {
        fcount += face_lbp[s].count;
    }
    return fcount;
}

int lut_check(const int length, const int id, const int print, const int kill)
{
    int f, z;
    int errors = 0;
    int max_feat = cascade_max_feature();
    /* max_feat = 2; */
    int nop_iterations = 1;


    char lut[max_feat][length];
    char *out = (char *)vbx_shared_malloc(length * sizeof(char));
    unsigned char *in = (unsigned char *)vbx_shared_malloc(length * sizeof(unsigned char));
    unsigned char *pattern = (unsigned char *)vbx_shared_malloc(length * sizeof(unsigned char));

    vbx_ubyte_t *v_out   = (vbx_ubyte_t*)vbx_sp_malloc(length);
    vbx_ubyte_t *v_in    = (vbx_ubyte_t*)vbx_sp_malloc(length);

#if 0
    for (z = 0; z < length; z++) {
        in[z] = (z & 0xff);
    }
#else
    unsigned char lfsr = 0x73;
    unsigned char bit;
    for (z = 0; z < length; z++) {
        in[z] = lfsr;
        pattern[z] = lfsr;
        bit  = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 4) ) & 1;
        lfsr =  (lfsr >> 1) | (bit << 7);
    }
#endif

#if 0
    vbx_dma_to_vector(v_in, in, length);

    vbx_set_vl(length);

    /* vbx(SEBU, VADD, v_in, 0, 0); */

    for (f = 256; f < 257; f++) {

        for (z = 0; z < nop_iterations; z++) {
            vbx(SVBU, VADD, v_in, 0, v_in);
        }

        if (f < 256) {
            vbx(SVBU, VLBPLUT, v_out, (f & 0xff), v_in);
        } else {
            vbx(SVBS, VLBPLUT, v_out, (f & 0xff), v_in);
        }

        vbx_dma_to_host(out, v_out, length);
        vbx_dma_to_host(in, v_in, length);
        vbx_sync();
        for (z = 0; z < length; z++) {
            lut[f][z] = out[z];
        }
    }

#if 0
    vbx_dma_to_host(in, v_in, length);
#elif 0
    for (z = 0; z < length; z++) {
        in[z] = z;
    }
#endif

    vbx_sync();
    vbx_sp_free();

#if 1
    for (f = 0; f < max_feat; f++) {
        errors += cascade_check_feature(face_lbp, face_lbp_max_stage, f, in, lut[f], length, print);
    }
#endif
#else
    vbx_dma_to_vector(v_in, in, length);
    vbx_set_vl(length);
    vbx(SVBU, VADD, v_in, 0, v_in);
    vbx(SVBS, VLBPLUT, v_out, 0, v_in);
    vbx_dma_to_host(in, v_in, length);
    vbx_sync();
    vbx_sp_free();
#endif

#if 1
    //check input
    for (z = 0; z < length; z++) {
        if (in[z] != pattern[z]) {
            errors++;
            /* printf("in @ %d\t%d != %d\n", z, in[z], pattern[z]); */
        }
    }
#endif

    if (kill) {
        printf("lut check %d: errors %d\n", id, errors);
        putchar(4);
        while(1);
    } 

    return errors;
}

int lut_check2(const int input_length, const int lut_length, const int lut_iterations, const int print)
{
    /* fails w/ input_length = lut_length (try 128) and lut_iterations = 13 */
    int i, z, errors = 0;
    int use_test_pattern = 0;
    unsigned char lfsr = 0x73;
    unsigned char bit;

    unsigned char      *in = (unsigned char *)vbx_shared_malloc(input_length*sizeof(unsigned char));
    unsigned char *pattern = (unsigned char *)vbx_shared_malloc(input_length*sizeof(unsigned char));
    vbx_ubyte_t *v_in    = (vbx_ubyte_t*)vbx_sp_malloc(input_length);

    if (in == NULL || pattern == NULL || v_in == NULL) {
        printf("error allocating\n");
    }

    /* generate test pattern or enumerated pattern */
    for (z = 0; z < input_length; z++) {
        if (use_test_pattern) {
            in[z] = lfsr;
            pattern[z] = lfsr;
            bit  = ((lfsr >> 0) ^ (lfsr >> 2) ^ (lfsr >> 3) ^ (lfsr >> 4) ) & 1;
            lfsr =  (lfsr >> 1) | (bit << 7);
        } else {
            in[z] = z;
            pattern[z] = z;
        }
    }

    /* run lbp lut custom instruction */
    vbx_set_vl(lut_length);
    for (i = 0; i < lut_iterations; i++) {
#if 0
        vbx(SVBU, VCUSTOM0, 0, 0, 0);
#else
        vbx(SVBU, VCUSTOM1, 0, 0, 0);
#endif
    }
    vbx_sync();

    /* dma in, add 0, dma out */
    vbx_dma_to_vector(v_in, in, input_length);

    vbx_set_vl(input_length);
    vbx(SVB, VADD, v_in, 0, v_in);

    vbx_dma_to_host(in, v_in, input_length);

    vbx_sync();
    vbx_sp_free();

    /* check input/pattern */
    for (z = 0; z < input_length; z++) {
        if (in[z] != pattern[z]) {
            errors++;
            if (print) {
                printf("in @ %d\t%d != %d\n", z, in[z], pattern[z]);
            }
        }
    }

    return errors;
}
