/* VECTORBLOX MXP SOFTWARE DEVELOPMENT KIT
 *
 * Copyright (C) 2012-2016 VectorBlox Computing Inc., Vancouver, British Columbia, Canada.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 *     * Neither the name of VectorBlox Computing Inc. nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * This agreement shall be governed in all respects by the laws of the Province
 * of British Columbia and by the laws of Canada.
 *
 * This file is part of the VectorBlox MXP Software Development Kit.
 *
 */


#include "vbx_copyright.h"
VBXCOPYRIGHT(test_haar)

#include "test.h"
#include "auxillary.h"

//FIXME stride for match not implemented
int compare_LBPPassStage_to_restricted(unsigned short *vbx_img, int log, lbp_stage_t lbp_stage, int window, int width, int height, int max_print_errors)
{
    int l, i, j, cell, errors = 0;

    unsigned char** scalar_patterns = test_scalar_patterns(vbx_img, log, width, height);

    unsigned char *pass, *vbx_pass;
    pass = (unsigned char*)vbx_shared_malloc(width*height*sizeof(unsigned char));
    vbx_pass = (unsigned char*)vbx_shared_malloc(width*height*sizeof(unsigned char));
    
    vbx_byte_t** v_lbp =(vbx_byte_t**)vbx_shared_malloc((log+1)*sizeof(vbx_byte_t*));
    for (l=0; l<log+1; l++) {
        v_lbp[l] = (vbx_byte_t*)vbx_sp_malloc((window+1)*width*sizeof(vbx_byte_t)); 
    }
    vbx_byte_t* v_lut = (vbx_byte_t*)vbx_sp_malloc(width*sizeof(vbx_byte_t)); 
    vbx_byte_t* v_stage = (vbx_byte_t*)vbx_sp_malloc(width*sizeof(vbx_byte_t)); 
    vbx_byte_t* v_pattern;
    lbp_feat_t feat;
    int dx, dy, dw, f;

    for (l=0; l<log+1; l++) {
        vbx_dma_to_vector(v_lbp[l]+width, scalar_patterns[l], (window)*width*sizeof(unsigned char));
    }
    vbx_sync();
    for(j=0; j < height-(window+1); j++) {
        for (l=0; l<log+1; l++) {
            vbx_set_vl(width * window);
            vbx(VVB, VMOV, v_lbp[l], v_lbp[l]+width, NULL);
            vbx_dma_to_vector(v_lbp[l] + window*width, scalar_patterns[l]+(j+window)*width, width*sizeof(unsigned char));
        }

        vbx_set_vl(width-(window+1));
        vbx(SVB, VMOV, v_stage, 0, NULL);
        for (f = 0; f < lbp_stage.count; f++) {
            feat = lbp_stage.feats[f];
            dx = feat.pos.src.x;
            dy = feat.pos.src.y;
            dw = feat.pos.size.x;
            v_pattern = v_lbp[dw>>1]+(dy*width+dx);

            vbx(SVBU, VLBPLUT, v_lut, f, v_pattern);
            vbx(VVB, VADD, v_stage, v_stage, v_lut);
        }
        vbx(SVB, VMOV, v_lut, 0, NULL);
        vbx(SVB, VCMV_GEZ, v_lut, 1, v_stage);
        vbx_dma_to_host(vbx_pass + j*width, v_lut, (width-(window+1))*sizeof(unsigned char));
        vbx_sync();
    }


    unsigned int *iImg, *iiImg;
    iImg = (unsigned int *)vbx_shared_malloc(width*height*sizeof(unsigned int));
    iiImg = (unsigned int *)vbx_shared_malloc(width*height*sizeof(unsigned int));

    gen_integrals(vbx_img, iImg, iiImg, width, height);

    image_t lbp_img = {iImg, {width, height}};
    for (j = 0; j < height - (window + 1); j++) {
        for (i = 0; i < width - (window + 1); i++) {
            pair_t lbp_p = {i, j};
            pass[j*width+i] = LBPPassStage(lbp_img, lbp_stage, lbp_p);
        }
    }

    /* test pass vs vbx pass */
    for (j = 0; j < height - (window + 1); j++) {
        errors += match_array_byte(vbx_pass + j*width, pass + j*width, "pass stage", width - (window + 1), 1, 0, max_print_errors, 1, j);
        if (errors > max_print_errors){
            max_print_errors = 0;
        }
    }
    return errors;
}
int compare_scalar_BLIP2_to_vector_BLIP(unsigned short* img, pixel* vbx_input, int width, int height, int max_print_errors, int scale_factor)
{
    int j, errors = 0;
    int scaled_width, scaled_height;

    /* scale facetor v/v+1, v is between 1-10 */
    scaled_width = width*scale_factor/(scale_factor+1);
    scaled_height= height*scale_factor/(scale_factor+1);

    unsigned short *scaled_img, *vbx_img, *vbx_scaled_img; 
    unsigned char *vbx_img8, *vbx_scaled_img8;
    unsigned int *iImg, *iiImg, *vbx_iImg, *vbx_iiImg;

    scaled_img = (unsigned short*)vbx_shared_malloc(scaled_width*scaled_height*sizeof(unsigned short));

    iImg = (unsigned int*)vbx_shared_malloc(scaled_width*scaled_height*sizeof(unsigned int));
    iiImg = (unsigned int*)vbx_shared_malloc(scaled_width*scaled_height*sizeof(unsigned int));

    vbx_scaled_img = (unsigned short*)vbx_shared_malloc(scaled_width*scaled_height*sizeof(unsigned short));
    vbx_img = (unsigned short*)vbx_shared_malloc(width*height*sizeof(unsigned short));
    vbx_img8 = (unsigned char*)vbx_shared_malloc(width*height*sizeof(unsigned char));
    vbx_scaled_img8 = (unsigned char*)vbx_shared_malloc(scaled_width*scaled_height*sizeof(unsigned char));

    vbx_iImg = (unsigned int*)vbx_shared_malloc(width*height*sizeof(unsigned int));
    vbx_iiImg = (unsigned int*)vbx_shared_malloc(width*height*sizeof(unsigned int));


#if 0
    scalar_BLIP2(img, height, width, scaled_img, scaled_height, scaled_width, scale_factor);
#else
    float percent = 1.0 * (scale_factor+1) / scale_factor;
    scalar_BLIP(img, height, width, scaled_img, scaled_height, scaled_width, &percent);
#endif
    gen_integrals(scaled_img, iImg, iiImg, scaled_width, scaled_height);

	vector_get_img(vbx_img, vbx_iImg, vbx_iiImg, vbx_input, 1, width, height, width, 1);
    vector_BLIP(vbx_img, height, width, vbx_scaled_img, vbx_iImg, vbx_iiImg, scaled_height, scaled_width, scale_factor, 1);
	vector_get_img8(vbx_img8, vbx_input, 1, width, height, width);
    /* vector_BLIP8(vbx_img8, height, width, vbx_scaled_img8, scaled_height, scaled_width, scale_factor); */
	vbx_timestamp_start();
	vbx_timestamp_t time_start, time_stop;
    double vbx_time;
    time_start = vbx_timestamp();
#if 1
    vector_BLIP8F3(vbx_img8, height, width, vbx_scaled_img8, scaled_height, scaled_width, scale_factor);
#else
    vector_BLIP8F2(vbx_img8, height, width, vbx_scaled_img8, scaled_height, scaled_width, scale_factor);
#endif
    time_stop = vbx_timestamp();
    vbx_time = vbx_print_vector_time(time_start, time_stop, 0.0);

    /* test greyscale image */
    for (j = 0; j < height; j++) {
        errors += match_array_half(img+j*width, vbx_img+j*width, "greyscale", width, 1, 0, max_print_errors, j);
        if(errors > max_print_errors) {
            max_print_errors = 0;
        }
    }

    /* test scaled image */
    for (j = 0; j < scaled_height; j++) {
        errors += match_array_half(scaled_img+j*scaled_width, vbx_scaled_img+j*scaled_width, "scaled greyscale", scaled_width, 1, 1, max_print_errors, j);
        if(errors > max_print_errors) {
            max_print_errors = 0;
        }
    }

    for (j = 0; j < scaled_height; j++) {
        errors += match_array_half_byte(scaled_img+j*scaled_width, vbx_scaled_img8+j*scaled_width, "scaled greyscale8", scaled_width, 1, 1, max_print_errors, j);
        if(errors > max_print_errors) {
            max_print_errors = 0;
        }
    }

#if 0
    /* test scaled_integral image */
    for (j = 0; j < scaled_height; j++) {
        errors += match_array_word(iImg+j*scaled_width, vbx_iImg+j*scaled_width, "scaled integral", scaled_width, 1, 0, max_print_errors, j);
        if(errors > max_print_errors) {
            max_print_errors = 0;
        }
    }

    /* test scaled squared integral image */
    for (j = 0; j < scaled_height; j++) {
        errors += match_array_word(iiImg+j*scaled_width, vbx_iiImg+j*scaled_width, "scaled squared", scaled_width, 1, 0, max_print_errors, j);
        if(errors > max_print_errors) {
            max_print_errors = 0;
        }
    }
#endif

    /* test scaled_integral image */
    return errors;
}
int compare_gen_integrals_to_vector_get_img(unsigned short* img, unsigned int* iImg, unsigned int* iiImg, unsigned short* vbx_img, unsigned int* vbx_iImg, unsigned int* vbx_iiImg, pixel* vbx_input, int width, int height, int max_print_errors)
{
    int j, errors = 0;

    gen_integrals(img, iImg, iiImg, width, height);
    /* bin is 1*/
	vector_get_img(vbx_img, vbx_iImg, vbx_iiImg, vbx_input, 1, width, height, width, 1);

    /* test greyscale image */
    for (j = 0; j < height; j++) {
        errors += match_array_half(img+j*width, vbx_img+j*width, "greyscale", width, 1, 0, max_print_errors, j);
        if(errors > max_print_errors) {
            max_print_errors = 0;
        }
    }

    /* test integral image */
    for (j = 0; j < height; j++) {
        errors += match_array_word(iImg+j*width, vbx_iImg+j*width, "integral", width, 1, 0, max_print_errors, j);
        if(errors > max_print_errors) {
            max_print_errors = 0;
        }
    }

    /* test squared integral image */
    for (j = 0; j < height; j++) {
        errors += match_array_word(iiImg+j*width, vbx_iiImg+j*width, "squared", width, 1, 0, max_print_errors, j);
        if(errors > max_print_errors) {
            max_print_errors = 0;
        }
    }
    return errors;
}

int compare_LBPRestricted_to_test_scalar_patterns(unsigned short* vbx_img, int log, int width, int height, int max_print_errors)
{
    int j, l, cell, errors = 0;

    /* generate patterns */
    unsigned char **patterns = LBPRestricted(vbx_img, width, height, log);
    unsigned char** scalar_patterns = test_scalar_patterns(vbx_img, log, width, height);


    /* test sums vs scalar sums */
    for (l=0; l<log+1; l++) {
        cell = 1 << l;
        for (j = 0; j < height - (3*cell-1); j++) {
            errors += match_array_byte(scalar_patterns[l]+j*width, patterns[l]+j*width, "restricted patterns", width - (3*cell-1), 1, 0, max_print_errors, 1, j);
            if(errors > max_print_errors) {
                max_print_errors = 0;
            }
        }
    }
    return errors;
}

int compare_LBPRestrictedCI_to_test_scalar_patterns(unsigned short* vbx_img, unsigned char* vbx_img8, int log, int width, int height, int max_print_errors)
{
    int i, j, l, cell, errors = 0, cell_errors, row_errors;

    /* generate patterns */
    unsigned char **patterns = (unsigned char**)malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        patterns[l] = (unsigned char*)vbx_shared_malloc(height*width*sizeof(unsigned char));
    }

    /* LBPRestrictedCI28(patterns, vbx_img8, width, height, log); */
    LBPRestricted_CI_column_8(patterns, vbx_img8, width, height, log);
    /* LBPRestricted_CI_column_8_scratch(patterns, vbx_img8, width, height, log); */
    unsigned char **scalar_patterns = test_scalar_patterns(vbx_img, log, width, height);

    /* test sums vs scalar sums */
    for (l=0; l<log+1; l++) {
        cell = 1 << l;
        cell_errors = 0;
        printf("Testing cell %d\n", cell);
        for (j = 0; j < height - (3*cell-1); j++) {
            row_errors = match_array_byte(scalar_patterns[l]+j*width, patterns[l]+j*width, "restricted patterns", width - (3*cell-1), 1, 0, max_print_errors, 1, j);
            if (row_errors) {
                printf("errors in row %d\n", j);
            }
            cell_errors += row_errors;
            errors += row_errors;
            if(errors > max_print_errors) {
                max_print_errors = 0;
            }
        }
        printf("Total errors: %d\n\n", cell_errors);
    }

    return errors;

}

int compare_LBPRestrictedPatterns2_to_test_scalar_patterns(unsigned short* vbx_img, int log, int width, int height, int max_print_errors)
{
    int cell, errors = 0;
    int i, j;

    /* generate patterns */
    unsigned short **sums = LBPRestrictedSums2(vbx_img, width, height, log);
    unsigned char **patterns = LBPRestrictedPatterns2(vbx_img, sums, width, height, log);

    unsigned char** scalar_patterns = test_scalar_patterns(vbx_img, log, width, height);

    /* test sums vs scalar sums */
    for (j=0; j<log+1; j++) {
        cell = 1 << j;
        for (i = 0; i < height - (3*cell-1); i++) {
            errors += match_array_byte(scalar_patterns[j]+i*width, patterns[j]+i*width, "patterns", width - (3*cell-1), 1, 0, max_print_errors, 1, i);
            if(errors > max_print_errors) {
                max_print_errors = 0;
            }
        }
    }
    return errors;
}

int compare_LBPRestrictedPatterns_to_test_scalar_patterns(unsigned short* vbx_img, int log, int width, int height, int max_print_errors)
{
    int cell, errors = 0;
    int i, j;

    /* generate patterns */
    unsigned char **sums = LBPRestrictedSums(vbx_img, width, height, log);
    unsigned char **patterns = LBPRestrictedPatterns(sums, width, height, log);

    unsigned char** scalar_patterns = test_scalar_patterns(vbx_img, log, width, height);

    /* test patterns vs scalar patterns */
    for (j=0; j<log+1; j++) {
        cell = (1 << j);
        for (i = 0; i < height - (3*cell-1); i++) {
            errors += match_array_byte(scalar_patterns[j]+i*width, patterns[j]+i*width, "patterns", width - (3*cell-1), 1, 0, max_print_errors, 1, i);
            if(errors > max_print_errors) {
                max_print_errors = 0;
            }
        }
    }
    return errors;
}

int compare_ScalarLBPRestrictedPatterns_to_SATBinaryPattern(unsigned short *vbx_img, int log, int width, int height, int max_print_errors)
{
    int l, i, j, cell, errors = 0;

    /* generate patterns */
    unsigned short **sums = ScalarLBPRestrictedSums(vbx_img, width, height, log);
    unsigned char **patterns = ScalarLBPRestrictedPatterns(sums, width, height, log);

    unsigned char **sat_patterns = (unsigned char**)vbx_shared_malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        sat_patterns[l] = (unsigned char*)vbx_shared_malloc(width*height*sizeof(unsigned char));
    }

    unsigned int *iImg, *iiImg;
    iImg = (unsigned int *)vbx_shared_malloc(width*height*sizeof(unsigned int));
    iiImg = (unsigned int *)vbx_shared_malloc(width*height*sizeof(unsigned int));
    gen_integrals(vbx_img, iImg, iiImg, width, height);

    image_t lbp_img = {iImg, {width, height}};
    for (l=0; l<log+1; l++) {
        cell = 1 << l;
        lbp_feat_t lbp_feat = {{{0, 0}, {cell, cell}}, 0, 0, {0, 0, 0, 0, 0, 0, 0, 0}};
        for (j = 0; j < height - (3*cell-1); j++) {
            for (i = 0; i < width - (3*cell-1); i++) {
                pair_t lbp_p = {i, j};
                sat_patterns[l][j*width+i] = SATBinaryPattern(lbp_img, &lbp_feat, lbp_p);
            }
        }
    }

    /* test patterns vs sat binary patterns */
    for (l=0; l<log+1; l++) {
        cell = 1 << l;
        for (j = 0; j < height - (3*cell-1); j++) {
            errors += match_array_byte(patterns[l] + j*width, sat_patterns[l] + j*width, "patterns", width - (3*cell-1), 1, 0, max_print_errors, 1, j);
            if (errors > max_print_errors){
                max_print_errors = 0;
            }
        }
    }
    return errors;
}

int compare_SATBinaryPattern_to_test_scalar_patterns(unsigned short *vbx_img, int log, int width, int height, int max_print_errors)
{
    int l, i, j, cell, errors = 0;

    unsigned char** scalar_patterns = test_scalar_patterns(vbx_img, log, width, height);

    unsigned char **sat_patterns = (unsigned char**)vbx_shared_malloc((log+1)*sizeof(unsigned char*));
    for (l=0; l<log+1; l++) {
        sat_patterns[l] = (unsigned char*)vbx_shared_malloc(width*height*sizeof(unsigned char));
    }

    unsigned int *iImg, *iiImg;
    iImg = (unsigned int *)vbx_shared_malloc(width*height*sizeof(unsigned int));
    iiImg = (unsigned int *)vbx_shared_malloc(width*height*sizeof(unsigned int));
    gen_integrals(vbx_img, iImg, iiImg, width, height);

    image_t lbp_img = {iImg, {width, height}};
    for (l=0; l<log+1; l++) {
        cell = 1 << l;
        lbp_feat_t lbp_feat = {{{0, 0}, {cell, cell}}, 0, 0, {0, 0, 0, 0, 0, 0, 0, 0}};
        for (j = 0; j < height - (3*cell-1); j++) {
            for (i = 0; i < width - (3*cell-1); i++) {
                pair_t lbp_p = {i, j};
                sat_patterns[l][j*width+i] = SATBinaryPattern(lbp_img, &lbp_feat, lbp_p);
            }
        }
    }

    /* test patterns vs sat binary patterns */
    for (l=0; l<log+1; l++) {
        cell = 1 << l;
        for (j = 0; j < height - (3*cell-1); j++) {
            errors += match_array_byte(scalar_patterns[l] + j*width, sat_patterns[l] + j*width, "patterns", width - (3*cell-1), 1, 0, max_print_errors, 1, j);
            if (errors > max_print_errors){
                max_print_errors = 0;
            }
        }
    }
    return errors;
}

int compare_ScalarLBPRestrictedPatterns_to_test_scalar_patterns(unsigned short* vbx_img, int log, int width, int height, int max_print_errors)
{
    int i, j, cell, errors = 0;

    /* generate patterns */
    unsigned short **sums = ScalarLBPRestrictedSums(vbx_img, width, height, log);
    unsigned char **patterns = ScalarLBPRestrictedPatterns(sums, width, height, log);

    unsigned char** scalar_patterns = test_scalar_patterns(vbx_img, log, width, height);

    /* test patterns vs scalar patterns */
    for (j=0; j<log+1; j++) {
        cell = 1 << j;
        for (i = 0; i < height - (3*cell-1); i++) {
            errors += match_array_byte(scalar_patterns[j] + i*width, patterns[j] + i*width, "patterns", width - (3*cell-1), 1, 0, max_print_errors, 1, i);
            if (errors > max_print_errors){
                max_print_errors = 0;
            }
        }
    }
    return errors;
}
int compare_ScalarLBPRestrictedSums_to_test_scalar_sums_half(unsigned short* vbx_img, int log, int width, int height, int max_print_errors)
{
    int i, j, cell, errors = 0;

    /* generate sums */
    unsigned short **sums = ScalarLBPRestrictedSums(vbx_img, width, height, log);
    unsigned short **scalar_sums = test_scalar_sums_half(vbx_img, log, width, height);

    /* test sums vs scalar sums */
    for (j=0; j<log+1; j++) {
        cell = 1 << j;
        for (i = 0; i < height - (cell-1); i++) {
            errors += match_array_half(scalar_sums[j] + i*width, sums[j] + i*width, "sum", width - (cell-1), 1, 0, max_print_errors, i);
            if (errors > max_print_errors){
                max_print_errors = 0;
            }
        }
    }
    return errors;
}

int compare_LBPRestrictedSums2_to_test_scalar_sums_half(unsigned short* vbx_img, int log, int width, int height, int max_print_errors)
{
    int cell, errors = 0;
    int i, j;

    /* generate sums */
    unsigned short **sums = LBPRestrictedSums2(vbx_img, width, height, log);
    unsigned short **scalar_sums = test_scalar_sums_half(vbx_img, log, width, height);

    /* test sums vs scalar sums */
    for (j=1; j<log+1; j++) {
        cell = 1 << j;
        for (i = 0; i < height - (cell - 1); i++) {
            errors += match_array_half(scalar_sums[j]+i*width, sums[j]+i*width, "sum", width - (cell - 1), 1, 0, max_print_errors, i);
            if (errors > max_print_errors){
                max_print_errors = 0;
            }
        }
    }
    return errors;
}

int compare_LBPRestrictedSums_to_test_scalar_sums_byte(unsigned short* vbx_img, int log, int width, int height, int max_print_errors)
{
    int cell, errors = 0;
    int i, j;

    /* generate sums */
    unsigned char **sums = LBPRestrictedSums(vbx_img, width, height, log);
    unsigned char **scalar_sums = test_scalar_sums_byte(vbx_img, log, width, height);

    /* test sums vs scalar sums */
    for (j=0; j<log+1; j++) {
        cell = 1 << j;
        for (i = 0; i < height - (cell-1); i++) {
            errors += match_array_byte(scalar_sums[j]+i*width, sums[j]+i*width, "sum", width - (cell-1), 1, 0, max_print_errors, 0, i);
            if (errors > max_print_errors){
                max_print_errors = 0;
            }
        }
    }
    return errors;
}

unsigned char** test_scalar_patterns(unsigned short* vbx_img, int log, int width, int height)
{
    int k, j, i, n, m;
    unsigned char** scalar_patterns = (unsigned char**)malloc((log+1) * sizeof(unsigned char*));
    for (k = 0; k < log+1; k++){
        scalar_patterns[k] = (unsigned char*)malloc(width*height*sizeof(unsigned char));
        
    }
    int p, p0, px, cell, pattern;
    int coeff_x[] = {0, 1, 2, 2, 2, 1, 0, 0};
    int coeff_y[] = {0, 0, 0, 1, 2, 2, 2, 1};
    for (k = 0; k < log+1; k++){
        cell = 1 << k;
        for (j = 0; j < height - (3*cell-1); j++){
            for (i = 0; i < width - (3*cell-1); i++){
                /* calculate center */
                p0 = 0;
                for (n = 0; n < cell; n++) {
                    for (m = 0; m < cell; m++){
                        p0 += vbx_img[(j+n+cell)*width+(i+m+cell)];
                    }
                }
                pattern = 0;
                for (p = 0; p < 8; p++){
                    px = 0;
                    for (n = 0; n < cell; n++){
                        for (m = 0; m < cell; m++){
                            px += vbx_img[(j+n+coeff_y[p]*cell)*width+(i+m+coeff_x[p]*cell)];
                        }
                    }
                    if(px >= p0){
                        pattern += 1 << (7-p);
                    }
                }
                scalar_patterns[k][j*width+i] = pattern;
            }
        }
    }
    return scalar_patterns;
}

unsigned char** test_scalar_sums_byte(unsigned short* vbx_img, int log, int width, int height)
{
    int k, j, i, n, m, cell;
    unsigned char** scalar_sums = (unsigned char**)malloc((log+1) * sizeof(unsigned char*));
    for (k = 0; k < log+1; k++){
        scalar_sums[k] = (unsigned char*)vbx_shared_malloc(width*height*sizeof(unsigned char));
    }
    for (k = 0; k < log+1; k++){
        cell = 1 << k;
        for (j = 0; j < height - (cell-1); j++){
            for (i = 0; i < width - (cell-1); i++){
                unsigned total = 0;
                for (n = 0; n < cell; n++){
                    for (m = 0; m < cell; m++){
                        total += vbx_img[(j+n)*width+(i+m)];
                    }
                }
                scalar_sums[k][j*width+i] = total >> (2*k);
            }
        }
    }
    return scalar_sums;
}

unsigned short** test_scalar_sums_half(unsigned short* vbx_img, int log, int width, int height)
{
    int k, j, i, n, m, cell;
    unsigned short** scalar_sums = (unsigned short**)malloc((log+1) * sizeof(unsigned short*));
    for (k = 0; k < log+1; k++){
        scalar_sums[k] = (unsigned short*)vbx_shared_malloc(width*height*sizeof(unsigned short));
    }
    for (k = 0; k < log+1; k++){
        cell = 1 << k;
        for (j = 0; j < height - (cell - 1); j++){
            for (i = 0; i < width - (cell - 1); i++){
                unsigned total = 0;
                for (n = 0; n < cell; n++){
                    for (m = 0; m < cell; m++){
                        total += vbx_img[(j+n)*width+(i+m)];
                    }
                }
                scalar_sums[k][j*width+i] = total;
            }
        }
    }
    return scalar_sums;
}

int compare_scalar_rgb2luma_to_vbw_rgb2luma16(unsigned short *img, unsigned short *vbx_img, pixel *vbx_input, int width, int height, int stride, int max_print_errors)
{
    int errors;
    scalar_rgb2luma(img, vbx_input, width, height, stride);
    vbw_rgb2luma16(vbx_img, vbx_input, width, height, stride);
    vbx_sync();

    errors = match_array_half(img, vbx_img, "greyscale", width, height, 0, max_print_errors, 0);
    return errors;
}

int test_lbp_ci(unsigned short* img, int width, int height)
{

    vbx_uhalf_t* v_a1  = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_b1  = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_1h = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));

    vbx_uhalf_t* v_a2  = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_b2  = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_2h  = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));

    vbx_uhalf_t* v_a4  = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_b4  = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));
    vbx_uhalf_t* v_4h = (vbx_uhalf_t*)vbx_sp_malloc(width*sizeof(vbx_uhalf_t));

    vbx_ubyte_t* v_1b  = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_2b  = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));
    vbx_ubyte_t* v_4b  = (vbx_ubyte_t*)vbx_sp_malloc(width*sizeof(vbx_ubyte_t));

    unsigned short* lbp1h = (unsigned short*)vbx_shared_malloc(width*sizeof(unsigned short));
    unsigned short* lbp2h = (unsigned short*)vbx_shared_malloc(width*sizeof(unsigned short));
    unsigned short* lbp4h = (unsigned short*)vbx_shared_malloc(width*sizeof(unsigned short));

    unsigned char* lbp1b = (unsigned char*)vbx_shared_malloc(width*sizeof(unsigned char));
    unsigned char* lbp2b = (unsigned char*)vbx_shared_malloc(width*sizeof(unsigned char));
    unsigned char* lbp4b = (unsigned char*)vbx_shared_malloc(width*sizeof(unsigned char));

    img = img + width;

    vbx_dma_to_vector(v_a1, img,         width*sizeof(unsigned short));
    vbx_dma_to_vector(v_b1, img + width, width*sizeof(unsigned short));
    vbx_dma_to_vector(v_a2, img,         width*sizeof(unsigned short));
    vbx_dma_to_vector(v_b2, img + width, width*sizeof(unsigned short));
    vbx_dma_to_vector(v_a4, img,         width*sizeof(unsigned short));
    vbx_dma_to_vector(v_b4, img + width, width*sizeof(unsigned short));
    vbx_sync();

    int i;
    int m = 48;
    for(i=0; i<m; i++){
        v_a1[i] = 0;
        v_b1[i] = 0;
        v_a2[i] = 0;
        v_b2[i] = 0;
        v_a4[i] = 0;
        v_b4[i] = 0;
    }
    int n = 12;
    int src_a1[] = {0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int src_b1[] = {0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    int src_a2[] = {0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int src_b2[] = {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    int src_a4[] = {0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0};
    int src_b4[] = {0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0};
    
    for(i=0; i<16; i++){
        v_a1[i] = src_a1[i];
        v_b1[i] = src_b1[i];
        v_a2[i] = src_a2[i];
        v_b2[i] = src_b2[i];
        v_a4[i] = src_a4[i];
        v_b4[i] = src_b4[i];
    }

    vbx_set_vl(width);
    vbx(VVHU, VCUSTOM1, v_1h, v_a1, v_b1); 
    vbx(VVHU, VCUSTOM2, v_2h, v_a2, v_b2); 
    vbx(VVHU, VCUSTOM3, v_4h, v_a4, v_b4); 
    vbx(VVHB, VADD, v_1b, v_1h, ((vbx_byte_t*)v_1h) + 1);
    vbx(VVHB, VADD, v_2b, v_2h, ((vbx_byte_t*)v_2h) + 1);
    vbx(VVHB, VADD, v_4b, v_4h, ((vbx_byte_t*)v_4h) + 1);
    vbx_dma_to_host(lbp1h, v_1h, width*sizeof(unsigned short));
    vbx_dma_to_host(lbp2h, v_2h, width*sizeof(unsigned short));
    vbx_dma_to_host(lbp4h, v_4h, width*sizeof(unsigned short));
    vbx_dma_to_host(lbp1b, v_1b, width*sizeof(unsigned char));
    vbx_dma_to_host(lbp2b, v_2b, width*sizeof(unsigned char));
    vbx_dma_to_host(lbp4b, v_4b, width*sizeof(unsigned char));
    vbx_sync();

    test_print_array_half(v_a1, n);
    test_print_array_half(v_b1, n);
    test_print_hex_array_half(lbp1h, n);
    test_print_hex_array_byte(lbp1b, n);

    test_print_array_half(v_a2, n);
    test_print_array_half(v_b2, n);
    test_print_hex_array_half(lbp2h, n);
    test_print_hex_array_byte(lbp2b, n);

    test_print_array_half(v_a4, n);
    test_print_array_half(v_b4, n);
    test_print_hex_array_half(lbp4h, n);
    test_print_hex_array_byte(lbp4b, n);

    vbx_sp_free();
    vbx_shared_free(lbp1h);
    vbx_shared_free(lbp2h);
    vbx_shared_free(lbp4h);
    vbx_shared_free(lbp1b);
    vbx_shared_free(lbp2b);
    vbx_shared_free(lbp4b);
    return 0;
}

int test_lbp_lut_ci(unsigned short* img, int log, int width, int height)
{
    assign_lbp_lut_ci();
    LBPRestrictedCI(img, width, height, log);

    return 0;
}

int compare_vbx_lbp_ci_to_scalar_patterns(unsigned short* img, int log, int width, int height, int max_print_errors)
{
    int j, l, cell, max_cell, errors = 0;
    unsigned char** scalar_patterns = test_scalar_patterns(img, log, width, height);

    max_cell = 1<<log;
    vbx_uhalf_t* v_in = (vbx_uhalf_t*)vbx_sp_malloc((1+2*max_cell)*width*sizeof(vbx_half_t));
    vbx_uhalf_t* v_top = (vbx_half_t*)vbx_sp_malloc(width*sizeof(vbx_half_t));
    vbx_uhalf_t* v_bot = (vbx_half_t*)vbx_sp_malloc(width*sizeof(vbx_half_t));
    vbx_ubyte_t* v_lbp = (vbx_ubyte_t*)v_bot;

    unsigned char* lbp = (unsigned char*)vbx_shared_malloc(width*sizeof(unsigned char));

    vbx_set_vl(width);
    for(l = 0; l < 1; l++){
        cell = 1<<l;
        for(j=0; j < height - 2*cell; j++){
            vbx_dma_to_vector(v_in, img+j*width, (1+2*cell)*width*sizeof(unsigned short));
            vbx(VVHU, VCUSTOM1, v_top, v_in, v_in+(1*cell)*width); 
            vbx(VVHU, VCUSTOM1, v_bot, v_in+(1*cell)*width, v_in+(2*cell)*width); 
            vbx(SVHBU, VAND, (vbx_ubyte_t*)v_top, 0xf0, v_top);
            vbx(SVHBU, VAND, (vbx_ubyte_t*)v_bot, 0x0f, v_bot);
            vbx(VVBU, VADD, v_lbp, v_bot, v_top); 
            vbx_dma_to_host(lbp, v_lbp, width*sizeof(unsigned char));
            vbx_sync();

            errors += match_array_byte(lbp, scalar_patterns[l]+j*width, "custom_lbp", width-2*cell, 1, 0, max_print_errors, 1, j);
            if (errors > max_print_errors){
                max_print_errors = 0;
            }

        }
    }
    vbx_sp_free();
    vbx_shared_free(lbp);
    return errors;
}

int compare_vbx_lut_to_vbx_lut_ci(int stage, int max_print_errors)
{
	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
	int vci_lanes = this_mxp->vcustom0_lanes;
    int sz = this_mxp->scratchpad_size/(16*sizeof(vbx_ubyte_t));

    vbx_byte_t* v_pass = (vbx_byte_t*)vbx_sp_malloc(sz*sizeof(vbx_byte_t));
    vbx_ubyte_t* v_pattern = (vbx_ubyte_t*)vbx_sp_malloc(sz*sizeof(vbx_byte_t));
    vbx_ubyte_t* v_lutc = (vbx_ubyte_t*)vbx_sp_malloc(sz*sizeof(vbx_byte_t));
    vbx_ubyte_t* v_group = (vbx_ubyte_t*)vbx_sp_malloc(sz*sizeof(vbx_byte_t));
    vbx_ubyte_t* v_sel = (vbx_ubyte_t*)vbx_sp_malloc(sz*sizeof(vbx_byte_t));
    vbx_ubyte_t* v_lut = (vbx_ubyte_t*)vbx_sp_malloc(sz*sizeof(vbx_word_t));
    vbx_ubyte_t* v_idx = (vbx_ubyte_t*)vbx_sp_malloc(sz*sizeof(vbx_word_t));
    if(v_idx == NULL) {
        printf("failed to allocate in compare_vbx_lut_to_vbx_lut_ci\n");
    }

    unsigned char* lut = (unsigned char*)vbx_shared_malloc(sz*sizeof(unsigned char));
    unsigned char* lut_c = (unsigned char*)vbx_shared_malloc(sz*sizeof(unsigned char));

    int f, n, s, errors = 0;
    for (n = 0; n < sz; n++) {
        v_pattern[n] = (n & 0xff);
    }

    for (f = 0; f < face_lbp[stage].count; f++) {
        lbp_feat_t feat = face_lbp[stage].feats[f];

        vbx_set_vl(sz);
        int total = f;
        s = 0;
        while(s < stage){
            total += face_lbp[s].count;
            s++;
        }

        if(total < 256) {
            vbx(SVBU, VLBPLUT, v_lutc, total, v_pattern);
        } else {
            vbx(SVBS, VLBPLUT, v_lutc, total-256, v_pattern);
        }

        vbx(SVB, VMOV, v_pass, feat.fail, 0);
        /* check if pattern is in lut */
        vbx(SVBU, VSHR, v_group, 5, v_pattern);
        for (n = 0; n < 8; n++) {
            vbx(SVB, VADD, v_sel, -n, v_group);
            vbx(SVBW, VCMV_Z, v_lut, feat.lut[n], v_sel);
        }

        vbx(SVBWU, VAND, v_idx, 0x1f, v_pattern);
        vbx(VVWB, VSHR, v_lut, v_idx, v_lut);
        vbx(SVB, VAND, v_lut, 1, v_lut);
        vbx(SVB, VCMV_LEZ, v_pass, feat.pass, v_lut);

        vbx_dma_to_host(lut_c, v_lutc, sz*sizeof(unsigned char));
        vbx_dma_to_host(lut, v_pass, sz*sizeof(unsigned char));
        vbx_sync();

        errors += match_array_byte(lut, lut_c, "custom_lut", sz, 1, 0, max_print_errors, 0, 0);

    }
    vbx_sp_free();
    vbx_shared_free(lut);
    vbx_shared_free(lut_c);
    return errors;
}

int main(void)
{

	vbx_timestamp_t time_start, time_stop;
	double scalar_time, vbx_time, vbx_time_masked;
	int i, j, k, l, m, n;
	int errors = 0;

	vbx_test_init();
	vbx_mxp_print_params();
    pixel *input, *scalar_input, *vbx_input, *vbx_input_masked;
    uint16_t *scalar_short;

	input         = (pixel *)vbx_shared_malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(pixel));
	scalar_input  = (pixel *)vbx_remap_cached(input, IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(pixel));
	scalar_short  = (uint16_t *)malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(uint16_t));
	vbx_input    = (pixel *)vbx_shared_malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(pixel));
	vbx_input_masked  = (pixel *)vbx_shared_malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(pixel));

#if UNIT
    unsigned char *vbx_img8;
    unsigned short *img, *vbx_img;
    unsigned int *iImg, *vbx_iImg;
    unsigned int *iiImg, *vbx_iiImg;
    img = (unsigned short*)malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(unsigned short));
    vbx_img = (unsigned short*)vbx_shared_malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(unsigned short));
    vbx_img8 = (unsigned char*)vbx_shared_malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(unsigned char));

    iImg = (unsigned int*)malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(unsigned int));
    vbx_iImg = (unsigned int*)vbx_shared_malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(unsigned int));

    iiImg = (unsigned int*)malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(unsigned int));
    vbx_iiImg = (unsigned int*)vbx_shared_malloc(IMAGE_WIDTH*IMAGE_HEIGHT*sizeof(unsigned int));
#endif//UNIT

	printf("Resolution = %dx%d\n", IMAGE_WIDTH, IMAGE_HEIGHT);
    printf("Initializing data\n");
	vbx_timestamp_start();
    for(l = 0; l < 1; l++){
        char *src;
        char *sdst;
        char *vdst;
        char *mdst;
        if(l == 0){
            load_lenna(input, IMAGE_WIDTH, IMAGE_HEIGHT);
            load_lenna(vbx_input, IMAGE_WIDTH, IMAGE_HEIGHT);
            load_lenna(vbx_input_masked, IMAGE_WIDTH, IMAGE_HEIGHT);
            printf("\nLenna\n");
            src = "lenna";
            sdst = "s_lenna";
            vdst = "v_lenna";
            mdst = "m_lenna";
        }else if(l == 1){
            load_ms(input, IMAGE_WIDTH, IMAGE_HEIGHT);
            load_ms(vbx_input, IMAGE_WIDTH, IMAGE_HEIGHT);
            load_ms(vbx_input_masked, IMAGE_WIDTH, IMAGE_HEIGHT);
            printf("\nMicrosoft\n");
            src = "ms";
            sdst = "s_ms";
            vdst = "v_ms";
            mdst = "m_ms";
        }else if(l == 2){
            load_blank(input, IMAGE_WIDTH, IMAGE_HEIGHT);
            load_blank(vbx_input, IMAGE_WIDTH, IMAGE_HEIGHT);
            load_blank(vbx_input_masked, IMAGE_WIDTH, IMAGE_HEIGHT);
            printf("\nblank\n");
            src = "blank";
            sdst = "s_blank";
            vdst = "v_blank";
            mdst = "m_blank";
        }
#if UNIT
    int window = 20;
    int log=0;
    while(((window/3)>>log) >= 2) log++;


    errors += compare_scalar_rgb2luma_to_vbw_rgb2luma16(img, vbx_img, vbx_input, IMAGE_WIDTH, IMAGE_HEIGHT, IMAGE_WIDTH, MAX_PRINT_ERRORS);
    vbw_rgb2luma8(vbx_img8, vbx_input, IMAGE_WIDTH, IMAGE_HEIGHT, IMAGE_WIDTH);


    int s;
#if LUT_CI
#if DOUBLE_LUT
    printf("Testing double lut\n");

    printf("Assign lbp double lut\n");
    assign_lbp_lut_ci2();
    int prev = errors;
    printf("Cascade check\n");
    /* errors += cascade_check_2w(face_lbp, face_lbp_max_stage, 256); */
    /* errors += cascade_check_2h(face_lbp, face_lbp_max_stage, 256); */
    errors += cascade_check_2b(face_lbp, face_lbp_max_stage, 256);
    if (errors) {
        printf("errors %d\n", errors-prev);
    }
#else
    assign_lbp_lut_ci();

    printf("Testing cascade\n");

    int prev = errors;

    printf("lut check\n");

#if 0
#if 0
    errors += lut_check(256, 0, 0, 0);
    if (errors) {
        printf("errors %d\n", errors-prev);
    }
#elif 1

    int print_errors = 0;

	vbx_mxp_t *this_mxp = VBX_GET_THIS_MXP();
	int vci_lanes = this_mxp->vcustom0_lanes;
    int num_features = cascade_max_feature();
    int input_length = 10;
    int lut_length = num_features*vci_lanes;
    int lut_iterations = 15;
#if 1
    lut_length = input_length = 128;
    lut_iterations = 13;
    print_errors = 0;
    errors += lut_check2(input_length, lut_length, lut_iterations, print_errors);
    if (errors) {
        printf("errors %d\n", errors-prev);
    }
#elif 1
    input_length = 64;
    lut_length = input_length;
    lut_iterations = 13;
    print_errors = 1;
    errors += lut_check2(input_length, lut_length, lut_iterations, print_errors);
    if (errors) {
        printf("errors %d\n", errors-prev);
    }
#else
    for(s = 2; s < 100; s=s+10){
        errors += lut_check2(s, lut_length, lut_iterations, print_errors);
        if (errors - prev > 0) {
            printf("%d\terrors %d\n", s, errors-prev);
        } else {
            printf("%d\n", s);
        }
        prev = errors;
    }
#endif
#else
    for(s = 0; s < 2000; s=s+100){
        errors += lut_check(s, 0, 0, 0);
        if (errors - prev > 0) {
            printf("%d\terrors %d\n", s, errors-prev);
        } else {
            printf("%d\n", s);
        }
        prev = errors;
    }
#endif

#elif 1

#else
    printf("check cascade\n");
    prev = errors;
    errors += cascade_check(face_lbp, face_lbp_max_stage, 256);
    if (errors) {
        printf("errors %d\n", errors-prev);
    }

    printf("Testing LBP LUT CI\n");
    prev = errors;
    for(s = 0; s < face_lbp_max_stage; s++){
        errors += compare_vbx_lut_to_vbx_lut_ci(s, MAX_PRINT_ERRORS);
    }
    if (errors) {
        printf("errors %d\n", errors-prev);
        prev = errors;
    }
#endif
#endif
#endif

#if 0
    printf("Printing grey scale img\n");
    printf("grey = [");
    for (j = 0; j < IMAGE_HEIGHT; j++) {
        printf("[");
        for (i = 0; i < IMAGE_WIDTH; i++) {
            printf("%d, ", vbx_img8[j*IMAGE_WIDTH+i]);
        }
        printf("],\n");
    }
    printf("]\n");
#endif
#if LBP_CI
    printf("Testing LBP Pattern CI\n");
    errors += compare_LBPRestrictedCI_to_test_scalar_patterns(vbx_img, vbx_img8, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);
#endif

#if BLIP
    printf("Testing BLIP\n");
    for(s = 1; s < 10; s++){
        errors += compare_scalar_BLIP2_to_vector_BLIP(img, vbx_input, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS, s);
    }
#endif
#if 0
    errors += compare_LBPRestrictedSums_to_test_scalar_sums_byte(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);
    errors += compare_LBPRestrictedSums2_to_test_scalar_sums_half(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);
    errors += compare_ScalarLBPRestrictedSums_to_test_scalar_sums_half(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);
    errors += compare_ScalarLBPRestrictedPatterns_to_test_scalar_patterns(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);
    errors += compare_LBPRestrictedPatterns2_to_test_scalar_patterns(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);
    errors += compare_LBPRestricted_to_test_scalar_patterns(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);
    /* overflow issues -- using bytes changes lbp pattern */
    errors += compare_LBPRestrictedPatterns_to_test_scalar_patterns(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);

    /* requires SKIP_INTEGRALS 0 */
    errors += compare_gen_integrals_to_vector_get_img(img, iImg, iiImg, vbx_img, vbx_iImg, vbx_iiImg, vbx_input, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);


    /* redundant test, compare to test_scalar_patterns instead */
    errors += compare_ScalarLBPRestrictedPatterns_to_SATBinaryPattern(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);

    errors += compare_SATBinaryPattern_to_test_scalar_patterns(vbx_img, log, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);

    errors += compare_LBPPassStage_to_restricted(vbx_img, log, face_lbp[0], window, IMAGE_WIDTH, IMAGE_HEIGHT, MAX_PRINT_ERRORS);
#endif
#else // UNIT

#if PRINT
        print_python_pixel(scalar_input, src, IMAGE_WIDTH, IMAGE_HEIGHT);
#endif

        time_start = vbx_timestamp();
        scalar_rgb2luma(scalar_short, input, IMAGE_WIDTH, IMAGE_HEIGHT, IMAGE_WIDTH);
        scalar_face_detect_luma(scalar_short, input, IMAGE_WIDTH, IMAGE_HEIGHT, IMAGE_WIDTH, sdst);
        time_stop = vbx_timestamp();
        scalar_time = vbx_print_scalar_time(time_start, time_stop);
#if PRINT
        print_python_pixel(scalar_input, sdst, IMAGE_WIDTH, IMAGE_HEIGHT);
#endif
        printf("\nVector");
        time_start = vbx_timestamp();
        vector_face_detect((pixel *)vbx_input, IMAGE_WIDTH, IMAGE_HEIGHT, IMAGE_WIDTH, 0, vdst);
        time_stop = vbx_timestamp();
        vbx_time = vbx_print_vector_time(time_start, time_stop, scalar_time);
#if PRINT
        print_python_pixel(vbx_input, vdst, IMAGE_WIDTH, IMAGE_HEIGHT);
#endif

        printf("\nVector Masked");
        time_start = vbx_timestamp();
        vector_face_detect((pixel *)vbx_input_masked, IMAGE_WIDTH, IMAGE_HEIGHT, IMAGE_WIDTH, 1, mdst);
        time_stop = vbx_timestamp();
        vbx_time_masked = vbx_print_vector_time(time_start, time_stop, scalar_time);
#if PRINT
        print_python_pixel(vbx_input_masked, mdst, IMAGE_WIDTH, IMAGE_HEIGHT);
#endif
        /* errors += match_array_pixel(input, vbx_input, "vector", IMAGE_WIDTH, IMAGE_HEIGHT, 0, MAX_PRINT_ERRORS, 0); */
        /* errors += match_array_pixel(input, vbx_input_masked, "masked", IMAGE_WIDTH, IMAGE_HEIGHT, 0, MAX_PRINT_ERRORS, 0); */
        errors += match_array_pixel(vbx_input, vbx_input_masked, "masked", IMAGE_WIDTH, IMAGE_HEIGHT, 0, MAX_PRINT_ERRORS, 0);
#endif // UNIT
    }
	VBX_TEST_END(errors);
	return errors;
}
